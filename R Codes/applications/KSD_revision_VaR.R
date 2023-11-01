library(ccgarch) 

library(mvnTest) # Test for normality
library(MVN) # Test for normality
library(mvtnorm)
library(Matrix)

library(progress)
library(pryr)
library(compositions)
# library(MVN)  
# library(energy)
library(foreach)
library(doSNOW)

library(ccgarch)
library(DEoptim)
library(fMultivar) # rcauchy2d
library(statmod) # gauss.quad() used

library(sn) # skew normal

# library(plotly) # To create interactive charts
# library(timetk) # To manipulate the data series
# library(tibble)
################### Functions ################### 
# Tools
getDim <- function(x){
  if(is.array(x)){
    n <- dim(x)[1]; dimen <- dim(x)[2]
  }else{
    x <- array(x)
    n <- dim(x)[1]; dimen <- 1
  }
  result <- list("n" = n, "dim" = dimen)
  return(result)
}
find_median_distance <- function(Z){
  
  if(is.data.frame(Z)){
    Z = data.matrix(Z)
  }else{
    Z = as.array(Z)
  }
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]
  
  # if size of Z is greater than 100, randomly sample 100 points
  if(size1 > 100){
    if(is.na(size2)){
      Zmed <- Z[sample(size1,100)]
    }else{
      Zmed <- Z[sample(size1,100),]
    }
    size1 = 100
  }else{
    Zmed <- Z
  }
  
  Zmedsq <- Zmed * Zmed;
  if(is.na(dim(Z)[2]))
    G <- Zmedsq
  else
    G <- rowSums(Zmedsq)
  
  # Create row/col repeated matrices
  Q <- rep.col(G,size1)
  R <- rep.row(t(G),size1)
  
  dists <- Q + R - 2 * Zmed %*% t(Zmed)
  dists[lower.tri(dists, diag = TRUE)] = 0
  dists <- array(dists,dim=c(size1^2,1))
  median_dist <- median(dists[dists > 0 ])
  
  return(median_dist)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
repmat <- function(X,a=1,b=1){
  rows <- dim(X)[1]
  cols <- dim(X)[2]
  if(is.null(cols))
    cols <- 1
  rowRep <- matrix(rep(t(X),a),ncol = cols, byrow = TRUE)
  newX <- matrix(rep(rowRep, b),ncol=cols*b)
  return(newX)
}

msqrt <- function(M){
  # computes the square-root of a positive definite matrix
  if(!is.matrix(M))M=as.matrix(M)
  n1=nrow(M)
  if(n1 == 1){
    Mh=sqrt(M)
    Mhinv=1/Mh
  }
  if(n1 > 1){
    M=(M+t(M))/2
    m1=eigen(M)
    V=m1$vectors
    eiv=sqrt(m1$values)
    L=diag(eiv)
    Linv=diag(1/eiv)
    Mh=V%*%L%*%t(V)
    Mhinv=V%*%Linv%*%t(V)
  }
  msqrt <- list(mtxsqrt=Mh,invsqrt=Mhinv)
}
# KSD_new
# Input: 
#   x: Standardized Residual Matrix (n x d)
#   Y: Original Data Matrix (d x n)
#   theta: theta in H0 (1 x p)  [use invsqrt of sigma]
#   theta_hat: estimated theta (1 x p)
#   score_function: "gaussian"
#   width: -1
# Output: Sn (New KSD Statistic after adjusted for estimation effect)
KSD_new = function(x, score_function = "gaussian", width = -1){
  h = width
  # Set up score function
  if(typeof(score_function) == 'character'){
    if(tolower(score_function) == 'gaussian'){
      score_function <- function(val){
        return (-1 * val)
      }
    }
  }else{
    
  }
  # Set up the bandwidth of RBF Kernel
  if((typeof(h)== 'character' & tolower(h) == 'median') | h == -1){
    h = sqrt(0.5 * find_median_distance(x))
  }else if(h == -2){
    #h = 1
    medInfo = ratio_median_heuristic(x, score_function)
    h = medInfo$h
  }
  # Get parameters
  dimensions <-getDim(x)
  n <- dimensions$n; dimen <- dimensions$dim
  if(is.numeric(score_function) | is.vector(score_function)){
    Sqx = score_function
  }else{
    Sqx = score_function(x)
  }
  ############### S0 ###############
  XY = x %*% t(x)
  if(is.null(dim(x)[2])){
    x2 = x^2
    sumx = x * Sqx
  } else{
    x2 = array(rowSums(x^2),dim=c(n,1))
    sumx = rowSums(x * Sqx)
  }
  X2e = repmat(x2,1,n)
  H = (X2e + t(X2e) - 2*XY)
  Kxy = exp(-H/(2*h^2))
  sqxdy = -(Sqx %*% t(x) - repmat(sumx,1,n))/h^2
  dxsqy = t(sqxdy)
  dxdy = (-H/h^4 + dimen/h^2)
  # n x n matrix with (i,j)-th entry being u_q(xi, xj)
  M = (Sqx %*% t(Sqx) + sqxdy + dxsqy + dxdy) * Kxy
  M_d <- M - diag(diag(M))
  S0 <- sum(M_d)/(n*(n-1)) 
  
  return(S0)
}
KSD <- function(x, score_function, kernel='rbf', width=-1, nboot=1000){
  #cleanup()
  
  ######## NEED TO IMPLEMENT#######
  #list[kernel, h, nboot] = process_varargin(varargin, 'kernel', 'rbf', 'width', -1, 'nboot', 0)
  
  h <- width
  
  # If score_function = 'gaussian', calculate the score function for standard norm
  if(typeof(score_function) == 'character'){
    if(tolower(score_function) == 'gaussian'){
      score_function <- function(val){
        return (-1 * val)
      }
    }
  }else{
    
  }
  
  # Set up the bandwidth of RBF Kernel
  if((typeof(h)== 'character' & tolower(h) == 'median') | h == -1){
    h = sqrt(0.5 * find_median_distance(x))
  }else if(h == -2){
    #h = 1
    medInfo = ratio_median_heuristic(x, score_function)
    h = medInfo$h
  }
  
  # Get parameters
  dimensions <-getDim(x)
  n <- dimensions$n; dimen <- dimensions$dim
  if(is.numeric(score_function) | is.vector(score_function)){
    Sqx = score_function
  }else{
    Sqx = score_function(x)
  }
  
  
  if(tolower(kernel) == 'rbf'){
    XY <- x %*% t(x)
    if(is.null(dim(x)[2])){
      x2 <- x^2
      sumx <- x * Sqx
    } else{
      x2 <- array(rowSums(x^2),dim=c(n,1))
      sumx <- rowSums(x * Sqx)
    }
    X2e <- repmat(x2,1,n)
    
    H <- (X2e + t(X2e) - 2*XY)
    Kxy <- exp(-H/(2*h^2))
    
    sqxdy <- -(Sqx %*% t(x) - repmat(sumx,1,n))/h^2
    dxsqy <- t(sqxdy)
    dxdy <- (-H/h^4 + dimen/h^2)
    
    M = (Sqx %*% t(Sqx) + sqxdy + dxsqy + dxdy) * Kxy
    
  }else{
    warning('Wrong Kernel')
  }
  
  
  M2 <- M - diag(diag(M))
  ksd <- sum(M2)/(n*(n-1))
  ksdV <- sum(M) / (n^2)
  
  ##---------------------Bootstrap methods---------------------------##
  
  bootstrapSamples <-  rep(NA,nboot)
  
  ##Process arguments
  bootmethod <- 'weighted'
  switch(tolower(bootmethod),
         'weighted'={
           for(i in 1:nboot){
             wtsboot <- rmultinom(1,size=n,prob=rep(1/n,n))/n
             bootstrapSamples[i] <- (t(wtsboot)-1/n)%*%M2%*%(wtsboot-1/n)
           }
           p <- mean(bootstrapSamples >= ksd)
           sprintf("Weighted : %.5f",p)
         },
         'efron'={
           for(i in 1:nboot){
             xdx <- array(as.integer(runif(n=n,max=n))+1)
             bootstrapSamples[i] <- sum(M2[xdx,xdx]) / (n*(n-1))
           }
           p <- mean(bootstrapSamples <= 0)
         },{
           warning('Wrong bootmethod')
         }
  )
  
  info <- list("bandwidth"=h, "M" = M, "nboot" = nboot, "ksd_V" = ksdV)
  result <- list("ksd" = ksd, "p"=p, "bootStrapSamples"=bootstrapSamples, "info"=info)
  
  return(result)
}

d = 3

GARCH = function(dist){
  # Initialization
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = c(0, 0), sigma = diag(2), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = c(0,0), sigma = (nu-2)/nu * diag(2),  df = nu))
  }
  y = array(rep(0,n*2),dim = c(2,n))
  eps_t = c(0.5,0.5)
  sigma_ii = c(0.3, 0.3)
  # Iteration
  # library(MTS)
  for (t in 1:(n+ncut)){
    sigma_ii = a + B %*% eps_t^2 + G %*% sigma_ii
    sigma_21 = 0.5 * sqrt(sigma_ii[1]) * sqrt(sigma_ii[2])
    sigma_t = matrix(c(sigma_ii[1], sigma_21,
                       sigma_21, sigma_ii[2]),2,2)
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    if (t > ncut){y[,t-ncut] = eps_t}
  }
  # detach("package:MTS", unload=TRUE)
  
  return(y)
}
GARCH_d = function(dist){
  # library(MTS)
  # Initialization
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d),  df = nu))
  } else if (dist == 'N+t'){
    ita = t(((1 - lambda) * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol") +
               lambda * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)) / 
              sqrt(lambda^2 + (1 - lambda)^2))
  } else if (dist == 'N+logN'){
    temp = t(matrix(rlnorm.rplus(n + ncut, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n + ncut))
    temp = t(msqrt(var(t(temp)))$invsqrt %*% (temp - rowMeans(temp)))
    ita = t(((1 - lambda) * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol") +
               lambda * temp) / 
              sqrt(lambda^2 + (1 - lambda)^2))
  }
  y = array(rep(0,d*n),dim = c(d,n))
  eps_t = rep(0.5, d)
  sigma_ii = rep(0.3, d)
  C_t = R # Constant Correlation Matrix
  # Iteration
  for (t in 1:(n+ncut)){
    sigma_ii = a + B %*% eps_t^2 + G %*% sigma_ii
    D_t = diag(c(sqrt(sigma_ii)))
    sigma_t = D_t %*% C_t %*% D_t
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    if (t > ncut){y[,t-ncut] = eps_t}
  }
  # detach("package:MTS", unload=TRUE)
  
  return(y)
}
GARCH. = function(dist, a., B., G., R.){
  # library(MTS)
  # Initialization
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = c(0, 0), sigma = diag(2), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = c(0,0), sigma = (nu-2)/nu * diag(2),  df = nu))
  }
  y = array(rep(0,n*2),dim = c(2,n))
  eps_t = c(0.5,0.5)
  sigma_ii = c(0.3, 0.3)
  r = R.[1,2]
  # Iteration
  for (t in 1:(n+ncut)){
    sigma_ii = a. + B. %*% eps_t^2 + G. %*% sigma_ii
    sigma_21 = r * sqrt(sigma_ii[1]) * sqrt(sigma_ii[2])
    sigma_t = matrix(c(sigma_ii[1], sigma_21,
                       sigma_21, sigma_ii[2]),2,2)
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    if (t > ncut){y[,t-ncut] = eps_t}
  }
  # detach("package:MTS", unload=TRUE)
  return(y)
}
library(sn)
# gamma_hat = c(0.3, 0.1, 0.1)
# gamma_hat = c(-0.42, 0,-0.25)
CP = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat)
DP = cp2dp(CP, "SN")
# delta = c(0.2, 0.6, 0.4)
# alpha = delta / sqrt(1 - sum(delta^2))

# cp <- list(mean=c(0,0), var.cov=matrix(c(3,2,2,3)/3, 2, 2), gamma1=c(0.8, 0.4))
# dp <- cp2dp(cp, "SN")
# rnd <- rmsn(5, dp=dp)

# In use
GARCH_d. = function(dist, a., B., G., R.){
  # library(MTS)
  # Initialization
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d),  df = nu))
  } else if (dist == 'SN'){
    ita = t(rmsn(n = n + ncut, dp = DP))
    # ita = t(rmsn(n = n + ncut, xi = rep(0, d), Omega = diag(3), alpha = alpha))
  }
  y = array(rep(0,d*(n+ncut)),dim = c(d,n+ncut))
  eps_t = rep(0.5, d)
  sigma_ii = rep(0.3, d)
  C_t = R.
  # Iteration
  for (t in 1:(n+ncut)){
    sigma_ii = a. + B. %*% eps_t^2 + G. %*% sigma_ii
    D_t = diag(c(sqrt(sigma_ii)))
    sigma_t = D_t %*% C_t %*% D_t
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    # if (t > ncut){y[,t-ncut] = eps_t}
    y[,t] = eps_t
  }
  # detach("package:MTS", unload=TRUE)
  
  return(y)
}


# Specify initial parameters for estimation
a_init <<- rep(0.005, d)
B_init <<- diag(rep(0.2 - 0.01, d)) + matrix(0.01, d, d)
G_init <<- diag(rep(0.4 - 0.01, d)) + matrix(0.01, d, d)
R_init <<- diag(rep(1 - 0.1, d)) + matrix(0.1, d, d)

# Repitition Number
nrep = 1000
# Burn-in Obs Generating GARCH Process
ncut = 1000
# d.f. for t dist.
# nu = 8
level = 0.05
nB = 1000

# w = c(1/3,1/3,1/3)
# w = c(1.4, -0.2, -0.2)
# w = c(0.5,0.2,0.3)




test_size = 1200
######################################################


# Import data
data_AAPL = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD revision/AAPL.csv", header = TRUE)
data_CSCO = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD revision/CSCO.csv", header = TRUE)
data_INTC = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD revision/INTC.csv", header = TRUE)
data_MSFT = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD revision/MSFT.csv", header = TRUE)

AAPL = data_AAPL$Adj.Close
CSCO = data_CSCO$Adj.Close
INTC = data_INTC$Adj.Close
MSFT = data_MSFT$Adj.Close

# Daily log return in percentages
# Daily log return in percentages
lr_AAPL = (log(AAPL[2:2850]) - log(AAPL[1:2849])) * 100
lr_CSCO = (log(CSCO[2:2850]) - log(CSCO[1:2849])) * 100
lr_INTC = (log(INTC[2:2850]) - log(INTC[1:2849])) * 100
lr_MSFT = (log(MSFT[2:2850]) - log(MSFT[1:2849])) * 100

data_all = rbind(lr_MSFT, lr_AAPL, lr_INTC)

# lr_SP5t = (log(SP5t[2:7577]) - log(SP5t[1:7576])) * 100
# lr_CSCOt = (log(CSCOt[2:7577]) - log(CSCOt[1:7576])) * 100
# lr_INTCt = (log(INTCt[2:7577]) - log(INTCt[1:7576])) * 100

# data_all = rbind(lr_SP5t, lr_CSCOt, lr_INTCt)

# data_test = data_all[,2275:(2274+test_size)]

data_test = data_all[,(2849-test_size+1):2849]
data = data_all[,1:(2849-test_size)]

# mu_hat = rowMeans(data)
# sigma_hat = c(sd(lr_SP5),sd(lr_CSCO),sd(lr_INTC))
# rho_hat = cor(t(data))

# skewness(lr_SP5)
# ljung.box.test(lr_SP5)
# detach("package:MTS", unload=TRUE)
library(vars)

# Specify initial parameters for estimation
d = 3
a_init <<- rep(0.005, d)
B_init <<- diag(rep(0.2 - 0.01, d)) + matrix(0.01, d, d)
G_init <<- diag(rep(0.4 - 0.01, d)) + matrix(0.01, d, d)
R_init <<- diag(rep(1 - 0.1, d)) + matrix(0.1, d, d)

mu_list = NULL
sigma_list = list()

# data = data_all[,1:2274]

# fit1 = VAR(t(data), type = "const", lag.max = 5, ic = "SC")
fit1 = VAR(t(data), p = 3, type = "const")
fit2 = restrict(x = fit1, method = "ser", thresh = 1.5)
Yt = residuals(fit2)
order = fit2$p[[1]]
coef = Acoef(fit2)
A1 = coef[[1]]
A2 = coef[[2]]
A3 = coef[[3]]
M = Bcoef(fit2)[,10]

est_s1 = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                  R = R_init, dvar = Yt, model = "extended"),
                  error = function(e) {return(NA)},
                  warning = function(w) NA)
Yt.res = est_s1$std.resid
ht = est_s1$h
a = est_s1$para.mat$a
B = est_s1$para.mat$A
G = est_s1$para.mat$B
R = est_s1$para.mat$R

for (t in 1:test_size){
  mu_t_l = M + A1 %*% data_all[,(2849-test_size+t-1)] + 
    A2 %*% data_all[,(2849-test_size+t-2)] + 
    A3 %*% data_all[,(2849-test_size+t-3)]
  mu_list = cbind(mu_list, mu_t_l)
  
  sigma2_t_l = a + B %*% (Yt[2849-test_size-order,])^2 +
    G %*% ht[2849-test_size-order,]
  sigma_list[[t]] = diag(as.numeric(sqrt(sigma2_t_l))) %*% R %*% diag(as.numeric(sqrt(sigma2_t_l)))
}

############################################
data_w = t(data)
mean_ret <- mu_list[,1]
print(round(mean_ret, 5))

cov_mat <- sigma_list[[1]] * 252
print(round(cov_mat,4))

num_port <- 50000

# Creating a matrix to store the weights
all_wts <- matrix(nrow = num_port,
                  ncol = 3)

# Creating an empty vector to store Portfolio returns
port_returns <- vector('numeric', length = num_port)

# Creating an empty vector to store Portfolio Standard deviation
port_risk <- vector('numeric', length = num_port)

# Creating an empty vector to store Portfolio Sharpe Ratio
sharpe_ratio <- vector('numeric', length = num_port)

for (i in seq_along(port_returns)) {
  
  wts <- runif(3)
  wts <- wts/sum(wts)
  
  # Storing weight in the matrix
  all_wts[i,] <- wts
  
  # Portfolio returns
  port_ret <- sum(wts * mean_ret)
  port_ret <- ((port_ret + 1)^252) - 1
  
  # Storing Portfolio Returns values
  port_returns[i] <- port_ret
  
  # Creating and storing portfolio risk
  port_sd <- sqrt(t(wts) %*% (cov_mat  %*% wts))
  port_risk[i] <- port_sd
  
  # Creating and storing Portfolio Sharpe Ratios
  # Assuming 0% Risk free rate
  sr <- port_ret/port_sd
  sharpe_ratio[i] <- sr
  
}

# Storing the values in the table
# portfolio_values <- tibble(Return = port_returns,
#                            Risk = port_risk,
#                            SharpeRatio = sharpe_ratio)
# 
# head(portfolio_values)
portfolio_values <- cbind(port_returns, port_risk, sharpe_ratio)
# Converting matrix to a tibble and changing column names
# all_wts <- tk_tbl(all_wts)
portfolio_values <- cbind(all_wts, portfolio_values)
head(portfolio_values)

min_var <- portfolio_values[which.min(portfolio_values[,5]),]
max_sr <- portfolio_values[which.max(portfolio_values[,6]),]

# w = min_var[1:3]
w = max_sr[1:3]


########################################
# For normal distribution
z_norm_l = qnorm(level)
z_norm_s = qnorm(1-level)
VaR_list_N_l = NULL
VaR_list_N_s = NULL

for (t in 1:test_size){
  temp_l = t(w) %*% mu_list[,t] + z_norm_l * sqrt(t(w) %*% sigma_list[[t]] %*% w)
  temp_s = t(w) %*% mu_list[,t] + z_norm_s * sqrt(t(w) %*% sigma_list[[t]] %*% w)
  VaR_list_N_l = c(VaR_list_N_l, temp_l)
  VaR_list_N_s = c(VaR_list_N_s, temp_s)
}

# For t distribution
nu = 9
z_t_l = qt(level, df = nu) * sqrt((nu-2)/nu)
z_t_s = qt(1-level, df = nu) * sqrt((nu-2)/nu)
VaR_list_t_l = NULL
VaR_list_t_s = NULL

for (t in 1:test_size){
  temp_l = t(w) %*% mu_list[,t] + z_t_l * sqrt(t(w) %*% sigma_list[[t]] %*% w)
  temp_s = t(w) %*% mu_list[,t] + z_t_s * sqrt(t(w) %*% sigma_list[[t]] %*% w)
  VaR_list_t_l = c(VaR_list_t_l, temp_l)
  VaR_list_t_s = c(VaR_list_t_s, temp_s)
}

# For SN distribution
VaR_list_SN_l = NULL
VaR_list_SN_s = NULL
r_list = NULL
for (t in 1:test_size){
  z = rsn(n=1000, xi = 0, omega = 1, alpha = 0)
  temp = t(w) %*% mu_list[,t] + sqrt(t(w) %*% sigma_list[[t]] %*% w) * z
  r_list = c(temp)
  VaR_list_SN_l = c(VaR_list_SN_l,sort(r_list)[1000*(level)])
  VaR_list_SN_s = c(VaR_list_SN_s,sort(r_list)[1000*(1-level)])
}


# Actual Portforlio Return
rp = w %*% data_all[,(2849-test_size+1):(2849)]

# # failure rate
fr_N_l = rp < VaR_list_N_l
fr_N_s = rp > VaR_list_N_s

fr_t_l = rp < VaR_list_t_l
fr_t_s = rp > VaR_list_t_s

fr_SN_l = rp < VaR_list_SN_l
fr_SN_s = rp > VaR_list_SN_s

temp = function(fr){
  J = test_size
  N = sum(fr)
  f_hat = N / J
  return(f_hat)
}

temp(fr_N_l)
temp(fr_N_s)
temp(fr_t_l)
temp(fr_t_s)
temp(fr_SN_l)
temp(fr_SN_s)

# 
# # Kupiec LR Statistic
# LR = function(fr){
#   J = test_size
#   N = sum(fr)
#   f_hat = N / J
#   LR = 2 * log(f_hat^N * (1-f_hat)^(J-N)) - 2 * log(level^N * (1-level)^(J-N))
#   return(1-pchisq(LR, 1))
# }
# c(LR(fr_N_l),
# LR(fr_N_s),
# LR(fr_t_l),
# LR(fr_t_s),
# LR(fr_SN_l),
# LR(fr_SN_s))

library(rugarch)
round(c(VaRTest(alpha = level, rp, VaR_list_N_l, conf.level = 0.95)$cc.LRp,
        VaRTest(alpha = 1-level, rp, VaR_list_N_s, conf.level = 0.95)$cc.LRp,
        VaRTest(alpha = level, rp, VaR_list_t_l, conf.level = 0.95)$cc.LRp,
        VaRTest(alpha = 1-level, rp, VaR_list_t_s, conf.level = 0.95)$cc.LRp,
        VaRTest(alpha = level, rp, VaR_list_SN_l, conf.level = 0.95)$cc.LRp,
        VaRTest(alpha = 1-level, rp, VaR_list_SN_s, conf.level = 0.95)$cc.LRp),3)

# round(c(temp(fr_t_l),
#         temp(fr_t_s),
#         VaRTest(alpha = level, rp, VaR_list_t_l, conf.level = 0.95)$cc.LRp,
#         VaRTest(alpha = 1-level, rp, VaR_list_t_s, conf.level = 0.95)$cc.LRp),3)




############################################
round(w,3)

# test_size = 600
# 0.375 0.126 0.375 0.187 0.375 0.126
# test_size = 800
# 0.049 0.293 0.021 0.171 0.049 0.369
# test_size = 1200
# nu = 8
# 0.355 0.826 0.326 0.654 0.342 0.507
# nu = 5   0.218 0.067
# nu = 6   0.277 0.287
# nu = 7   0.326 0.654
# nu = 9   0.326 0.721
# test_size = 1500
# 0.177 0.834 0.078 0.480 0.041 0.972