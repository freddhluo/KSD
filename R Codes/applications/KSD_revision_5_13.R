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

nB = 1000

######################################################
# Import data
data_AAPL = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/AAPL.csv", header = TRUE)
data_CSCO = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/CSCO.csv", header = TRUE)
data_INTC = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/INTC.csv", header = TRUE)
data_MSFT = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/MSFT.csv", header = TRUE)

data_AA = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/AA.csv", header = TRUE)
data_CAT = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/CAT.csv", header = TRUE)
data_DIS = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/DIS.csv", header = TRUE)

AAPL = data_AAPL$Adj.Close
CSCO = data_CSCO$Adj.Close
INTC = data_INTC$Adj.Close
MSFT = data_MSFT$Adj.Close

AA = data_AA$Adj.Close
CAT = data_CAT$Adj.Close
DIS = data_DIS$Adj.Close
# Daily log return in percentages
lr_AAPL = (log(AAPL[2:2850]) - log(AAPL[1:2849])) * 100
lr_CSCO = (log(CSCO[2:2850]) - log(CSCO[1:2849])) * 100
lr_INTC = (log(INTC[2:2850]) - log(INTC[1:2849])) * 100
lr_MSFT = (log(MSFT[2:2850]) - log(MSFT[1:2849])) * 100

lr_AA = (log(AA[2:2850]) - log(AA[1:2849])) * 100
lr_CAT = (log(CAT[2:2850]) - log(CAT[1:2849])) * 100
lr_DIS = (log(DIS[2:2850]) - log(DIS[1:2849])) * 100
# data = rbind(lr_MSFT, lr_CSCO, lr_INTC)
# data = rbind(lr_MSFT, lr_AAPL, lr_INTC)
data = rbind(lr_AAPL, lr_CAT, lr_DIS)

data = data[,1:1849]

library(lubridate)

date = as.character(data_MSFT$Date[2:2850])
date = as.Date(date, "%Y-%m-%d")
plot.data = data.frame(date, lr_MSFT, lr_AAPL, lr_INTC)
library(scales)
par(mfcol = c(3,1))
p1 = ggplot(data = plot.data, aes(x = date, y = lr_MSFT)) +
  geom_line(color = "black", size = 0.5) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
  xlab('Year') + ylab('log returns') + ggtitle('(A) Microsoft Corporation')

p2 = ggplot(data = plot.data, aes(x = date, y = lr_AAPL)) +
  geom_line(color = "black", size = 0.5) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
  xlab('Year') + ylab('log returns') + ggtitle('(B) Apple Corporation')

p3 = ggplot(data = plot.data, aes(x = date, y = lr_INTC)) +
  geom_line(color = "black", size = 0.5) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
  xlab('Year') + ylab('log returns') + ggtitle('(C) Intel Corporation')
p1
p2
p3
library(cowplot)
p4 <- cowplot::plot_grid(p1, p2, p3, nrow = 3)#将p1-p4四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）
p4

mu_hat = rowMeans(data)
# sigma_hat = c(sd(lr_MSFT),sd(lr_AAPL),sd(lr_INTC))
sigma_hat = c(sd(lr_AAPL),sd(lr_CAT),sd(lr_DIS))

rho_hat = cor(t(data))

# skewness(lr_SP5)
# ljung.box.test(lr_SP5)
detach("package:MTS", unload=TRUE)
library(vars)
fit1 = VAR(t(data), p = 3, type = "const")
fit2 = restrict(x = fit1, method = "ser", thresh = 1.5)
Yt = residuals(fit2)
order = fit2$p[[1]]
coef = Acoef(fit2)

# Specify initial parameters for estimation
d = 3
a_init <<- rep(0.005, d)
B_init <<- diag(rep(0.2 - 0.01, d)) + matrix(0.01, d, d)
G_init <<- diag(rep(0.4 - 0.01, d)) + matrix(0.01, d, d)
R_init <<- diag(rep(1 - 0.1, d)) + matrix(0.1, d, d)

est_s1 = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                  R = R_init, dvar = Yt, model = "extended"),
                  error = function(e) {return(NA)},
                  warning = function(w) NA)


Yt.res = est_s1$std.resid

a = est_s1$para.mat$a
B = est_s1$para.mat$A
G = est_s1$para.mat$B
R = est_s1$para.mat$R

R.invsqrt = tryCatch(msqrt(R)$invsqrt,
                     warning = function(w) matrix(NA, d, d),
                     error = function(e) matrix(NA, d, d))

Xt = R.invsqrt %*% t(Yt.res)
n = ncol(Xt)
p = 3
# rowMeans(Xt)
# var(t(Xt))
# skewness(Xt[3,])

######### estimate t dist. d.f. parameter nu #########
# method 1: ML
logL = function(nu){
  sum = 0
  for (i in 1:n){
    sum = sum + log(nu + sum(Xt[,i]^2))
  }
  temp = n * log(gamma((nu + p) / 2)) - n * log(gamma(nu / 2)) +
    n * nu / 2 * log(nu) - (nu + p) / 2 * sum
  return(-temp)
}
opt = optim(par = 5, fn = logL, method = 'L-BFGS-B', lower = 2 + 0.00001, upper = 100)
nu = opt$par

# method 2 ECME
eqt = function(nu){
  sum = 0
  for (i in 1:n){
    delta_i = sum(Xt[,i]^2) * nu / (nu - 2)
    wi = (nu + p) / (nu + delta_i)
    sum = sum + log(wi) - wi
  }
  temp = - digamma(nu / 2) + log(nu / 2) + 1 / n * sum + 1 + digamma((nu + p) / 2) - log((nu + p) / 2)
  return(temp)
}
list = NULL
for (i in 3:30){
  list = rbind(list, c(i, eqt(i)))
}
plot(list[,1], list[,2])

library(pracma)
sol = bisect(eqt, 4, 100)
nu = sol$root

# eqt = function(nu){
#   sum = 0
#   for (i in 1:n){
#     delta_i = sum(Xt[,i]^2) * nu / (nu - 2)
#     wi = (nu + p) / (nu + delta_i)
#     sum = sum + log(wi) - wi
#   }
#   temp = - digamma(nu / 2) + log(nu / 2) + 1 / n * sum + 1 + digamma((nu + p) / 2) - log((nu + p) / 2)
#   return(temp)
# }
# list = NULL
# for (i in 3:30){
#   list = rbind(list, c(i, eqt(i)))
# }
# plot(list[,1], list[,2])
#
# library(pracma)
# sol = bisect(eqt, 4, 100)
# nu = sol$root

########## estimate skew normal dist. parameter gamma #########
# logL = function(gamma){
#   CP = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma)
#   DP = cp2dp(CP, "SN")
#   Omega = DP$Omega
#   xi = DP$beta
#   alpha = DP$alpha
#   c = sign(gamma_hat)*(2*abs(gamma_hat)/(4-pi))^(1/3)
#   mu_z = c / (sqrt(1+c^2))
#   Sigma_z = diag((1-mu_z^2)^(1/2))
#   sum = 0
#   for (i in 1:n){
#     sum = sum - 1 / 2 * t(Xt[,i] - xi) %*% solve(Omega) %*% (Xt[,i] - xi) + 
#       log(pnorm(t(alpha) %*% Sigma_z %*% (Xt[,i] - xi)))
#   }
#   temp = - n / 2 * log(det(Omega)) + sum
#   return(-temp)
# }
# opt = optim(par = c(0,0,0), fn = logL, method = 'BFGS')
# gamma = opt$par

m = selm(formula = t(Xt) ~ 1, family = 'SN')
gamma_hat = coef(m, 'CP',vector = F)$gamma1

# m2 = selm(formula = t(Xt) ~ 1, family = 'SN', method = "MPLE")
# gamma_hat = coef(m2, 'CP',vector = F)$gamma1

# coef(m, 'DP', vector = F)
# summary(m)
########## Score function #########
# normal score function
sqx = function(x){
  return(-x)
}

nu  = 8
# t score function
sqx = function(x){
  return(t(apply(x, 1, function(x) -(nu + d) / (nu - 2) * 1 / (1 + 1 / (nu - 2) * sum(x^2)) * x)))
}


# skew-normal score function
Omega = DP$Omega
xi = DP$beta
alpha = DP$alpha
# c = (2 * gamma_hat / (4 - pi))^(1/3)
c = sign(gamma_hat)*(2*abs(gamma_hat)/(4-pi))^(1/3)
mu_z = c / (sqrt(1+c^2))
Sigma_z = diag((1-mu_z^2)^(1/2))
sqx = function(x){
  return(t(apply(x, 1, function(x) - solve(Omega) %*% (x - xi) +
                   as.numeric(dnorm(t(alpha) %*% Sigma_z %*% (x - xi)) /
                                pnorm(t(alpha) %*% Sigma_z %*% (x - xi))) *
                   (Sigma_z %*% alpha)
  )))
  # temp = - solve(Omega) %*% (x - xi) +
  #   as.numeric(dnorm(t(alpha) %*% Sigma_z %*% (x - xi)) /
  #   pnorm(t(alpha) %*% Sigma_z %*% (x - xi))) * (Sigma_z %*% alpha)
}

# Observed Test Statistics
T_obs = KSD_new(x = t(Xt), score_function = sqx, width = -1)
T_obs = as.numeric(T_obs)

# Bootstrap
bootssample_s = NULL
for (b in 1:nB){
  # Resample
  # set.seed(b+1)
  # Recover Ys_b (n x d)
  ita_b = t(GARCH_d.('t', a, B, G, R))
  sim.VAR.par0 = function(ita, coef, order){
    Y = matrix(0, nrow = d, ncol = n + ncut) # d x n
    for (t in (order+1):(n+ncut)){
      sum = rep(0, d)
      for (i in 1:order){
        sum = sum + coef[[i]] %*% Y[,t-i]
      }
      Y[,t] = sum + ita[,t]
    }
    return(Y[,(ncut+1):(n+ncut)])
  }
  Yt_b = t(sim.VAR.par0(t(ita_b), coef, order))
  # Estimation
  # est_s.b = VAR(y = data.frame(Yt_b), p = 3, type = "none")
  
  est_s.b = VAR(data.frame(Yt_b), p = 3, type = "const")
  # est_s.b = restrict(x = est_s.b0, method = "ser", thresh = qnorm(p = 0.99))
  # Yt_b.res = residuals(est_s.b)
  # order = est_s.b$p[[1]]
  # coef = Acoef(est_s.b)
  
  
  # Residual
  Yt_b.res = residuals(est_s.b)
  # Estimation
  est_s1.b = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                      R = R_init, dvar = Yt_b, model = "extended"),
                      error = function(e){return(NA)},
                      warning = function(w) NA)
  if (is.na(est_s1.b[1])) {next}
  
  Yt_b.res1 = est_s1.b$std.resid
  
  ################# M2 ################
  Rs_b = est_s1.b$para.mat$R
  # library(MTS)
  Rs_b.invsqrt = tryCatch(msqrt(Rs_b)$invsqrt,
                          warning = function(w) matrix(NA,d,d),
                          error = function(e)  matrix(NA,d,d))
  if (is.na(Rs_b.invsqrt[1,1])) {next}
  
  Xs_b = Rs_b.invsqrt %*% t(Yt_b.res1)
  # detach("package:MTS", unload=TRUE)
  
  # KSD Test Statistic for bootstrap data
  Ts_b = KSD_new(x = t(Xs_b), score_function = sqx, width = -1)
  Ts_b = as.numeric(Ts_b)
  bootssample_s = c(bootssample_s, Ts_b)
}

pval_s_KSD = sum((bootssample_s - T_obs) > 0) / length(bootssample_s)
# bootssample_s[0.9*1000]

# N
# 0
# t, nu = 7.526442
# 0
# t, nu = 8
# 0
# SN, gamma_hat = (4.873530e-06 -1.807968e-01  3.782506e-03)
# 0