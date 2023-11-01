# Comparsion_ST_GARCH.R
# Multivariate Tests for CCC-GARCH(1,1) model with ST innovations
# Time varying and stochastic variance

library(mvtnorm) # Generate random multivariate normal dist
library(Matrix)
library(MTS)
# library(progress)
# library(KSD)
library(pryr)
library(compositions) # rlnorm.rplus
library(foreach)
library(doSNOW)
library(DEoptim)
# library(micompr) # micompr, cmpoutput()
library(fMultivar) # rcauchy2d
library(statmod) # gauss.quad() used
library(ccgarch)
library(sn)
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

###############   Settings   ###############
# Dimension of data
d = 2
#Mean
mu = rep(0,d)
#Variance-Covariance Matrix
sigma = matrix(c(1,0.5,0.5,1), ncol = 2)
# sigma = matrix(c(1, 0.5, 0.25, 0.125, 0.0625,
#                  0.5, 1, 0.5, 0.25, 0.125,
#                  0.25, 0.5, 1, 0.5, 0.25,
#                  0.125, 0.25, 0.5, 1, 0.5,
#                  0.0625, 0.125, 0.25, 0.5, 1), d, d)
# Sample Size
n = 500
# Repitition Number
nrep = 10000
# Bootstrap Number
nB = 1
# Burn-in Samples Number
ncut = 1000
# d.f. of t dist
nu = 5
# tolerance
tol = .Machine$double.eps^0.25
# Gaussian quadrature weight
weight. = "chebyshev2"

Spec_GARCH = function(N){
  if (N == 2.1){
    # Scenario 1 for d = 2
    a <<- c(0.1, 0.1)
    B <<- matrix(c(0.3, 0.1,
                   0.1, 0.2),2,2)
    G <<- matrix(c(0.2, 0.1,
                   0.01, 0.3),2,2)
    R <<- matrix(c(1, 0.5,
                   0.5, 1),2,2)  
  } else if (N == 2.2){
    # Scenario 2 for d = 2 (More persistence for G)
    a <<- c(0.1, 0.1)
    B <<- matrix(c(0.01, 0.01,
                   0.01, 0.01),2,2)
    G <<- matrix(c(0.85, 0.01,
                   0.01, 0.85),2,2)
    R <<- matrix(c(1, 0.5,
                   0.5, 1),2,2)  
  } else if (N == 3.1){
    # Scenario 1 for d = 3
    a <<- rep(0.1, 3)
    B <<- matrix(c(0.3, 0.1, 0.1,
                   0.1, 0.2, 0.1,
                   0.1, 0.1, 0.1), 3, 3)
    G <<- matrix(c(0.2, 0.1, 0.01,
                   0.01, 0.3, 0.1,
                   0.01, 0.01, 0.1), 3, 3)
    R <<- matrix(0.5, d, d) + diag(1-0.5, d)
  }
}

# Specify initial parameters for estimation
a_init <<- rep(0.005, d)
B_init <<- diag(rep(0.2 - 0.01, d)) + matrix(0.01, d, d)
G_init <<- diag(rep(0.4 - 0.01, d)) + matrix(0.01, d, d)
R_init <<- diag(rep(1 - 0.1, d)) + matrix(0.1, d, d)

Spec_GARCH(2.1)
# distribution mixing factor
lambda = 0.2

GARCH_d = function(dist){
  signal = 0
  while (signal == 0){
    # Initialization
    if (dist == 'N'){
      ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
      # ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 't'){
      ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d), df = nu))
      # ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 'logN'){
      ita = t(matrix(rlnorm.rplus(n + ncut, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n + ncut))
      ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 'mchisq'){
      temp = rWishart(n = n + ncut, df = nu, Sigma = diag(d))
      ita = NULL
      for (i in 1:dim(temp)[1]){
        ita = rbind(ita, temp[i,i,])
      }
      ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 't+N'){
      ita = ((1 - lambda) * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu) +
               lambda * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol")) /
        sqrt(lambda^2 + (1 - lambda)^2)
      ita = t(ita)
    } else if (dist == 't+logN'){
      temp = t(matrix(rlnorm.rplus(n + ncut, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n + ncut))
      temp = t(msqrt(var(t(temp)))$invsqrt %*% (temp - rowMeans(temp)))
      ita = ((1 - lambda) * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu) +
               lambda * temp) / sqrt(lambda^2 + (1 - lambda)^2)
      ita = t(ita)
    } else if (dist == 'SN'){
      ita = t(rmsn(n + ncut, dp = DP))
    } else if (dist == 'ST'){
      ita = t(rmst(n + ncut, dp = DP_ST))
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
      sigma_t.mtxsqrt = tryCatch(msqrt(sigma_t)$mtxsqrt,
                                 warning = function(w) matrix(NA, d, d),
                                 error = function(e) matrix(NA, d, d))
      if (is.na(sigma_t.mtxsqrt[1,1])) break
      eps_t = sigma_t.mtxsqrt %*% ita[,t]
      # eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
      if (t > ncut){y[,t-ncut] = eps_t}
      if (t == (n+ncut)) {signal = 1}
    }
  }
  return(y)
}
GARCH_d. = function(dist, a., B., G., R.){
  signal = 0
  while (signal == 0){
    # Initialization
    if (dist == 'N'){
      ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
      # ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 't'){
      ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d),  df = nu))
      # ita = msqrt(var(t(ita)))$invsqrt %*% (ita - rowMeans(ita))
    } else if (dist == 'SN'){
      ita = t(rmsn(n + ncut, dp = DP))
    } else if (dist == 'ST'){
      ita = t(rmst(n + ncut, dp = DP_ST))
    }
    y = array(rep(0,d*n),dim = c(d,n))
    eps_t = rep(0.5, d)
    sigma_ii = rep(0.3, d)
    C_t = R.
    # Iteration
    for (t in 1:(n+ncut)){
      sigma_ii = a. + B. %*% eps_t^2 + G. %*% sigma_ii
      D_t = diag(c(sqrt(sigma_ii)))
      sigma_t = D_t %*% C_t %*% D_t
      sigma_t.mtxsqrt = tryCatch(msqrt(sigma_t)$mtxsqrt,
                                 warning = function(w) matrix(NA, d, d),
                                 error = function(e) matrix(NA, d, d))
      if (is.na(sigma_t.mtxsqrt[1,1])) break
      eps_t = sigma_t.mtxsqrt %*% ita[,t]
      # eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
      if (t > ncut){y[,t-ncut] = eps_t}
      if (t == (n+ncut)) {signal = 1}
    }
  }
  return(y)
}

# skew-normal parameters
gamma_hat = c(0, -0.6)
# gamma_hat = c(0, 0.2, -0.2, 0, -0.1)
CP = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat)
DP = cp2dp(CP, "SN")
# skew-normal score function
# Omega = DP$Omega
# xi = DP$beta
# alpha = DP$alpha
# # c = (2 * gamma_hat / (4 - pi))^(1/3)
# c = sign(gamma_hat)*(2*abs(gamma_hat)/(4-pi))^(1/3)
# mu_z = c / (sqrt(1+c^2))
# Sigma_z = diag((1-mu_z^2)^(1/2))
# sqx = function(x){
#   return(t(apply(x, 1, function(x) - solve(Omega) %*% (x - xi) +
#                    as.numeric(dnorm(t(alpha) %*% Sigma_z %*% (x - xi)) /
#                                 pnorm(t(alpha) %*% Sigma_z %*% (x - xi))) *
#                    (Sigma_z %*% alpha)
#   )))
# }


# skew-t parameters
# d = 2
CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 17.22322)
DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# d = 5
# CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 70.55935)
# DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# skew-t score function
sqx = function(x){
  h = 0.00001
  df = rep(0, d)
  for (i in 1:d){
    vec01 = rep(0,d)
    vec01[i] = 1
    # five-point method
    df[i] = (-dmst(x = x + 2 * vec01 * h, dp = DP_ST, log = T) +
               8 * dmst(x = x + vec01 * h, dp = DP_ST, log = T) -
               8 * dmst(x = x - vec01 * h, dp = DP_ST, log = T) +
               dmst(x = x - 2 * vec01 * h, dp = DP_ST, log = T)) / 12 / h
  }

  return(df)
}

# Sample Size & Block length
Spec_Nsample = function(N, w = 0.3){
  if (N == 100){
    n <<- 100
    m1 <<- n * w
    m2 <<- n * (1 - w)
  } else if (N == 200){
    n <<- 200
    m1 <<- n * w
    m2 <<- n * (1 - w)
  } else if (N == 500){
    n <<- 500
    m1 <<- n * w
    m2 <<- n * (1 - w)
  }
}
Spec_Nsample(500, 0.3)
# weight
w = 0.3
################### Execution ###################
# pb <- progress_bar$new(
#   format = "  Percent [:bar] :percent Time :elapsed",
#   total = nrep, clear = FALSE, width= 60)

# plist_size_KSD = NULL
# plist_power_KSD = NULL

func = function(cnt){
  sign = 0
  while (sign == 0){
    # for (count in (1:nrep)){
    ################# Generate Data #################
    # set.seed(count + 666)
    # Generate Sample Data
    # ccc_s = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
    #                  d.f = Inf, cut = ncut, model = "diagonal")
    # Ys = ccc_s$eps
    Ys = t(GARCH_d('ST'))
    # ccc_p = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
    #                  d.f = nu, cut = ncut, model = "diagonal")
    # Yp = ccc_p$eps
    Yp = t(GARCH_d('N'))
    ################# Estimation #################
    # est_s = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
    #                             R = R_init, dvar = Ys, model = "extended"),
    #             silent = TRUE)
    # if ('try-error' %in% class(est_s)) {return(c(NA, NA))}
    est_s = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                     R = R_init, dvar = Ys, model = "extended"),
                     error = function(e) {return(NA)},
                     warning = function(w) NA)
    if (is.na(est_s[1])){next}
    # if (is.na(est_s[1])){next}
    # est_p = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
    #                             R = R_init, dvar = Yp, model = "extended"),
    #             silent = TRUE)
    # if ('try-error' %in% class(est_p)) {return(c(NA, NA))}
    est_p = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                     R = R_init, dvar = Yp, model = "extended"),
                     error = function(e) {return(NA)},
                     warning = function(w) {return(NA)})
    if (is.na(est_p[1])) {next}
    # if (is.na(est_p[1])) {next}
    
    Ys.res = est_s$std.resid
    Yp.res = est_p$std.resid
    
    ########################## M2 #########################
    Rs = est_s$para.mat$R
    Rp = est_p$para.mat$R
    
    Rs.invsqrt = tryCatch(msqrt(Rs)$invsqrt,
                          warning = function(w) matrix(NA, d, d),
                          error = function(e) matrix(NA, d, d))
    if (is.na(Rs.invsqrt[1,1])) next
    Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                          warning = function(w) matrix(NA, d, d),
                          error = function(e) matrix(NA, d, d))
    if (is.na(Rp.invsqrt[1,1])) next
    
    Xs = Rs.invsqrt %*% t(Ys.res)
    Xp = Rp.invsqrt %*% t(Yp.res)
    
    ################## Several Monte Carlo Tests ###################
    as = est_s$para.mat$a
    Bs = est_s$para.mat$A
    Gs = est_s$para.mat$B
    ap = est_p$para.mat$a
    Bp = est_p$para.mat$A
    Gp = est_p$para.mat$B
    
    # Observed KSD Test Statistics
    Ts_obs = KSD_new(x = t(Xs), score_function = t(apply(t(Xs), 1, sqx)), width = -1)
    Ts_obs = as.numeric(Ts_obs)
    Tp_obs = KSD_new(x = t(Xp), score_function = t(apply(t(Xp), 1, sqx)), width = -1)
    Tp_obs = as.numeric(Tp_obs)
    
    # Simulate Null dist of T
    # bootssample_s = NULL
    # bootssample_p = NULL
    b = 1
    while (b <= nB){
      ######## Generate residual under H0 ########
      # Recover Ys_b (n x d)
      Ys_b = t(GARCH_d.('ST', as, Bs, Gs, Rs))
      Yp_b = t(GARCH_d.('ST', ap, Bp, Gp, Rp))
      ###### Estimation ######
      # est_s.b = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
      #                             R = R_init, dvar = Ys_b, model = "extended"),
      #             silent = TRUE)
      # if ('try-error' %in% class(est_s.b)) {next}
      est_s.b = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                         R = R_init, dvar = Ys_b, model = "extended"),
                         error = function(e){return(NA)},
                         warning = function(w) NA)
      if (is.na(est_s.b[1])) {next}
      # est_p.b = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
      #                             R = R_init, dvar = Yp_b, model = "extended"),
      #             silent = TRUE)
      # if ('try-error' %in% class(est_p.b)) {next}
      est_p.b = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                         R = R_init, dvar = Yp_b, model = "extended"),
                         error = function(e){return(NA)},
                         warning = function(w) NA)
      if (is.na(est_p.b[1])) {next}
      
      Ys_b.res = est_s.b$std.resid
      Yp_b.res = est_p.b$std.resid
      
      ################# M2 ################
      Rs_b = est_s.b$para.mat$R
      Rp_b = est_p.b$para.mat$R
      
      Rs_b.invsqrt = tryCatch(msqrt(Rs_b)$invsqrt,
                              warning = function(w) matrix(NA,d,d),
                              error = function(e) matrix(NA,d,d))
      if (is.na(Rs_b.invsqrt[1,1])) {next}
      Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
                              warning = function(w) matrix(NA,d,d),
                              error = function(e) matrix(NA,d,d))
      if (is.na(Rp_b.invsqrt[1,1])) {next}
      
      Xs_b = Rs_b.invsqrt %*% t(Ys_b.res)
      Xp_b = Rp_b.invsqrt %*% t(Yp_b.res)
      #####################
      # KSD Test Statistic for bootstrap data
      Ts_b = KSD_new(x = t(Xs_b), score_function = t(apply(t(Xs_b), 1, sqx)), width = -1)
      Ts_b = as.numeric(Ts_b)
      Tp_b = KSD_new(x = t(Xp_b), score_function = t(apply(t(Xp_b), 1, sqx)), width = -1)
      Tp_b = as.numeric(Tp_b)
      
      b = b + 1
    }
    # bootssample_s = c(bootssample_s, Ts_b)
    # bootssample_p = c(bootssample_p, Tp_b)
    
    # sample_s = c(sample_s, Ts_obs)
    # sample_p = c(sample_p, Tp_obs)
    # pval_s_KSD = sum((bootssample_s - Ts_obs) > 0) / length(bootssample_s)
    # pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
    
    
    
    #####################################################
    ################### KSD Block Test ##################
    #####################################################
    # # Block
    # X1s = Xs[,1:m1]
    # X2s = Xs[,(n-m2+1):n]
    # 
    # X1p = Xp[,1:m1]
    # X2p = Xp[,(n-m2+1):n]
    # 
    # # KSD Package Results
    # result1s <- KSD(t(X1s),score_function = sqx, 'rbf',-1.0)
    # result2s <- KSD(t(X2s),score_function = sqx, 'rbf',-1.0)
    # Ts1obs = result1s[["ksd"]]
    # Ts2obs = result2s[["ksd"]]
    # Tsobs = w * Ts1obs + (1-w) * Ts2obs
    # bootssample_s = w * result1s[["bootStrapSamples"]] + (1-w) * result2s[["bootStrapSamples"]]
    # pval_s_KSD_b = sum((bootssample_s - Tsobs) > 0) / 1000 # nboots=1000
    # 
    # result1p <- KSD(t(X1p),score_function = sqx,'rbf',-1.0)
    # result2p <- KSD(t(X2p),score_function = sqx,'rbf',-1.0)
    # Tp1obs = result1p[["ksd"]]
    # Tp2obs = result2p[["ksd"]]
    # Tpobs = w * Tp1obs + (1-w) * Tp2obs
    # bootssample_p = w * result1p[["bootStrapSamples"]] + (1-w) * result2p[["bootStrapSamples"]]
    # pval_p_KSD_b = sum((bootssample_p - Tpobs) > 0) / 1000
    
    
    ######################
    sign = 1
  }
  # ############################ Bai Test ################################
  # #### for t VAR / GARCH
  # sigma_ii.s = est_s[["h"]]
  # Omega_ii.s = (nu - 2) / nu * sigma_ii.s
  # Omega_12.s = sqrt(Omega_ii.s[,1]) * sqrt(Omega_ii.s[,2]) * Rs[1,2]
  # a_s = (nu + Ys[,1]^2 / Omega_ii.s[,1]) / (nu + 1)
  # Omega2_1.s = a_s * (Omega_ii.s[,2] - Omega_12.s^2 / Omega_ii.s[,1])
  # X1s = Ys[,1] / sqrt(Omega_ii.s[,1])
  # X2s = (Ys[,2] - Omega_12.s / Omega_ii.s[,1] * Ys[,1]) / sqrt(Omega2_1.s)
  # U1s = pt(X1s, df = nu)
  # U2s = pt(X2s, df = nu+1)
  # U1s[which(U1s == 1)] = 0.9999999
  # U2s[which(U2s == 1)] = 0.9999999
  # 
  # sigma_ii.p = est_p[["h"]]
  # Omega_ii.p = (nu - 2) / nu * sigma_ii.p
  # Omega_12.p = sqrt(Omega_ii.p[,1]) * sqrt(Omega_ii.p[,2]) * Rp[1,2]
  # a_p = (nu + Yp[,1]^2 / Omega_ii.p[,1]) / (nu + 1)
  # Omega2_1.p = a_p * (Omega_ii.p[,2] - Omega_12.p^2 / Omega_ii.p[,1])
  # X1p = Yp[,1] / sqrt(Omega_ii.p[,1])
  # X2p = (Yp[,2] - Omega_12.p / Omega_ii.p[,1] * Yp[,1]) / sqrt(Omega2_1.p)
  # U1p = pt(X1p, df = nu)
  # U2p = pt(X2p, df = nu+1)
  # U1p[which(U1p == 1)] = 0.9999999
  # U2p[which(U2p == 1)] = 0.9999999
  # # Omega_hats = var(Ys.res) * (nu - 2) / nu
  # # # Omega_hats = Rs * (nu - 2) / nu
  # # a_s = (nu + Ys.res[,1]^2 / Omega_hats[1,1]) / (nu + 1)
  # # Omega2_1s = a_s * (Omega_hats[2,2] - Omega_hats[2,1]^2 / Omega_hats[1,1])
  # # X1s = Ys.res[,1] / sqrt(Omega_hats[1,1])
  # # X2s = (Ys.res[,2] - Omega_hats[2,1] / Omega_hats[1,1] * Ys.res[,1]) / sqrt(Omega2_1s)
  # # U1s = pt(X1s, df = nu)
  # # U2s = pt(X2s, df = nu+1)
  # #
  # # Omega_hatp = var(Yp.res) * (nu - 2) / nu
  # # # Omega_hatp = Rp * (nu - 2) / nu
  # # a_p = (nu + Yp.res[,1]^2 / Omega_hatp[1,1]) / (nu + 1)
  # # Omega2_1p = a_p * (Omega_hatp[2,2] - Omega_hatp[2,1]^2 / Omega_hatp[1,1])
  # # X1p = Yp.res[,1] / sqrt(Omega_hatp[1,1])
  # # X2p = (Yp.res[,2] - Omega_hatp[2,1] / Omega_hatp[1,1] * Yp.res[,1]) / sqrt(Omega2_1p)
  # # U1p = pt(X1p, df = nu)
  # # U2p = pt(X2p, df = nu+1)
  # ################# Functions #################
  # w = function(t, weight, alpha = 0, beta = 0){
  #   if (weight == 'legendre'){
  #     return(1)
  #   } else if (weight == 'chebyshev1'){
  #     return(1 / sqrt(1 - t^2))
  #   } else if (weight == 'chebyshev2'){
  #     return(sqrt(1 - t^2))
  #   } else if (weight == 'jacobi'){
  #     return((1 - t)^alpha * (1 + t)^beta)
  #   }
  # }
  # integral = function(g, a, b, weight = 'legendre', N = 50, alpha. = 0, beta. = 0){
  #   quad = gauss.quad(N, weight, alpha., beta.)
  #   w = as.array(quad$weights)
  #   t = as.array(quad$nodes)
  #   ft = apply(t, 1, function(t) g(0.5 * ((b - a) * t + a + b)) / w(t, weight, alpha = alpha., beta = beta.))
  #   result = (b - a) / 2 * sum(w * ft)
  #   return(result)
  # }
  # J = function(r, n = n, U1 = U1, U2 = U2){
  #   J2n = (sum(U1 - r <= 0) + sum(U2 - r <= 0)) / (2 * n)
  #   return(J2n)
  # }
  # dg = function(r){
  #   return(c(1,
  #            -(nu + 1) * qt(r, nu) / (nu + (qt(r, nu))^2),
  #            1 - (nu + 1) * (qt(r, nu))^2 / (nu + (qt(r, nu))^2),
  #            - (nu + 2) * qt(r, nu + 1) / (nu + 1 + (qt(r, nu + 1))^2),
  #            1 - (nu + 2) * (qt(r, nu + 1))^2 / (nu + 1 + (qt(r, nu + 1))^2)))
  # }
  # C = function(s){
  #   mat = diag(5)
  #   for (i in 1:5){
  #     for (j in 1:5){
  #       mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
  #     }
  #   }
  #   return(mat)
  # }
  # int_dgdJ. = function(s, U1, U2){
  #   ind1 = (U1 >= s)
  #   ind2 = (U2 >= s)
  #   A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 5, 1)
  #   A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 5, 1)
  #   B1 = apply(A1, 1, sum) / (2 * n)
  #   B2 = apply(A2, 1, sum) / (2 * n)
  #   return(B1 + B2)
  # }
  # W_J = function(r, U1, U2){
  #   A = integral(function(s) t(dg(s)) %*% solve(C(s)) %*% int_dgdJ.(s, U1, U2), 0, r)
  #   return(sqrt(2 * n) * (J(r, n, U1, U2) - A))
  # }
  # W_J.s = function(r){
  #   return(- abs(W_J(r, U1s, U2s)))
  # }
  # W_J.p = function(r){
  #   return(- abs(W_J(r, U1p, U2p)))
  # }
  # ctrl = DEoptim.control(itermax = 1, trace = F)
  # Sn.s = -DEoptim(W_J.s, lower = 0, upper = 1, control = ctrl)$optim$bestval
  # Sn.p = -DEoptim(W_J.p, lower = 0, upper = 1, control = ctrl)$optim$bestval
  #######################################
  return(c(
    Ts_obs, Tp_obs,
    Ts_b, Tp_b))
  # Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
  # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
}
# func(2)
# for (cnt in 1:nrep){
#   cat(func(cnt))
# }
################### Size Evaluation ###################
library(pbmcapply)
set.seed(666)
# Sample Size
n = 500
# Repitition Number
nrep = 10000
# Bootstrap Number
nB = 1
lambda = 0.4
result = pbmclapply(1:nrep, func, mc.cores = 40)
# typeof(result)
# bind_rows(as.data.frame(result))
res = NULL
for (i in 1:nrep){
  res = rbind(res, result[[i]])
}
cat("n =", n, "lambda =", lambda)
CV_s = quantile(res[,3], probs = c(0.99, 0.95, 0.90))
CV_p = quantile(res[,4], probs = c(0.99, 0.95, 0.90))
c(sum((res[,1] - CV_s[1]) > 0) / nrep,
  sum((res[,1] - CV_s[2]) > 0) / nrep,
  sum((res[,1] - CV_s[3]) > 0) / nrep,
  sum((res[,2] - CV_p[1]) > 0) / nrep,
  sum((res[,2] - CV_p[2]) > 0) / nrep,
  sum((res[,2] - CV_p[3]) > 0) / nrep)

# write.csv(res, file = paste('GARCH_ST_N', '0', lambda, sep = '_'), row.names = F)



###################################################
###################### POWER ######################
###################################################
func_p = function(cnt){
  sign = 0
  while (sign == 0){
    # for (count in (1:nrep)){
    ################# Generate Data #################
    # set.seed(count + 666)
    # Generate Sample Data
    # ccc_s = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
    #                  d.f = Inf, cut = ncut, model = "diagonal")
    # Ys = ccc_s$eps
    # ccc_p = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
    #                  d.f = nu, cut = ncut, model = "diagonal")
    # Yp = ccc_p$eps
    Yp = t(GARCH_d('t+N'))
    ################# Estimation #################
    # est_s = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
    #                             R = R_init, dvar = Ys, model = "extended"),
    #             silent = TRUE)
    # if ('try-error' %in% class(est_s)) {return(c(NA, NA))}
    # if (is.na(est_s[1])){next}
    # est_p = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
    #                             R = R_init, dvar = Yp, model = "extended"),
    #             silent = TRUE)
    # if ('try-error' %in% class(est_p)) {return(c(NA, NA))}
    est_p = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                     R = R_init, dvar = Yp, model = "extended"),
                     error = function(e) {return(NA)},
                     warning = function(w) {return(NA)})
    if (is.na(est_p[1])) {next}
    # if (is.na(est_p[1])) {next}
    
    Yp.res = est_p$std.resid
    
    ########################## M1 ##########################
    # #Estimated Mean
    # mu_hats = colMeans(Ys.res)
    # mu_hatp = colMeans(Yp.res)
    # #Estimated Sigma
    # sigma_hats = var(Ys.res)
    # sigma_hatp = var(Yp.res)
    # 
    # sigma_hats.invsqrt = tryCatch(msqrt(sigma_hats)$invsqrt,
    #                       warning = function(w) matrix(NA, d, d),
    #                       error = function(e) matrix(NA, d, d))
    # if (is.na(sigma_hats.invsqrt[1,1])) return(c(NA, NA))
    # sigma_hatp.invsqrt = tryCatch(msqrt(sigma_hatp)$invsqrt,
    #                       warning = function(w) matrix(NA, d, d),
    #                       error = function(e) matrix(NA, d, d))
    # if (is.na(sigma_hatp.invsqrt[1,1])) return(c(NA, NA))
    # 
    # Xs = sigma_hats.invsqrt %*% (t(Ys.res) - mu_hats)
    # Xp = sigma_hatp.invsqrt %*% (t(Yp.res) - mu_hatp)
    
    ########################## M2 #########################
    Rp = est_p$para.mat$R
    
    Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                          warning = function(w) matrix(NA, d, d),
                          error = function(e) matrix(NA, d, d))
    if (is.na(Rp.invsqrt[1,1])) next
    
    Xp = Rp.invsqrt %*% t(Yp.res)
    
    ################### Several Monte Carlo Tests ###################
    ap = est_p$para.mat$a
    Bp = est_p$para.mat$A
    Gp = est_p$para.mat$B
    
    # Observed KSD Test Statistics
    Tp_obs = KSD_new(x = t(Xp), score_function = t(apply(t(Xp), 1, sqx)), width = -1)
    Tp_obs = as.numeric(Tp_obs)
    
    b = 1
    while (b <= nB){
      ######## Generate residual under H0 ########
      # Recover Ys_b (n x d)
      Yp_b = t(GARCH_d.('ST', ap, Bp, Gp, Rp))
      ###### Estimation ######
      # est_s.b = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
      #                             R = R_init, dvar = Ys_b, model = "extended"),
      #             silent = TRUE)
      # if ('try-error' %in% class(est_s.b)) {next}
      # est_p.b = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
      #                             R = R_init, dvar = Yp_b, model = "extended"),
      #             silent = TRUE)
      # if ('try-error' %in% class(est_p.b)) {next}
      est_p.b = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                         R = R_init, dvar = Yp_b, model = "extended"),
                         error = function(e){return(NA)},
                         warning = function(w) NA)
      if (is.na(est_p.b[1])) {next}
      
      Yp_b.res = est_p.b$std.resid
      
      ################# M2 ################
      Rp_b = est_p.b$para.mat$R
      
      Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
                              warning = function(w) matrix(NA,d,d),
                              error = function(e) matrix(NA,d,d))
      if (is.na(Rp_b.invsqrt[1,1])) {next}
      Xp_b = Rp_b.invsqrt %*% t(Yp_b.res)
      #####################
      # KSD Test Statistic for bootstrap data
      Tp_b = KSD_new(x = t(Xp_b), score_function = t(apply(t(Xp_b), 1, sqx)), width = -1)
      Tp_b = as.numeric(Tp_b)
      
      b = b + 1
    }
    
    #####################################################
    ################### KSD Block Test ##################
    #####################################################
    # # Block
    # X1p = Xp[,1:m1]
    # X2p = Xp[,(n-m2+1):n]
    # 
    # # KSD Package Results
    # result1p <- KSD(t(X1p),score_function = sqx,'rbf',-1.0)
    # result2p <- KSD(t(X2p),score_function = sqx,'rbf',-1.0)
    # Tp1obs = result1p[["ksd"]]
    # Tp2obs = result2p[["ksd"]]
    # Tpobs = w * Tp1obs + (1-w) * Tp2obs
    # bootssample_p = w * result1p[["bootStrapSamples"]] + (1-w) * result2p[["bootStrapSamples"]]
    # pval_p_KSD_b = sum((bootssample_p - Tpobs) > 0) / 1000
    
    
    ######################
    sign = 1
  }
  ############################ Bai Test ################################
  #### for t VAR / GARCH
  # sigma_ii.p = est_p[["h"]]
  # Omega_ii.p = (nu - 2) / nu * sigma_ii.p
  # Omega_12.p = sqrt(Omega_ii.p[,1]) * sqrt(Omega_ii.p[,2]) * Rp[1,2]
  # a_p = (nu + Yp[,1]^2 / Omega_ii.p[,1]) / (nu + 1)
  # Omega2_1.p = a_p * (Omega_ii.p[,2] - Omega_12.p^2 / Omega_ii.p[,1])
  # X1p = Yp[,1] / sqrt(Omega_ii.p[,1])
  # X2p = (Yp[,2] - Omega_12.p / Omega_ii.p[,1] * Yp[,1]) / sqrt(Omega2_1.p)
  # U1p = pt(X1p, df = nu)
  # U2p = pt(X2p, df = nu+1)
  # U1p[which(U1p == 1)] = 0.9999999
  # U2p[which(U2p == 1)] = 0.9999999
  # # Omega_hats = var(Ys.res) * (nu - 2) / nu
  # # # Omega_hats = Rs * (nu - 2) / nu
  # # a_s = (nu + Ys.res[,1]^2 / Omega_hats[1,1]) / (nu + 1)
  # # Omega2_1s = a_s * (Omega_hats[2,2] - Omega_hats[2,1]^2 / Omega_hats[1,1])
  # # X1s = Ys.res[,1] / sqrt(Omega_hats[1,1])
  # # X2s = (Ys.res[,2] - Omega_hats[2,1] / Omega_hats[1,1] * Ys.res[,1]) / sqrt(Omega2_1s)
  # # U1s = pt(X1s, df = nu)
  # # U2s = pt(X2s, df = nu+1)
  # # 
  # # Omega_hatp = var(Yp.res) * (nu - 2) / nu
  # # # Omega_hatp = Rp * (nu - 2) / nu
  # # a_p = (nu + Yp.res[,1]^2 / Omega_hatp[1,1]) / (nu + 1)
  # # Omega2_1p = a_p * (Omega_hatp[2,2] - Omega_hatp[2,1]^2 / Omega_hatp[1,1])
  # # X1p = Yp.res[,1] / sqrt(Omega_hatp[1,1])
  # # X2p = (Yp.res[,2] - Omega_hatp[2,1] / Omega_hatp[1,1] * Yp.res[,1]) / sqrt(Omega2_1p)
  # # U1p = pt(X1p, df = nu)
  # # U2p = pt(X2p, df = nu+1)
  # # ################# Functions #################
  # w = function(t, weight, alpha = 0, beta = 0){
  #   if (weight == 'legendre'){
  #     return(1)
  #   } else if (weight == 'chebyshev1'){
  #     return(1 / sqrt(1 - t^2))
  #   } else if (weight == 'chebyshev2'){
  #     return(sqrt(1 - t^2))
  #   } else if (weight == 'jacobi'){
  #     return((1 - t)^alpha * (1 + t)^beta)
  #   }
  # }
  # integral = function(g, a, b, weight = 'legendre', N = 50, alpha. = 0, beta. = 0){
  #   quad = gauss.quad(N, weight, alpha., beta.)
  #   w = as.array(quad$weights)
  #   t = as.array(quad$nodes)
  #   ft = apply(t, 1, function(t) g(0.5 * ((b - a) * t + a + b)) / w(t, weight, alpha = alpha., beta = beta.))
  #   result = (b - a) / 2 * sum(w * ft)
  #   return(result)
  # }
  # J = function(r, n = n, U1 = U1, U2 = U2){
  #   J2n = (sum(U1 - r <= 0) + sum(U2 - r <= 0)) / (2 * n)
  #   return(J2n)
  # }
  # dg = function(r){
  #   return(c(1,
  #            -(nu + 1) * qt(r, nu) / (nu + (qt(r, nu))^2),
  #            1 - (nu + 1) * (qt(r, nu))^2 / (nu + (qt(r, nu))^2),
  #            - (nu + 2) * qt(r, nu + 1) / (nu + 1 + (qt(r, nu + 1))^2),
  #            1 - (nu + 2) * (qt(r, nu + 1))^2 / (nu + 1 + (qt(r, nu + 1))^2)))
  # }
  # C = function(s){
  #   mat = diag(5)
  #   for (i in 1:5){
  #     for (j in 1:5){
  #       mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
  #     }
  #   }
  #   return(mat)
  # }
  # int_dgdJ. = function(s, U1, U2){
  #   ind1 = (U1 >= s)
  #   ind2 = (U2 >= s)
  #   A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 5, 1)
  #   A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 5, 1)
  #   B1 = apply(A1, 1, sum) / (2 * n)
  #   B2 = apply(A2, 1, sum) / (2 * n)
  #   return(B1 + B2)
  # }
  # W_J = function(r, U1, U2){
  #   A = integral(function(s) t(dg(s)) %*% solve(C(s)) %*% int_dgdJ.(s, U1, U2), 0, r)
  #   return(sqrt(2 * n) * (J(r, n, U1, U2) - A))
  # }
  # W_J.p = function(r){
  #   return(- abs(W_J(r, U1p, U2p)))
  # }
  # ctrl = DEoptim.control(itermax = 1, trace = F)
  # Sn.p = -DEoptim(W_J.p, lower = 0, upper = 1, control = ctrl)$optim$bestval
  #######################################
  return(c(
    # pval_s_skew, pval_p_skew,
    # pval_s_kurt, pval_p_kurt,
    # pval_s_skku, pval_p_skku,
    # pval_s_sk, pval_p_sk,
    # pval_s_ku, pval_p_ku,
    # pval_s_JB, pval_p_JB,
    # pval_p_KSD,
    # pval_p_KSD_b,
    Tp_obs,
    Tp_b))
  # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
}

library(pbmcapply)
set.seed(666)
# Sample Size
n = 500
# Repitition Number
nrep = 10000
# Bootstrap Number
nB = 1
lambda = 1
result = pbmclapply(1:nrep, func_p, mc.cores = 40)
# typeof(result)
# bind_rows(as.data.frame(result))
res = NULL
for (i in 1:nrep){
  res = rbind(res, result[[i]])
}
cat("n =", n, "lambda =", lambda)
CV_p = quantile(res[,2], probs = c(0.99, 0.95, 0.90))
c(sum((res[,1] - CV_p[1]) > 0) / nrep,
  sum((res[,1] - CV_p[2]) > 0) / nrep,
  sum((res[,1] - CV_p[3]) > 0) / nrep)

# df = data.frame(result)
# cat("n = ",n, "lambda =", lambda)
# cbind(t(t(apply(df[1,] < 0.01, 1, sum) / nrep)),
#       t(t(apply(df[1,] < 0.05, 1, sum) / nrep)),
#       t(t(apply(df[1,] < 0.10, 1, sum) / nrep)))
# apply(df[2:4,], 1, sum) / nrep

write.csv(res, file = paste('GARCH_ST_logN', lambda, sep = "_"))


########################################################
# ---------- Use foreach and doSNOW package ---------- #
eval_test = function(plist_power){
  #Remove NAs
  np = sum(is.na(plist_power))
  plist_power[is.na(plist_power)] = 1
  power1 = sum(plist_power - 0.01 < 0) / (length(plist_power) - np)
  power5 = sum(plist_power - 0.05 < 0) / (length(plist_power) - np)
  power10 = sum(plist_power - 0.10 < 0) / (length(plist_power) - np)
  print(c(power1, power5, power10))
}
cl <- makeCluster(40)
registerDoSNOW(cl)
lambda = 1
# clusterSetRNGStream(iseed = 666)
pb = txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
result = foreach(cnt = 1:nrep,
                 .combine = 'rbind',
                 .options.snow = opts,
                 .packages = c("mvtnorm",
                               "Matrix",
                               "MTS",
                               "pryr",
                               "compositions",
                               "ccgarch",
                               "DEoptim",
                               "statmod",
                               "fMultivar",
                               "sn")) %dopar% {
                                 tryCatch(func_p(cnt),error = function(e) return(NA))
                               }
# func(cnt)
# {
#                                a = tryCatch(func(cnt),error = function(e) return(e))
#                                print(a)
#                              }


close(pb)
stopCluster(cl)
# stopImplicitCluster()
write.csv(result, file = paste('GARCH_ST_logN', lambda, sep = "_"))
# write.csv(result, file = paste('GARCH_t_N', lambda, sep = '_'))
plist_power_KSD = result[,1]
plist_power_Bai = result[,2:4]

eval_test(plist_power_KSD)
np = colSums(is.na(plist_power_Bai))[1]
plist_power_Bai[is.na(plist_power_Bai)] = 0
colSums(plist_power_Bai) / (length(plist_power_Bai[,1]) - np)



