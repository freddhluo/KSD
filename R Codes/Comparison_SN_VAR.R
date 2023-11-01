# Comparsion_VAR.R
# Multivariate Tests for VAR(p) model
# Time varying and stochastic mean

library(mvnTest) # Test for normality
library(MVN) # Test for normality
library(energy) # Test for normality
library(mvtnorm) # Generate random multivariate normal dist
library(Matrix)
library(MTS)
library(progress)
# library(KSD)
library(pryr)
library(compositions)
# library(MVN)  
# library(energy)
library(foreach)
library(doSNOW)
library(vars)

library(DEoptim)
library(fMultivar) # rcauchy2d
library(statmod) # gauss.quad() used

library(sn) # SN distribution

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
n = 100
# Repitition Number
nrep = 1000       
# Bootstrap Number
nB = 1000
# Burn-in Samples Number
ncut = 1000
# VAR order
VAR.p = 3
# d.f. of t dist
nu = 5

# gamma_hat = c(0, 0.2, -0.2, 0, -0.1)
gamma_hat = c(0, -0.6)
CP = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat)
DP = cp2dp(CP, "SN")
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
}

# skew-t parameters
# d = 2
CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 17.22322)
DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# d = 5
# CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 70.55935)
# DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu

# tolerance
tol = .Machine$double.eps^0.25
# quadrature rule
# weight. = "legendre"
# weight. = "chebyshev1"
weight. = "chebyshev2"
# weight. = "jacobi"
# alpha & beta if weight. = "jacobi"
# alpha = 0
# beta = 0
lambda = 1

# Specify VAR Parameters
A = array(c(0.3,-0.2,0.65,-0.4,
               -0.4,-0.6,0.4,0.4,
               0.5,0.1,0.1,0.5), c(d, d, VAR.p))

# Generate VAR parameters
# A = array(0, c(d, d, VAR.p))
# Q = array(0, c(d, d, VAR.p))
# lambda0 = array(0, c(d, VAR.p))
# # Specify eigenvalues
# lambda0[,1] = c(0.3, -0.6, 0.2, 0.35, 0.65)
# lambda0[,2] = c(0.2, 0.45, 0.1, -0.4, 0.5)
# lambda0[,3] = c(0.7, 0.3, 0.4, -0.2, -0.5)
# 
# # Get Q (eigenvectors)
# for (k in 1:VAR.p){
#   Q[,,k] = runif(d*d, 0.4, 0.6)
#   Q[,,k] = apply(Q[,,k], 2, function(v) v / norm(v))
# }
# # Get A (AR parameters)
# for (i in 1:VAR.p){
#   A[,,i] = Q[,,i] %*% diag(lambda0[,i]) %*% solve(Q[,,i])
# }
# A
# apply(A,3,eigen)



# Function to simulate VAR data
# current version
# Generate sample data for evaluation
sim.VAR0 = function(dist, A, sigma_u){
  # Simulation
  if (dist == 'N'){
    ita = rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol")
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 't'){
    ita = rmvt(n + ncut, delta = rep(0,2), sigma = (5 - 2) / 5 * diag(2), df = 5)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 't+N'){
    ita = ((1 - lambda) * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu) +
             lambda * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol")) /
      sqrt(lambda^2 + (1 - lambda)^2)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 't+logN'){
    temp = t(matrix(rlnorm.rplus(n + ncut, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n + ncut))
    temp = t(msqrt(var(t(temp)))$invsqrt %*% (temp - rowMeans(temp)))
    ita = ((1 - lambda) * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu) +
             lambda * temp) / sqrt(lambda^2 + (1 - lambda)^2)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 'SN'){
    ita = rmsn(n + ncut, dp = DP)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 'ST'){
    ita = rmst(n + ncut, dp = DP_ST)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  }
  Y = matrix(0, nrow = d, ncol = n + ncut) # d x n
  for (t in (VAR.p+1):(n+ncut)){
    sum = rep(0, d)
    for (i in 1:VAR.p){
      sum = sum + A[,,i] %*% Y[,t-i]
    }
    Y[,t] = sum + ita[,t]
  }
  return(Y[,(ncut+1):(n+ncut)])
}
# Simulate data based on observed sample to mimic the dist of test statistic
sim.VAR.par0 = function(dist, coef, mu_hat, R, order){
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d),  df = nu))
  } else if (dist == 'SN'){
    ita = t(rmsn(n + ncut, dp = DP))
  } else if (dist == 'ST'){
    ita = t(rmst(n + ncut, dp = DP_ST))
  }
  ita = mu_hat + msqrt(R)$mtxsqrt %*% ita
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

# tolerance
tol = .Machine$double.eps^0.25
# quadrature rule
# weight. = "legendre"
# weight. = "chebyshev1"
weight. = "chebyshev2"
# weight. = "jacobi"
# alpha & beta if weight. = "jacobi"
# alpha = 0
# beta = 0

################### Execution ###################
# pb <- progress_bar$new(
#   format = "  Percent [:bar] :percent Time :elapsed",
#   total = nrep, clear = FALSE, width= 60)

# plist_size_KSD = NULL
# plist_power_KSD = NULL

func = function(cnt){
  # for (count in (1:nrep)){
  ################# Generate Data #################
  # set.seed(cnt + 12345)
  # Generate Sample Data
  # Ys = t(sim.VAR('N', coef))
  Ys = t(sim.VAR0('SN', A, sigma))
  # Yp = t(sim.VAR('t', coef))
  Yp = t(sim.VAR0('t', A, sigma))
  # plot(x = 1:n, y = Yp[,1])
  ################Estimation#################
  # est_s = VAR(y = data.frame(Ys), lag.max = 5, ic = "AIC", type = "none")
  # est_p = VAR(y = data.frame(Yp), lag.max = 5, ic = "AIC", type = "none")
  est_s = VAR(y = data.frame(Ys), p = VAR.p, type = "none")
  est_p = VAR(y = data.frame(Yp), p = VAR.p, type = "none")
  
  order.s = est_s$p[[1]]
  order.p = est_p$p[[1]]
  
  coef.s = Acoef(est_s)
  coef.p = Acoef(est_p)
  
  Ys.res = residuals(est_s)
  Yp.res = residuals(est_p)
  
  #Estimated Mean
  mu_hats = colMeans(Ys.res)
  mu_hatp = colMeans(Yp.res)
  #Estimated Sigma
  Rs = var(Ys.res)
  Rp = var(Yp.res)
  # Invsqrt of R
  Rs.invsqrt = tryCatch(msqrt(Rs)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  # if (is.na(Rs.invsqrt)) {next}
  if (is.na(Rs.invsqrt[1,1])) return(NA)
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  # if (is.na(Rp.invsqrt)) {next}
  if (is.na(Rp.invsqrt[1,1])) return(NA)
  
  Xs = Rs.invsqrt %*% (t(Ys.res) - mu_hats)
  Xp = Rp.invsqrt %*% (t(Yp.res) - mu_hatp)
  
  # Xs = t(residuals(est_s))
  # Xp = t(residuals(est_p))
  
  # mu_hats = rowMeans(Xs)
  # mu_hatp = rowMeans(Xp)
  
  # Rs = var(t(Xs))
  # Rp = var(t(Xp))
  
  # Rs.invsqrt = tryCatch(msqrt(Rs)$invsqrt,
  #                       warning = function(w) NA,
  #                       error = function(e) NA)
  # if (is.na(Rs.invsqrt)) {next}
  # if (is.na(Rs.invsqrt)) return(c(NA, NA))
  # Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
  #                       warning = function(w) NA,
  #                       error = function(e) NA)
  # if (is.na(Rp.invsqrt)) {next}
  # if (is.na(Rp.invsqrt)) return(c(NA, NA))
  # Xs = Rs.invsqrt %*% Xs
  # Xp = Rp.invsqrt %*% Xp
  
  ################### KSD Test ###################
  # Observed Test Statistics
  Ts_obs = KSD_new(x = t(Xs), score_function = sqx, width = -1)
  Ts_obs = as.numeric(Ts_obs)
  Tp_obs = KSD_new(x = t(Xp), score_function = sqx, width = -1)
  Tp_obs = as.numeric(Tp_obs)
  
  # Simulate Null dist of T
  bootssample_s = NULL
  bootssample_p = NULL
  for (b in 1:nB){
    # Generate residual under H0
    # Recover Ys_b (n x d)
    Ys_b = t(sim.VAR.par0('SN', coef.s, mu_hats, Rs, order.s))
    Yp_b = t(sim.VAR.par0('SN', coef.p, mu_hatp, Rp, order.p))
    # Estimation
    # est_s.b = VAR(y = data.frame(Ys_b), lag.max = 5, ic = "SC", type = "none")
    # est_p.b = VAR(y = data.frame(Yp_b), lag.max = 5, ic = "SC", type = "none")
    est_s.b = VAR(y = data.frame(Ys_b), p = VAR.p, type = "none")
    est_p.b = VAR(y = data.frame(Yp_b), p = VAR.p, type = "none")
    # Residual
    Ys_b.res = residuals(est_s.b)
    Yp_b.res = residuals(est_p.b)
    #Estimated Mean
    mu_hats_b = colMeans(Ys_b.res)
    mu_hatp_b = colMeans(Yp_b.res)
    #Estimated Sigma
    Rs_b = var(Ys_b.res)
    Rp_b = var(Yp_b.res)
    # Invsqrt of R
    Rs_b.invsqrt = tryCatch(msqrt(Rs_b)$invsqrt,
                            warning = function(w) matrix(NA, d, d),
                            error = function(e) matrix(NA, d, d))
    if (is.na(Rs_b.invsqrt[1,1])) next
    Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
                            warning = function(w) matrix(NA, d, d),
                            error = function(e) matrix(NA, d, d))
    if (is.na(Rp_b.invsqrt[1,1])) next
    # Standardized Residual
    Xs_b = Rs_b.invsqrt %*% (t(Ys_b.res) - mu_hats_b)
    Xp_b = Rp_b.invsqrt %*% (t(Yp_b.res) - mu_hatp_b)
    # # Estimation
    # est_s.b = VAR(y = data.frame(Ys_b), lag.max = 5, ic = "SC", type = "none")
    # est_p.b = VAR(y = data.frame(Yp_b), lag.max = 5, ic = "SC", type = "none")
    # Xs_b = t(residuals(est_s.b))
    # Xp_b = t(residuals(est_p.b))
    # Rs_b = var(t(Xs_b))
    # Rp_b = var(t(Xp_b))
    # Rs_b.invsqrt = tryCatch(msqrt(Rs_b)$invsqrt,
    #                         warning = function(w) NA,
    #                         error = function(e) NA)
    # if (is.na(Rs_b.invsqrt)) {next}
    # Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
    #                         warning = function(w) NA,
    #                         error = function(e) NA)
    # if (is.na(Rp_b.invsqrt)) {next}
    # Xs_b = Rs_b.invsqrt %*% Xs_b
    # Xp_b = Rp_b.invsqrt %*% Xp_b
    
    # KSD Test Statistic for bootstrap data
    Ts_b = KSD_new(x = t(Xs_b), score_function = sqx, width = -1)
    Ts_b = as.numeric(Ts_b)
    Tp_b = KSD_new(x = t(Xp_b), score_function = sqx, width = -1)
    Tp_b = as.numeric(Tp_b)
    bootssample_s = c(bootssample_s, Ts_b)
    bootssample_p = c(bootssample_p, Tp_b)
  }
  
  pval_s_KSD = sum((bootssample_s - Ts_obs) > 0) / length(bootssample_s)
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
  
  # plist_size_KSD = c(plist_size_KSD, pval_s_KSD)
  # plist_power_KSD = c(plist_power_KSD, pval_p_KSD)
  #######################################
  return(c(pval_s_KSD, pval_p_KSD))
           # Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
  # pb$tick()
  # Sys.sleep(1 / 100)
}



n = 100
# lambda = 0.2
cl <- makeCluster(64)
registerDoSNOW(cl)
pb = txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
result = foreach(cnt = 1:nrep,
                 .combine = 'rbind',
                 .options.snow = opts,
                 .packages = c("mvtnorm",
                               "Matrix",
                               "MTS",
                               "progress",
                               "pryr",
                               "compositions",
                               "MVN",
                               "energy",
                               'mvnTest',
                               "DEoptim",
                               "fMultivar",
                               "statmod",
                               "vars",
                               "sn")) %dopar% {
                                 tryCatch(func(cnt),error = function(e) return(NA))
                               }

close(pb)
stopCluster(cl)
# result
cat("n =",n, "lambda =",lambda)
cbind(apply(result[,1:2]<0.01, 2, sum) / nrep,
      apply(result[,1:2]<0.05, 2, sum) / nrep,
      apply(result[,1:2]<0.10, 2, sum) / nrep)
# apply(result[,29:34], 2, sum) / nrep
# which(is.na(result))

