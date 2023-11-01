# Comparsion_VAR.R
# Multivariate Normality Tests for VAR(p) model
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

# Henze test statistic
henze = function(X){
  n = dim(X)[2]
  temp = 0
  for (j in 1:n){
    for (k in 1:n){
      temp = temp + exp((sum(X[,j]^2) - sum(X[,k]^2)) / 4 / gamma) *
        cos((t(X[,j]) %*% X[,k]) / 2 / gamma)
    }
  }
  Tn = sqrt(n) * (temp / n^2 - 1)
  return(Tn)
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
# A = array(c(0.5,-0.2,0.65,-0.4,
#             -0.4,-0.6,0.4,0.4,
#             0.5,0.1,0.1,0.5), c(d, d, VAR.p))
# A1 = matrix(c(0.2,0.1,-0.2,0, 0,
#                  0,-0.3,0.1,-0.1,0,
#                  0,0.05,0.15,0,0,
#                  -0.05,0,0.1,-0.2,0,
#                  0.05,-0.1,-0.1,0,0.3),d,d,byrow = T)
# A2 = matrix(c(0.25,0.05,0.1,0, 0,
#               -0.2,0.1,0.1,0,0,
#               0.1,0.1,-0.2,0,0,
#               0,0,0,-0.1,0.1,
#               0,0,0,0.2,0.3),d,d,byrow = T)
# A3 = matrix(c(-0.3,0.05,0.1,0, 0,
#               -0.2,0.2,0.1,0,0,
#               0.05,-0.1,0.2,0,0,
#               0,0,0,-0.15,-0.1,
#               0,0,0,0.05,0.2),d,d,byrow = T)
# A = array(0, c(d,d,VAR.p))
# A[,,1]=A1
# A[,,2]=A2
# A[,,3]=A3

# # Generate VAR parameters
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
# # A
# # apply(A,3,eigen)

# skew_normal parameters
# gamma_hat = c(0, 0.2, -0.2, 0, -0.1)
gamma_hat = c(0, -0.6)
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

# Function to simulate VAR data
# not in use
sim.VAR = function(dist, coef. = coef){
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = mu, sigma = sigma, method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = mu, sigma = (nu-2)/nu * sigma,  df = nu))
  }
  Y = matrix(0, nrow = d, ncol = n + ncut) # d x n
  for (t in (VAR.p+1):(n+ncut)){
    sum = rep(0, d)
    for (i in 1:VAR.p){
      sum = sum + coef.[,,i] %*% Y[,t-i]
    }
    Y[,t] = sum + ita[,t]
  }
  return(Y[,(ncut+1):(n+ncut)])
}
# current version
sim.VAR0 = function(dist, A, sigma_u){
  # # Generate mu
  # mu = 0
  # # mu = runif(1, -10, 10)
  # # Generate sigma_u
  # scale = 2 * runif(d, 0, 1) + 1e-30
  # sigma_u = diag(scale) + matrix(runif(d^2, 0, 1), d, d)
  # sigma_u = t(sigma_u) %*% sigma_u
  # while (any(eigen(sigma_u, only.values = T)$values < 0)){
  #   scale = 2 * runif(d, 0, 1) + 1e-30
  #   sigma_u = diag(scale) + matrix(runif(d^2, 0, 1), d, d)
  #   sigma_u = t(sigma_u) %*% sigma_u
  # }
  # # Generate A (d x d x VAR.p)
  # A = array(0, c(d,d,VAR.p))
  # lambda = 2.5
  # for (i in 1:VAR.p){
  #   A[,,i] = lambda^(-i) * matrix(runif(d^2, 0, 1), d, d) - 
  #     (2 * lambda)^(-i) * matrix(1, d, d)
  # }
  # eigenval = apply(A, 3, function(M) eigen(M, only.values = T)$values)
  # while(max(abs(eigenval)) >= 1){
  #   A = array(0, c(d,d,VAR.p))
  #   lambda = 2.5
  #   for (i in 1:VAR.p){
  #     A[,,i] = lambda^(-i) * matrix(runif(d^2, 0, 6), d, d) - 
  #       (2 * lambda)^(-i) * matrix(3, d, d)
  #   }
  # }
  # Simulation
  if (dist == 'N'){
    ita = rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol")
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
    # var(t(ita))
    # ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = sigma_u, method = "chol"))
  } else if (dist == 't'){
    ita = rmvt(n + ncut, delta = rep(0,d), sigma = (5 - 2) / 5 * diag(d), df = 5)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
    # ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * sigma_u,  df = nu))
  } else if (dist == 'N+t'){
    ita = ((1 - lambda) * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol") +
            lambda * rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)) / 
      sqrt(lambda^2 + (1 - lambda)^2)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
  } else if (dist == 'N+logN'){
    temp = t(matrix(rlnorm.rplus(n + ncut, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n + ncut))
    temp = t(msqrt(var(t(temp)))$invsqrt %*% (temp - rowMeans(temp)))
    ita = ((1 - lambda) * rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol") +
            lambda * temp) / 
      sqrt(lambda^2 + (1 - lambda)^2)
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
# e.g. Y = sim.VAR(coef)
# not in use
sim.VAR.par = function(dist, coef, order){
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = rep(0,d), sigma = (nu-2)/nu * diag(d),  df = nu))
  }
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
# current version
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
  Ys = t(sim.VAR0('N', A, sigma))
  # Yp = t(sim.VAR('t', coef))
  Yp = t(sim.VAR0('SN', A, sigma))
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
  
  ################### Skew / Kurt Test ###################
  # b1_hats = apply(Xs^3, 1, mean)
  # b1_hatp = apply(Xp^3, 1, mean)
  # b2_hats = apply(Xs^4, 1, mean)
  # b2_hatp = apply(Xp^4, 1, mean)
  # 
  # lambda_s.sk = (n - order.s) * sum(b1_hats^2) / 6 
  # lambda_p.sk = (n - order.p) * sum(b1_hatp^2) / 6
  # lambda_s.ku = (n - order.s) * sum((b2_hats - 3)^2) / 24
  # lambda_p.ku = (n - order.p) * sum((b2_hatp - 3)^2) / 24
  # 
  # pval_s_skew = 1 - pchisq(lambda_s.sk, d)
  # pval_s_kurt = 1 - pchisq(lambda_s.ku, d)
  # pval_p_skew = 1 - pchisq(lambda_p.sk, d)
  # pval_p_kurt = 1 - pchisq(lambda_p.ku, d)
  # pval_s_skku = 1 - pchisq(lambda_s.sk + lambda_s.ku, 2 * d)
  # pval_p_skku = 1 - pchisq(lambda_p.sk + lambda_p.ku, 2 * d)
  
  ################### DH Test ###################
  # result_s_DH = DH.test(data = t(Xs), qqplot = FALSE)
  # result_p_DH = DH.test(data = t(Xp), qqplot = FALSE)
  ################### AD Test ###################
  # result_s_AD = AD.test(data = t(Xs), qqplot = FALSE)
  # result_p_AD = AD.test(data = t(Xp), qqplot = FALSE)
  ################### CM Test ###################
  # result_s_CM = CM.test(data = t(Xs), qqplot = FALSE)
  # result_p_CM = CM.test(data = t(Xp), qqplot = FALSE)
  ################### HZ Test ###################
  # result_s_HZ = HZ.test(data = t(Xs), qqplot = FALSE)
  # result_p_HZ = HZ.test(data = t(Xp), qqplot = FALSE)
  ################### R Test ###################
  # result_s_R = R.test(data = t(Xs), qqplot = FALSE)
  # result_p_R = R.test(data = t(Xp), qqplot = FALSE)
  ################### Mardia Test ###################
  # result_s_mardia = mvn(data = t(Xs), mvnTest = "mardia")
  # result_p_mardia = mvn(data = t(Xp), mvnTest = "mardia")
  ################### Henze Test ###################
  # result_s_henze = mvn(data = t(Xs), mvnTest = "hz")
  # result_p_henze = mvn(data = t(Xp), mvnTest = "hz")
  ################### Roy Test ###################
  # result_s_roy = mvn(data = t(Xs), mvnTest = "royston")
  # result_p_roy = mvn(data = t(Xp), mvnTest = "royston")
  ################### doh Test ###################
  # result_s_doh = mvn(data = t(Xs), mvnTest = "dh")
  # result_p_doh = mvn(data = t(Xp), mvnTest = "dh")
  ################### Energy Test ###################
  # result_s_energy = mvn(data = t(Xs), mvnTest = "energy")
  # result_p_energy = mvn(data = t(Xp), mvnTest = "energy")
  
  # result_s_energy = mvnorm.etest(t(Xs),1000)
  # result_p_energy = mvnorm.etest(t(Xp),1000)
  ################### JB Test ###################
  # g.Xs = Xs^3 - 3 * Xs
  # JB1_s = sum((apply(g.Xs, 1, sum))^2) / (6 * (n - order.s))
  # h.Xs = Xs^4 - 6 * Xs^2 + 3
  # JB2_s = sum(h.Xs) / sqrt(24 * (n - order.s) * d)
  # pval_s_JB = 1 - pchisq(JB1_s + JB2_s^2, d + 1)
  # 
  # g.Xp = Xp^3 - 3 * Xp
  # JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * (n - order.p))
  # h.Xp = Xp^4 - 6 * Xp^2 + 3
  # JB2_p = sum(h.Xp) / sqrt(24 * (n - order.p) * d)
  # pval_p_JB = 1 - pchisq(JB1_p + JB2_p^2, d + 1)
  # 
  # # demeaned version of JB Test
  # # Xs.mu = apply(Xs, 1, mean)
  # # Xs_m = Xs - Xs.mu
  # # g.Xs_m = Xs_m^3
  # # JB1_s_m = sum((apply(g.Xs_m, 1, sum))^2) / (6 * (n - order.s))
  # # h.Xs_m = Xs_m^4 - 6 * Xs_m^2 + 3
  # # JB2_s_m = sum(h.Xs_m) / sqrt(24 * (n - order.s) * d)
  # # pval_s_JB_m = 1 - pchisq(JB1_s_m + JB2_s_m^2, d + 1)
  # # 
  # # Xp.mu = apply(Xp, 1, mean)
  # # Xp_m = Xp - Xp.mu
  # # g.Xp_m = Xp_m^3
  # # JB1_p_m = sum((apply(g.Xp_m, 1, sum))^2) / (6 * (n - order.p))
  # # h.Xp_m = Xp_m^4 - 6 * Xp_m^2 + 3
  # # JB2_p_m = sum(h.Xp_m) / sqrt(24 * (n - order.p) * d)
  # # pval_p_JB_m = 1 - pchisq(JB1_p_m + JB2_p_m^2, d + 1)
  ################### KSD Test ###################
  # Observed Test Statistics
  Ts_obs = KSD_new(x = t(Xs), score_function = "gaussian", width = -1)
  Ts_obs = as.numeric(Ts_obs)
  Tp_obs = KSD_new(x = t(Xp), score_function = "gaussian", width = -1)
  Tp_obs = as.numeric(Tp_obs)
  HZs_obs = henze.(Xs)
  HZp_obs = henze.(Xp)
  # Simulate Null dist of T
  # bootssample_s = NULL
  # bootssample_p = NULL
  # sample_HZs = NULL
  # sample_HZp = NULL
  for (b in 1:nB){
    # Generate residual under H0
    # Recover Ys_b (n x d)
    Ys_b = t(sim.VAR.par0('N', coef.s, mu_hats, Rs, order.s))
    Yp_b = t(sim.VAR.par0('N', coef.p, mu_hatp, Rp, order.p))
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
    Ts_b = KSD_new(x = t(Xs_b), score_function = "gaussian", width = -1)
    Ts_b = as.numeric(Ts_b)
    Tp_b = KSD_new(x = t(Xp_b), score_function = "gaussian", width = -1)
    Tp_b = as.numeric(Tp_b)
    # bootssample_s = c(bootssample_s, Ts_b)
    # bootssample_p = c(bootssample_p, Tp_b)
    
    HZs_b = henze.(Xs_b)
    HZp_b = henze.(Xp_b)
    # sample_HZs = c(sample_HZs, HZs_b)
    # sample_HZp = c(sample_HZp, HZp_b)
  }

  # pval_s_KSD = sum((bootssample_s - Ts_obs) > 0) / length(bootssample_s)
  # pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)

  # plist_size_KSD = c(plist_size_KSD, pval_s_KSD)
  # plist_power_KSD = c(plist_power_KSD, pval_p_KSD)
  # ################### Bai Test ###################
  # ### for normal VAR
  # R2_1s = Rs[2,2] - Rs[1,2]^2 / Rs[1,1]
  # X1s = Ys.res[,1] / sqrt(Rs[1,1])
  # X2s = (Ys.res[,2] - Rs[1,2] / Rs[1,1] * Ys.res[,1]) / sqrt(R2_1s)
  # U1s = pnorm(X1s, 0, 1)
  # U2s = pnorm(X2s, 0, 1)
  # 
  # U1s[which(U1s == 1)] = 0.9999999
  # U2s[which(U2s == 1)] = 0.9999999
  # 
  # R2_1p = Rp[2,2] - Rp[1,2]^2 / Rp[1,1]
  # X1p = Yp.res[,1] / sqrt(Rp[1,1])
  # X2p = (Yp.res[,2] - Rp[1,2] / Rp[1,1] * Yp.res[,1]) / sqrt(R2_1p)
  # U1p = pnorm(X1p, 0, 1)
  # U2p = pnorm(X2p, 0, 1)
  # 
  # U1p[which(U1p == 1)] = 0.9999999
  # U2p[which(U2p == 1)] = 0.9999999
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
  # g = function(r){
  #   return(c(r,
  #            dnorm(qnorm(r,0,1),0,1),
  #            dnorm(qnorm(r,0,1),0,1) * qnorm(r,0,1)))
  # }
  # dg = function(r){
  #   return(c(1,
  #            -qnorm(r,0,1),
  #            1 - qnorm(r,0,1)^2))
  # }
  # C = function(s){
  #   mat = diag(3)
  #   for (i in 1:3){
  #     for (j in 1:3){
  #       mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
  #     }
  #   }
  #   return(mat)
  # }
  # int_dgdJ. = function(s, U1, U2){
  #   ind1 = (U1 >= s)
  #   ind2 = (U2 >= s)
  #   A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 3, 1)
  #   A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 3, 1)
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
  # Sn.s = -DEoptim(W_J.s, lower = tol, upper = 1 - tol, control = ctrl)$optim$bestval
  # Sn.p = -DEoptim(W_J.p, lower = tol, upper = 1 - tol, control = ctrl)$optim$bestval
  
  #######################################
  return(c(
           # pval_s_skew, pval_p_skew,
           # pval_s_kurt, pval_p_kurt,
           # pval_s_skku, pval_p_skku,
           # result_s_DH@p.value, result_p_DH@p.value,
           # result_s_AD@p.value, result_p_AD@p.value,
           # result_s_CM@p.value, result_p_CM@p.value,
           # result_s_HZ@p.value, result_p_HZ@p.value,
           # result_s_R@p.value,  result_p_R@p.value,
           # as.numeric(levels(result_s_mardia[["multivariateNormality"]][["p value"]]))[1], as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[1],
           # as.numeric(levels(result_s_mardia[["multivariateNormality"]][["p value"]]))[2], as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[2],
           # result_s_henze[["multivariateNormality"]][["p value"]], result_p_henze[["multivariateNormality"]][["p value"]],
           # result_s_roy[["multivariateNormality"]][["p value"]], result_p_roy[["multivariateNormality"]][["p value"]],
           # result_s_doh[["multivariateNormality"]][["p value"]], result_p_doh[["multivariateNormality"]][["p value"]],
           # result_s_energy[["multivariateNormality"]][["p value"]], result_p_energy[["multivariateNormality"]][["p value"]],
           # pval_s_JB, pval_p_JB,
           # pval_s_JB_m, pval_p_JB_m,
           # pval_s_KSD, pval_p_KSD,
           # Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
            Ts_obs, Tp_obs,
            Ts_b, Tp_b,
            HZs_obs, HZp_obs,
            HZs_b, HZp_b))
  # pb$tick()
  # Sys.sleep(1 / 100)
}



n = 500
lambda = 1
nrep = 10000
nB = 1
cl <- makeCluster(30)
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
                               'sn')) %dopar% {
                                 tryCatch(func(cnt),error = function(e) return(NA))
                               }

close(pb)
stopCluster(cl)
# result
# cat("n =",n, "lambda =",lambda)
# cbind(apply(result[,1:28]<0.01, 2, sum) / nrep,
#       apply(result[,1:28]<0.05, 2, sum) / nrep,
#       apply(result[,1:28]<0.10, 2, sum) / nrep)
# apply(result[,29:34], 2, sum) / nrep
# # which(is.na(result))


CV_s_KSD = quantile(result[,3], probs = c(0.99, 0.95, 0.90))
CV_p_KSD = quantile(result[,4], probs = c(0.99, 0.95, 0.90))
c(sum((result[,1] - CV_s_KSD[1]) > 0) / nrep,
  sum((result[,1] - CV_s_KSD[2]) > 0) / nrep,
  sum((result[,1] - CV_s_KSD[3]) > 0) / nrep,
  sum((result[,2] - CV_p_KSD[1]) > 0) / nrep,
  sum((result[,2] - CV_p_KSD[2]) > 0) / nrep,
  sum((result[,2] - CV_p_KSD[3]) > 0) / nrep)
CV_s_HZ = quantile(result[,7], probs = c(0.99, 0.95, 0.90))
CV_p_HZ = quantile(result[,8], probs = c(0.99, 0.95, 0.90))
c(sum((result[,5] - CV_s_HZ[1]) > 0) / nrep,
  sum((result[,5] - CV_s_HZ[2]) > 0) / nrep,
  sum((result[,5] - CV_s_HZ[3]) > 0) / nrep,
  sum((result[,6] - CV_p_HZ[1]) > 0) / nrep,
  sum((result[,6] - CV_p_HZ[2]) > 0) / nrep,
  sum((result[,6] - CV_p_HZ[3]) > 0) / nrep)
###################################################
###################### POWER ######################
###################################################

func_p = function(cnt){
  # for (count in (1:nrep)){
  ################# Generate Data #################
  # set.seed(cnt + 12345)
  # Generate Sample Data
  # Ys = t(sim.VAR('N', coef))
  # Yp = t(sim.VAR('t', coef))
  Yp = t(sim.VAR0('N+logN', A, sigma))
  # plot(x = 1:n, y = Yp[,1])
  ################Estimation#################
  # est_s = VAR(y = data.frame(Ys), lag.max = 5, ic = "AIC", type = "none")
  # est_p = VAR(y = data.frame(Yp), lag.max = 5, ic = "AIC", type = "none")
  est_p = VAR(y = data.frame(Yp), p = VAR.p, type = "none")
  
  order.p = est_p$p[[1]]
  
  coef.p = Acoef(est_p)
  
  Yp.res = residuals(est_p)
  
  #Estimated Mean
  mu_hatp = colMeans(Yp.res)
  #Estimated Sigma
  Rp = var(Yp.res)
  # Invsqrt of R
  # if (is.na(Rs.invsqrt)) {next}
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  # if (is.na(Rp.invsqrt)) {next}
  if (is.na(Rp.invsqrt[1,1])) return(NA)
  
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
  
  ################### Skew / Kurt Test ###################
  b1_hatp = apply(Xp^3, 1, mean)
  b2_hatp = apply(Xp^4, 1, mean)
  
  lambda_p.sk = (n - order.p) * sum(b1_hatp^2) / 6
  lambda_p.ku = (n - order.p) * sum((b2_hatp - 3)^2) / 24
  
  pval_p_skew = 1 - pchisq(lambda_p.sk, d)
  pval_p_kurt = 1 - pchisq(lambda_p.ku, d)
  pval_p_skku = 1 - pchisq(lambda_p.sk + lambda_p.ku, 2 * d)
  
  ################### DH Test ###################
  result_p_DH = DH.test(data = t(Xp), qqplot = FALSE)
  ################### AD Test ###################
  # result_s_AD = AD.test(data = t(Xs), qqplot = FALSE)
  # result_p_AD = AD.test(data = t(Xp), qqplot = FALSE)
  ################### CM Test ###################
  # result_s_CM = CM.test(data = t(Xs), qqplot = FALSE)
  # result_p_CM = CM.test(data = t(Xp), qqplot = FALSE)
  ################### HZ Test ###################
  result_p_HZ = HZ.test(data = t(Xp), qqplot = FALSE)
  ################### R Test ###################
  result_p_R = R.test(data = t(Xp), qqplot = FALSE)
  ################### Mardia Test ###################
  result_p_mardia = mvn(data = t(Xp), mvnTest = "mardia")
  ################### Henze Test ###################
  result_p_henze = mvn(data = t(Xp), mvnTest = "hz")
  ################### Roy Test ###################
  result_p_roy = mvn(data = t(Xp), mvnTest = "royston")
  ################### doh Test ###################
  result_p_doh = mvn(data = t(Xp), mvnTest = "dh")
  ################### Energy Test ###################
  result_p_energy = mvn(data = t(Xp), mvnTest = "energy")
  
  # result_s_energy = mvnorm.etest(t(Xs),1000)
  # result_p_energy = mvnorm.etest(t(Xp),1000)
  ################### JB Test ###################
  g.Xp = Xp^3 - 3 * Xp
  JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * (n - order.p))
  h.Xp = Xp^4 - 6 * Xp^2 + 3
  JB2_p = sum(h.Xp) / sqrt(24 * (n - order.p) * d)
  pval_p_JB = 1 - pchisq(JB1_p + JB2_p^2, d + 1)
  
  # demeaned version of JB Test
  # Xs.mu = apply(Xs, 1, mean)
  # Xs_m = Xs - Xs.mu
  # g.Xs_m = Xs_m^3
  # JB1_s_m = sum((apply(g.Xs_m, 1, sum))^2) / (6 * (n - order.s))
  # h.Xs_m = Xs_m^4 - 6 * Xs_m^2 + 3
  # JB2_s_m = sum(h.Xs_m) / sqrt(24 * (n - order.s) * d)
  # pval_s_JB_m = 1 - pchisq(JB1_s_m + JB2_s_m^2, d + 1)
  # 
  # Xp.mu = apply(Xp, 1, mean)
  # Xp_m = Xp - Xp.mu
  # g.Xp_m = Xp_m^3
  # JB1_p_m = sum((apply(g.Xp_m, 1, sum))^2) / (6 * (n - order.p))
  # h.Xp_m = Xp_m^4 - 6 * Xp_m^2 + 3
  # JB2_p_m = sum(h.Xp_m) / sqrt(24 * (n - order.p) * d)
  # pval_p_JB_m = 1 - pchisq(JB1_p_m + JB2_p_m^2, d + 1)
  ################### KSD Test ###################
  # Observed Test Statistics
  Tp_obs = KSD_new(x = t(Xp), score_function = "gaussian", width = -1)
  Tp_obs = as.numeric(Tp_obs)
  
  # Simulate Null dist of T
  bootssample_p = NULL
  for (b in 1:nB){
    # Generate residual under H0
    # Recover Ys_b (n x d)
    Yp_b = t(sim.VAR.par0('N', coef.p, mu_hatp, Rp, order.p))
    # Estimation
    # est_s.b = VAR(y = data.frame(Ys_b), lag.max = 5, ic = "SC", type = "none")
    # est_p.b = VAR(y = data.frame(Yp_b), lag.max = 5, ic = "SC", type = "none")
    est_p.b = VAR(y = data.frame(Yp_b), p = VAR.p, type = "none")
    # Residual
    Yp_b.res = residuals(est_p.b)
    #Estimated Mean
    mu_hatp_b = colMeans(Yp_b.res)
    #Estimated Sigma
    Rp_b = var(Yp_b.res)
    # Invsqrt of R
    Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
                            warning = function(w) matrix(NA, d, d),
                            error = function(e) matrix(NA, d, d))
    if (is.na(Rp_b.invsqrt[1,1])) next
    # Standardized Residual
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
    Tp_b = KSD_new(x = t(Xp_b), score_function = "gaussian", width = -1)
    Tp_b = as.numeric(Tp_b)
    bootssample_p = c(bootssample_p, Tp_b)
  }
  
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
  
  # plist_size_KSD = c(plist_size_KSD, pval_s_KSD)
  # plist_power_KSD = c(plist_power_KSD, pval_p_KSD)
  # ################### Bai Test ###################
  ### for normal VAR
  R2_1p = Rp[2,2] - Rp[1,2]^2 / Rp[1,1]
  X1p = Yp.res[,1] / sqrt(Rp[1,1])
  X2p = (Yp.res[,2] - Rp[1,2] / Rp[1,1] * Yp.res[,1]) / sqrt(R2_1p)
  U1p = pnorm(X1p, 0, 1)
  U2p = pnorm(X2p, 0, 1)
  
  U1p[which(U1p == 1)] = 0.9999999
  U2p[which(U2p == 1)] = 0.9999999
  ################# Functions #################
  w = function(t, weight, alpha = 0, beta = 0){
    if (weight == 'legendre'){
      return(1)
    } else if (weight == 'chebyshev1'){
      return(1 / sqrt(1 - t^2))
    } else if (weight == 'chebyshev2'){
      return(sqrt(1 - t^2))
    } else if (weight == 'jacobi'){
      return((1 - t)^alpha * (1 + t)^beta)
    }
  }
  integral = function(g, a, b, weight = 'legendre', N = 50, alpha. = 0, beta. = 0){
    quad = gauss.quad(N, weight, alpha., beta.)
    w = as.array(quad$weights)
    t = as.array(quad$nodes)
    ft = apply(t, 1, function(t) g(0.5 * ((b - a) * t + a + b)) / w(t, weight, alpha = alpha., beta = beta.))
    result = (b - a) / 2 * sum(w * ft)
    return(result)
  }
  J = function(r, n = n, U1 = U1, U2 = U2){
    J2n = (sum(U1 - r <= 0) + sum(U2 - r <= 0)) / (2 * n)
    return(J2n)
  }
  g = function(r){
    return(c(r,
             dnorm(qnorm(r,0,1),0,1),
             dnorm(qnorm(r,0,1),0,1) * qnorm(r,0,1)))
  }
  dg = function(r){
    return(c(1,
             -qnorm(r,0,1),
             1 - qnorm(r,0,1)^2))
  }
  C = function(s){
    mat = diag(3)
    for (i in 1:3){
      for (j in 1:3){
        mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
      }
    }
    return(mat)
  }
  int_dgdJ. = function(s, U1, U2){
    ind1 = (U1 >= s)
    ind2 = (U2 >= s)
    A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 3, 1)
    A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 3, 1)
    B1 = apply(A1, 1, sum) / (2 * n)
    B2 = apply(A2, 1, sum) / (2 * n)
    return(B1 + B2)
  }
  W_J = function(r, U1, U2){
    A = integral(function(s) t(dg(s)) %*% solve(C(s)) %*% int_dgdJ.(s, U1, U2), 0, r)
    return(sqrt(2 * n) * (J(r, n, U1, U2) - A))
  }
  W_J.p = function(r){
    return(- abs(W_J(r, U1p, U2p)))
  }
  ctrl = DEoptim.control(itermax = 1, trace = F)
  Sn.p = -DEoptim(W_J.p, lower = tol, upper = 1 - tol, control = ctrl)$optim$bestval
  
  #######################################
  return(c(pval_p_skew,
           pval_p_kurt,
           pval_p_skku,
           result_p_DH@p.value,
           # result_s_AD@p.value, result_p_AD@p.value,
           # result_s_CM@p.value, result_p_CM@p.value,
           result_p_HZ@p.value,
           result_p_R@p.value,
           as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[1],
           as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[2],
           result_p_henze[["multivariateNormality"]][["p value"]],
           result_p_roy[["multivariateNormality"]][["p value"]],
           result_p_doh[["multivariateNormality"]][["p value"]],
           result_p_energy[["multivariateNormality"]][["p value"]],
           pval_p_JB,
           # pval_s_JB_m, pval_p_JB_m,
           pval_p_KSD,
           Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
  # pb$tick()
  # Sys.sleep(1 / 100)
}

n = 500
lambda = 1
cl <- makeCluster(30)
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
                               "vars")) %dopar% {
                                 tryCatch(func_p(cnt),error = function(e) return(NA))
                               }

close(pb)
stopCluster(cl)
# result
cat("n =",n, "lambda =",lambda)
cbind(apply(result[,1:14]<0.01, 2, sum) / nrep,
      apply(result[,1:14]<0.05, 2, sum) / nrep,
      apply(result[,1:14]<0.10, 2, sum) / nrep)
apply(result[,15:17], 2, sum) / nrep


#################### not in use ########################
########################################################
library(pbmcapply)
# library(parallel)
# res_mclapply <- mclapply(1:nrep, func, mc.cores=80) #use all the cores on this machine
# res_mclapply <- dplyr::bind_rows(res_mclapply)
set.seed(666)
# 
# A = array(0, c(d, d, VAR.p))
# Q = array(0, c(d, d, VAR.p))
# lambda = array(0, c(d, VAR.p))
# # Specify eigenvalues
# lambda[,1] = c(0.9, -0.9)
# # lambda[,2] = c(0.6, 0.8)
# # lambda[,3] = c(0.7, 0.3)
# # Get Q (eigenvectors)
# for (k in 1:VAR.p){
#   Q[,,k] = runif(d*d, 0, 1)
#   Q[,,k] = apply(Q[,,k], 2, function(v) v / norm(v))
# }
# # Get A (AR parameters)
# for (i in 1:VAR.p){
#   A[,,i] = Q[,,i] %*% diag(lambda[,i]) %*% solve(Q[,,i])
# }
# A
# apply(A,3,eigen)

# Sample Size
n = 100
lambda = 1
# Repitition Number
nrep = 1000
# Bootstrap Number
nB = 1000
result = pbmclapply(1:nrep, func, mc.cores = 62)
# typeof(result)
# bind_rows(as.data.frame(result))
df = data.frame(result)
cat("n =",n, "lambda =", lambda)
cbind(t(t(apply(df < 0.01, 1, sum) / nrep)),
t(t(apply(df < 0.05, 1, sum) / nrep)),
t(t(apply(df < 0.10, 1, sum) / nrep)))