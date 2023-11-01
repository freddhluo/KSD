# Comparsion_VAR.R
# Multivariate Normality Tests for VAR(p) model
# Time varying and stochastic mean

# library(mvnTest) # Test for normality
# library(MVN) # Test for normality
# library(energy) # eqdist.test()
library(mvtnorm) # Generate random multivariate normal dist
library(Matrix)
library(MTS)
# library(progress)
# library(KSD)
library(pryr)
library(compositions) # rlnorm.rplus
# library(foreach)
# library(doSNOW)
library(vars)
library(DEoptim)
# library(micompr) # micompr, cmpoutput()
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

# score function for t(0, (nu-2)/nu*I, nu) distribution
# Input: data x (n x d); d.f. nu; dimension d
# Output: Sq(x) (n x d)
sqx = function(x){
  return(t(apply(x, 1, function(x) -(nu + d) / (nu - 2) * 1 / (1 + 1 / (nu - 2) * sum(x^2)) * x)))
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
# tolerance
tol = .Machine$double.eps^0.25
# Gaussian quadrature weight
weight. = "chebyshev2"
# Specify VAR Parameters
A = array(c(0.3,-0.2,0.65,-0.4,
            -0.4,-0.6,0.4,0.4,
            0.5,0.1,0.1,0.5), c(d, d, VAR.p))

# A = array(0, c(d, d, VAR.p))
# Q = array(0, c(d, d, VAR.p))
# lambda = array(0, c(d, VAR.p))
# # Specify eigenvalues
# lambda[,1] = c(0.3, -0.7, 0.2, 0.6, 0.8)
# # lambda[,2] = c(0.6, 0.8, 0.1, -0.4, 0.5)
# # lambda[,3] = c(0.7, 0.3, 0.4, 0.9, -0.5)
# 
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
# sim.VAR = function(dist, coef. = coef){
#   if (dist == 'N'){
#     ita = t(rmvnorm(n + ncut, mean = mu, sigma = sigma, method = "chol"))
#   } else if (dist == 't'){
#     ita = t(rmvt(n = n + ncut, delta = mu, sigma = (nu-2)/nu * sigma,  df = nu))
#   } else if (dist == 'logN'){
#     ita = t(matrix(rlnorm.rplus(n + ncut, meanlog = mu, varlog = sigma),nrow = n + ncut))
#   }
#   Y = matrix(0, nrow = d, ncol = n + ncut) # d x n
#   for (t in (VAR.p+1):(n+ncut)){
#     sum = rep(0, d)
#     for (i in 1:VAR.p){
#       sum = sum + coef.[,,i] %*% Y[,t-i]
#     }
#     Y[,t] = sum + ita[,t]
#   }
#   return(Y[,(ncut+1):(n+ncut)])
# }

# Generate sample data for evaluation
sim.VAR0 = function(dist, A, sigma_u){
  # Simulation
  if (dist == 'N'){
    ita = rmvnorm(n + ncut, mean = rep(0,d), sigma = diag(d), method = "chol")
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
    # ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = sigma_u, method = "chol"))
  } else if (dist == 't'){
    ita = rmvt(n + ncut, delta = rep(0,2), sigma = (5 - 2) / 5 * diag(2), df = 5)
    ita = mu + msqrt(sigma_u)$mtxsqrt %*% t(ita)
    # ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * sigma_u,  df = nu))
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

# sim.VAR.par = function(dist, coef, order){
#   if (dist == 'N'){
#     ita = t(rmvnorm(n + ncut, mean = c(0, 0), sigma = diag(2), method = "chol"))
#   } else if (dist == 't'){
#     ita = t(rmvt(n = n + ncut, delta = c(0,0), sigma = (nu-2)/nu * diag(2),  df = nu))
#   }
#   Y = matrix(0, nrow = d, ncol = n + ncut) # d x n
#   for (t in (order+1):(n+ncut)){
#     sum = rep(0, d)
#     for (i in 1:order){
#       sum = sum + coef[[i]] %*% Y[,t-i]
#     }
#     Y[,t] = sum + ita[,t]
#   }
#   return(Y[,(ncut+1):(n+ncut)])
# }

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
  Ys = t(sim.VAR0('t', A, sigma))
  # Yp = t(sim.VAR('t', coef))
  Yp = t(sim.VAR0('t+N', A, sigma))
  # plot(x = 1:n, y = Yp[,1])
  ################Estimation#################
  # est_s = VAR(y = data.frame(Ys), lag.max = 5, ic = "AIC", type = "none")
  # est_p = VAR(y = data.frame(Yp), lag.max = 5, ic = "AIC", type = "none")
  
  # estimation with AR order known
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
  if (is.na(Rs.invsqrt[1,1])) return(c(NA, NA))
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  # if (is.na(Rp.invsqrt)) {next}
  if (is.na(Rp.invsqrt[1,1])) return(c(NA, NA))
  
  Xs = Rs.invsqrt %*% (t(Ys.res) - mu_hats)
  Xp = Rp.invsqrt %*% (t(Yp.res) - mu_hatp)
  
  ################### Energy Test ###################
  # # Generate residual under H0
  # i = 1
  # while (i <= 10){
  #   # Recover Ys_b (n x d)
  #   Ys_b = t(sim.VAR.par0('t', coef.s, mu_hats, Rs, order.s))
  #   Yp_b = t(sim.VAR.par0('t', coef.p, mu_hatp, Rp, order.p))
  #   ###### Estimation ######
  #   est_s.b = VAR(y = data.frame(Ys_b), p = VAR.p, type = "none")
  #   est_p.b = VAR(y = data.frame(Yp_b), p = VAR.p, type = "none")
  #   # Residual
  #   Ys_b.res = residuals(est_s.b)
  #   Yp_b.res = residuals(est_p.b)
  #   #Estimated Mean
  #   mu_hats_b = colMeans(Ys_b.res)
  #   mu_hatp_b = colMeans(Yp_b.res)
  #   #Estimated Sigma
  #   Rs_b = var(Ys_b.res)
  #   Rp_b = var(Yp_b.res)
  #   # Invsqrt of R
  #   Rs_b.invsqrt = tryCatch(msqrt(Rs_b)$invsqrt,
  #                           warning = function(w) matrix(NA, d, d),
  #                           error = function(e) matrix(NA, d, d))
  #   if (is.na(Rs_b.invsqrt[1,1])) next
  #   Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
  #                           warning = function(w) matrix(NA, d, d),
  #                           error = function(e) matrix(NA, d, d))
  #   if (is.na(Rp_b.invsqrt[1,1])) next
  #   # Standardized Residual
  #   Xs_b = Rs_b.invsqrt %*% (t(Ys_b.res) - mu_hats_b)
  #   Xp_b = Rp_b.invsqrt %*% (t(Yp_b.res) - mu_hatp_b)
  #   if (i == 1){
  #     data_s = rbind(t(Xs), t(Xs_b))
  #     data_p = rbind(t(Xp), t(Xp_b))
  #     data2_s = rbind(Ys, Ys_b)
  #     data2_p = rbind(Yp, Yp_b)
  #   } else {
  #     data_s = rbind(data_s, t(Xs_b))
  #     data_p = rbind(data_p, t(Xp_b))
  #     data2_s = rbind(data2_s, Ys_b)
  #     data2_p = rbind(data2_p, Yp_b)
  #   }
  #   i = i + 1
  # }
  # result_s_energy = eqdist.etest(data_s, c(n - order.s, (n - order.s)*10), method = "original", R = nB)
  # pval_s_energy = result_s_energy[["p.value"]]
  # result_p_energy = eqdist.etest(data_p, c(n - order.p, (n - order.p)*10), method = "original", R = nB)
  # pval_p_energy = result_p_energy[["p.value"]]
  # 
  # result_s_micompr = cmpoutput(name = "micomprs", ve_npcs = 0.9, data = data_s, obs_lvls = factor(c(rep(1,n - order.s),rep(2,(n - order.s)*10))))
  # pval_s_micompr = mean(result_s_micompr[["p.values"]][["nonparametric"]])
  # result_p_micompr = cmpoutput(name = "micomprp", ve_npcs = 0.9, data = data_p, obs_lvls = factor(c(rep(1,n - order.p),rep(2,(n - order.p)*10))))
  # pval_p_micompr = mean(result_p_micompr[["p.values"]][["nonparametric"]])
  # 
  # # method 2
  # result_s_energy2 = eqdist.etest(data2_s, c(n, 10*n), method = "original", R = nB)
  # pval_s_energy2 = result_s_energy2[["p.value"]]
  # result_p_energy2 = eqdist.etest(data2_p, c(n, 10*n), method = "original", R = nB)
  # pval_p_energy2 = result_p_energy2[["p.value"]]
  # 
  # result_s_micompr2 = cmpoutput(name = "micomprs", ve_npcs = 0.9, data = data2_s, obs_lvls = factor(c(rep(1,n),rep(2,10*n))))
  # pval_s_micompr2 = mean(result_s_micompr2[["p.values"]][["nonparametric"]])
  # result_p_micompr2 = cmpoutput(name = "micomprp", ve_npcs = 0.9, data = data2_p, obs_lvls = factor(c(rep(1,n),rep(2,10*n))))
  # pval_p_micompr2 = mean(result_p_micompr2[["p.values"]][["nonparametric"]])
  # 
  # # result_s_peacock = peacock2(t(Xs), data.H0)
  # 
  # # result_s_cramer = cramer.test(t(Xs), data.H0)
  # # pval_s_cramer = result_s_cramer[["p.value"]]
  # # result_p_cramer = cramer.test(t(Xp), data.H0)
  # # pval_p_cramer = result_p_cramer[["p.value"]]
  
  ################### Several Monte Carlo Tests ###################
  # Observed skew/kurt Test 1 Statistics
  # b1_s = sum((t(Xs) %*% Xs)^3) / (n - order.s)^2
  # b1_p = sum((t(Xp) %*% Xp)^3) / (n - order.p)^2
  # b2_s = sum((colSums(Xs^2))^2) / (n - order.s)
  # b2_p = sum((colSums(Xp^2))^2) / (n - order.p)

  # Observed skew/kurt Test 2 Statistics
  # b1_hats = apply(Xs^3, 1, mean)
  # b1_hatp = apply(Xp^3, 1, mean)
  # b2_hats = apply(Xs^4, 1, mean)
  # b2_hatp = apply(Xp^4, 1, mean)
  # lambda_s.sk = (n - order.s) * sum(b1_hats^2) / 6
  # lambda_p.sk = (n - order.p) * sum(b1_hatp^2) / 6
  # lambda_s.ku = (n - order.s) * sum((b2_hats - 3)^2) / 24
  # lambda_p.ku = (n - order.p) * sum((b2_hatp - 3)^2) / 24

  # Observed JB Test Statistics
  # g.Xs = Xs^3 - 3 * Xs
  # JB1_s = sum((apply(g.Xs, 1, sum))^2) / (6 * (n - order.s))
  # h.Xs = Xs^4 - 6 * Xs^2 + 3
  # JB2_s = sum(h.Xs) / sqrt(24 * (n - order.s) * d)
  # JB_s.obs = JB1_s + JB2_s^2
  # 
  # g.Xp = Xp^3 - 3 * Xp
  # JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * (n - order.p))
  # h.Xp = Xp^4 - 6 * Xp^2 + 3
  # JB2_p = sum(h.Xp) / sqrt(24 * (n - order.p) * d)
  # JB_p.obs = JB1_p + JB2_p^2

  # Observed KSD Test Statistics
  Ts_obs = KSD_new(x = t(Xs), score_function = sqx, width = -1)
  Ts_obs = as.numeric(Ts_obs)
  Tp_obs = KSD_new(x = t(Xp), score_function = sqx, width = -1)
  Tp_obs = as.numeric(Tp_obs)

  # Simulate Null dist of T
  # skew_s_B = NULL
  # skew_p_B = NULL
  # kurt_s_B = NULL
  # kurt_p_B = NULL
  # skew_s = NULL
  # kurt_s = NULL
  # skku_s = NULL
  # skew_p = NULL
  # kurt_p = NULL
  # skku_p = NULL
  # JB_s = NULL
  # JB_p = NULL
  bootssample_s = NULL
  bootssample_p = NULL
  for (b in 1:nB){
    ######## Generate residual under H0 ########
    # Recover Ys_b (n x d)
    Ys_b = t(sim.VAR.par0('t', coef.s, mu_hats, Rs, order.s))
    Yp_b = t(sim.VAR.par0('t', coef.p, mu_hatp, Rp, order.p))
    ###### Estimation ######
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
    #####################
    order.s.b = est_s.b$p[[1]]
    order.p.b = est_p.b$p[[1]]
    # sk / ku (1)
    # b1_s_b = sum((t(Xs_b) %*% Xs_b)^3) / (n - order.s.b)^2
    # b1_p_b = sum((t(Xp_b) %*% Xp_b)^3) / (n - order.p.b)^2
    # b2_s_b = sum((colSums(Xs_b^2))^2) / (n - order.s.b)
    # b2_p_b = sum((colSums(Xp_b^2))^2) / (n - order.p.b)

    # skew_s_B = c(skew_s_B, b1_s_b)
    # skew_p_B = c(skew_p_B, b1_p_b)
    # kurt_s_B = c(kurt_s_B, b2_s_b)
    # kurt_p_B = c(kurt_p_B, b2_p_b)

    # sk / ku (2)
    # b1_hats.b = apply(Xs_b^3, 1, mean)
    # b1_hatp.b = apply(Xp_b^3, 1, mean)
    # b2_hats.b = apply(Xs_b^4, 1, mean)
    # b2_hatp.b = apply(Xp_b^4, 1, mean)

    # lambda_s.sk_b = (n - order.s.b) * sum(b1_hats.b^2) / 6
    # lambda_p.sk_b = (n - order.p.b) * sum(b1_hatp.b^2) / 6
    # lambda_s.ku_b = (n - order.s.b) * sum((b2_hats.b - 3)^2) / 24
    # lambda_p.ku_b = (n - order.p.b) * sum((b2_hatp.b - 3)^2) / 24

    # skew_s = c(skew_s, lambda_s.sk_b)
    # kurt_s = c(kurt_s, lambda_s.ku_b)
    # skku_s = c(skku_s, lambda_s.sk_b + lambda_s.ku_b)
    # skew_p = c(skew_p, lambda_p.sk_b)
    # kurt_p = c(kurt_p, lambda_p.ku_b)
    # skku_p = c(skku_p, lambda_p.sk_b + lambda_p.ku_b)

    # JB
    # g.Xs_b = Xs_b^3 - 3 * Xs_b
    # JB1_s_b = sum((apply(g.Xs_b, 1, sum))^2) / (6 * (n - order.s.b))
    # h.Xs_b = Xs_b^4 - 6 * Xs_b^2 + 3
    # JB2_s_b = sum(h.Xs_b) / sqrt(24 * (n - order.s.b) * d)
    # JB_s_b = JB1_s_b + JB2_s_b^2

    # g.Xp_b = Xp_b^3 - 3 * Xp_b
    # JB1_p_b  = sum((apply(g.Xp_b, 1, sum))^2) / (6 * (n - order.p.b))
    # h.Xp_b = Xp_b^4 - 6 * Xp_b^2 + 3
    # JB2_p_b = sum(h.Xp_b) / sqrt(24 * (n - order.p.b) * d)
    # JB_p_b = JB1_p_b + JB2_p_b^2

    # JB_s = c(JB_s, JB_s_b)
    # JB_p = c(JB_p, JB_p_b)

    # KSD Test Statistic for bootstrap data
    Ts_b = KSD_new(x = t(Xs_b), score_function = sqx, width = -1)
    Ts_b = as.numeric(Ts_b)
    Tp_b = KSD_new(x = t(Xp_b), score_function = sqx, width = -1)
    Tp_b = as.numeric(Tp_b)

    bootssample_s = c(bootssample_s, Ts_b)
    bootssample_p = c(bootssample_p, Tp_b)
  }

  # pval_s_sk = sum((skew_s_B - b1_s) > 0) / length(skew_s_B)
  # pval_p_sk = sum((skew_p_B - b1_p) > 0) / length(skew_p_B)
  # pval_s_ku = sum((kurt_s_B - b2_s) > 0) / length(kurt_s_B)
  # pval_p_ku = sum((kurt_p_B - b2_p) > 0) / length(kurt_p_B)
  # 
  # pval_s_skew = sum((skew_s - lambda_s.sk) > 0) / length(skew_s)
  # pval_s_kurt = sum((kurt_s - lambda_s.ku) > 0) / length(kurt_s)
  # pval_s_skku = sum((skku_s - lambda_s.sk - lambda_s.ku) > 0) / length(skku_s)
  # pval_p_skew = sum((skew_p - lambda_p.sk) > 0) / length(skew_p)
  # pval_p_kurt = sum((kurt_p - lambda_p.ku) > 0) / length(kurt_p)
  # pval_p_skku = sum((skku_p - lambda_p.sk - lambda_p.ku) > 0) / length(skku_p)
  # 
  # pval_s_JB = sum((JB_s - JB_s.obs) > 0) / length(JB_s)
  # pval_p_JB = sum((JB_p - JB_p.obs) > 0) / length(JB_p)

  pval_s_KSD = sum((bootssample_s - Ts_obs) > 0) / length(bootssample_s)
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)

  ############################ Bai Test ################################
  #### for normal VAR
  # R2_1s = Rs[2,2] - Rs[1,2]^2 / Rs[1,1]
  # X1s = Ys.res[,1] / sqrt(Rs[1,1])
  # X2s = (Ys.res[,2] - Rs[1,2] / Rs[1,1] * Ys.res[,1]) / sqrt(R2_1s)
  # U1s = pnorm(X1s, df = nu)
  # U2s = pnorm(X2s, df = nu+1)
  #### for t VAR / GARCH
  Omega_hats = Rs * (nu - 2) / nu
  a_s = (nu + Ys.res[,1]^2 / Omega_hats[1,1]) / (nu + 1)
  Omega2_1s = a_s * (Omega_hats[2,2] - Omega_hats[2,1]^2 / Omega_hats[1,1])
  X1s = Ys.res[,1] / sqrt(Omega_hats[1,1])
  X2s = (Ys.res[,2] - Omega_hats[2,1] / Omega_hats[1,1] * Ys.res[,1]) / sqrt(Omega2_1s)
  U1s = pt(X1s, df = nu)
  U2s = pt(X2s, df = nu+1)

  Omega_hatp = Rp * (nu - 2) / nu
  a_p = (nu + Yp.res[,1]^2 / Omega_hatp[1,1]) / (nu + 1)
  Omega2_1p = a_p * (Omega_hatp[2,2] - Omega_hatp[2,1]^2 / Omega_hatp[1,1])
  X1p = Yp.res[,1] / sqrt(Omega_hatp[1,1])
  X2p = (Yp.res[,2] - Omega_hatp[2,1] / Omega_hatp[1,1] * Yp.res[,1]) / sqrt(Omega2_1p)
  U1p = pt(X1p, df = nu)
  U2p = pt(X2p, df = nu+1)
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
  dg = function(r){
    return(c(1,
             -(nu + 1) * qt(r, nu) / (nu + (qt(r, nu))^2),
             1 - (nu + 1) * (qt(r, nu))^2 / (nu + (qt(r, nu))^2),
             - (nu + 2) * qt(r, nu + 1) / (nu + 1 + (qt(r, nu + 1))^2),
             1 - (nu + 2) * (qt(r, nu + 1))^2 / (nu + 1 + (qt(r, nu + 1))^2)))
  }
  C = function(s){
    mat = diag(5)
    for (i in 1:5){
      for (j in 1:5){
        mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
      }
    }
    return(mat)
  }
  int_dgdJ. = function(s, U1, U2){
    ind1 = (U1 >= s)
    ind2 = (U2 >= s)
    A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 5, 1)
    A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 5, 1)
    B1 = apply(A1, 1, sum) / (2 * n)
    B2 = apply(A2, 1, sum) / (2 * n)
    return(B1 + B2)
  }
  W_J = function(r, U1, U2){
    A = integral(function(s) t(dg(s)) %*% solve(C(s)) %*% int_dgdJ.(s, U1, U2), 0, r)
    return(sqrt(2 * n) * (J(r, n, U1, U2) - A))
  }
  W_J.s = function(r){
    return(- abs(W_J(r, U1s, U2s)))
  }
  W_J.p = function(r){
    return(- abs(W_J(r, U1p, U2p)))
  }
  ctrl = DEoptim.control(itermax = 1, trace = F)
  Sn.s = -DEoptim(W_J.s, lower = 0, upper = 1, control = ctrl)$optim$bestval
  Sn.p = -DEoptim(W_J.p, lower = 0, upper = 1, control = ctrl)$optim$bestval
  #######################################
  return(c(
           # pval_s_energy, pval_p_energy,
           # pval_s_energy2, pval_p_energy2,
           # pval_s_micompr, pval_p_micompr,
           # pval_s_micompr2, pval_p_micompr2,
           # pval_s_skew, pval_p_skew,
           # pval_s_kurt, pval_p_kurt,
           # pval_s_skku, pval_p_skku,
           # pval_s_sk, pval_p_sk,
           # pval_s_ku, pval_p_ku,
           # pval_s_JB, pval_p_JB,
           pval_s_KSD, pval_p_KSD,
           Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
  # pb$tick()
  # Sys.sleep(1 / 100)
}


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
n = 500
# Repitition Number
nrep = 1000
# Bootstrap Number
nB = 1000
lambda = 0.2
result = pbmclapply(1:nrep, func, mc.cores = 35)
# typeof(result)
# bind_rows(as.data.frame(result))
df = data.frame(result)
cat("n = ",n, "lambda =", lambda)
cbind(t(t(apply(df[1:2,] < 0.01, 1, sum) / nrep)),
      t(t(apply(df[1:2,] < 0.05, 1, sum) / nrep)),
      t(t(apply(df[1:2,] < 0.10, 1, sum) / nrep)))
apply(df[3:8,], 1, sum) / nrep


# temp = NULL
# for (i in 1:nrep){
#   if (typeof(result[[i]]) != "character") {
#     temp = rbind(temp, result[[i]])}
#   }
# typeof(result[[998]])

# n = 500 t~N
# [,1]  [,2]  [,3]
# 1  0.009 0.051 0.107
# 2  0.000 0.000 0.000
# 3  0.008 0.050 0.111
# 4  0.000 0.000 0.000
# 5  0.008 0.049 0.109
# 6  0.000 0.000 0.000
# 7  0.010 0.052 0.112
# 8  0.000 0.000 0.000
# 9  0.009 0.057 0.114
# 10 0.000 0.000 0.000
# 11 0.007 0.050 0.111
# 12 0.000 0.000 0.000
# 13 0.012 0.046 0.092
# 14 0.857 0.997 1.000
# 0.061 0.088 0.118 0.000 0.002 0.014 


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
  Yp = t(sim.VAR0('t+N', A, sigma))
  # plot(x = 1:n, y = Yp[,1])
  ################Estimation#################
  # est_s = VAR(y = data.frame(Ys), lag.max = 5, ic = "AIC", type = "none")
  # est_p = VAR(y = data.frame(Yp), lag.max = 5, ic = "AIC", type = "none")
  
  # estimation with AR order known
  est_p = VAR(y = data.frame(Yp), p = VAR.p, type = "none")
  order.p = est_p$p[[1]]
  coef.p = Acoef(est_p)
  Yp.res = residuals(est_p)
  
  #Estimated Mean
  mu_hatp = colMeans(Yp.res)
  #Estimated Sigma
  Rp = var(Yp.res)
  # Invsqrt of R
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  # if (is.na(Rp.invsqrt)) {next}
  if (is.na(Rp.invsqrt[1,1])) return(c(NA, NA))
  
  Xp = Rp.invsqrt %*% (t(Yp.res) - mu_hatp)
  ################### Several Monte Carlo Tests ###################
  # Observed skew/kurt Test 1 Statistics
  # b1_s = sum((t(Xs) %*% Xs)^3) / (n - order.s)^2
  # b1_p = sum((t(Xp) %*% Xp)^3) / (n - order.p)^2
  # b2_s = sum((colSums(Xs^2))^2) / (n - order.s)
  # b2_p = sum((colSums(Xp^2))^2) / (n - order.p)
  
  # Observed skew/kurt Test 2 Statistics
  # b1_hats = apply(Xs^3, 1, mean)
  # b1_hatp = apply(Xp^3, 1, mean)
  # b2_hats = apply(Xs^4, 1, mean)
  # b2_hatp = apply(Xp^4, 1, mean)
  # lambda_s.sk = (n - order.s) * sum(b1_hats^2) / 6
  # lambda_p.sk = (n - order.p) * sum(b1_hatp^2) / 6
  # lambda_s.ku = (n - order.s) * sum((b2_hats - 3)^2) / 24
  # lambda_p.ku = (n - order.p) * sum((b2_hatp - 3)^2) / 24
  
  # Observed JB Test Statistics
  # g.Xs = Xs^3 - 3 * Xs
  # JB1_s = sum((apply(g.Xs, 1, sum))^2) / (6 * (n - order.s))
  # h.Xs = Xs^4 - 6 * Xs^2 + 3
  # JB2_s = sum(h.Xs) / sqrt(24 * (n - order.s) * d)
  # JB_s.obs = JB1_s + JB2_s^2
  # 
  # g.Xp = Xp^3 - 3 * Xp
  # JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * (n - order.p))
  # h.Xp = Xp^4 - 6 * Xp^2 + 3
  # JB2_p = sum(h.Xp) / sqrt(24 * (n - order.p) * d)
  # JB_p.obs = JB1_p + JB2_p^2
  
  # Observed KSD Test Statistics
  Tp_obs = KSD_new(x = t(Xp), score_function = sqx, width = -1)
  Tp_obs = as.numeric(Tp_obs)
  
  # Simulate Null dist of T
  # skew_s_B = NULL
  # skew_p_B = NULL
  # kurt_s_B = NULL
  # kurt_p_B = NULL
  # skew_s = NULL
  # kurt_s = NULL
  # skku_s = NULL
  # skew_p = NULL
  # kurt_p = NULL
  # skku_p = NULL
  # JB_s = NULL
  # JB_p = NULL
  bootssample_p = NULL
  for (b in 1:nB){
    ######## Generate residual under H0 ########
    # Recover Ys_b (n x d)
    Yp_b = t(sim.VAR.par0('t', coef.p, mu_hatp, Rp, order.p))
    ###### Estimation ######
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
    #####################
    order.p.b = est_p.b$p[[1]]
    # sk / ku (1)
    # b1_s_b = sum((t(Xs_b) %*% Xs_b)^3) / (n - order.s.b)^2
    # b1_p_b = sum((t(Xp_b) %*% Xp_b)^3) / (n - order.p.b)^2
    # b2_s_b = sum((colSums(Xs_b^2))^2) / (n - order.s.b)
    # b2_p_b = sum((colSums(Xp_b^2))^2) / (n - order.p.b)
    
    # skew_s_B = c(skew_s_B, b1_s_b)
    # skew_p_B = c(skew_p_B, b1_p_b)
    # kurt_s_B = c(kurt_s_B, b2_s_b)
    # kurt_p_B = c(kurt_p_B, b2_p_b)
    
    # sk / ku (2)
    # b1_hats.b = apply(Xs_b^3, 1, mean)
    # b1_hatp.b = apply(Xp_b^3, 1, mean)
    # b2_hats.b = apply(Xs_b^4, 1, mean)
    # b2_hatp.b = apply(Xp_b^4, 1, mean)
    
    # lambda_s.sk_b = (n - order.s.b) * sum(b1_hats.b^2) / 6
    # lambda_p.sk_b = (n - order.p.b) * sum(b1_hatp.b^2) / 6
    # lambda_s.ku_b = (n - order.s.b) * sum((b2_hats.b - 3)^2) / 24
    # lambda_p.ku_b = (n - order.p.b) * sum((b2_hatp.b - 3)^2) / 24
    
    # skew_s = c(skew_s, lambda_s.sk_b)
    # kurt_s = c(kurt_s, lambda_s.ku_b)
    # skku_s = c(skku_s, lambda_s.sk_b + lambda_s.ku_b)
    # skew_p = c(skew_p, lambda_p.sk_b)
    # kurt_p = c(kurt_p, lambda_p.ku_b)
    # skku_p = c(skku_p, lambda_p.sk_b + lambda_p.ku_b)
    
    # JB
    # g.Xs_b = Xs_b^3 - 3 * Xs_b
    # JB1_s_b = sum((apply(g.Xs_b, 1, sum))^2) / (6 * (n - order.s.b))
    # h.Xs_b = Xs_b^4 - 6 * Xs_b^2 + 3
    # JB2_s_b = sum(h.Xs_b) / sqrt(24 * (n - order.s.b) * d)
    # JB_s_b = JB1_s_b + JB2_s_b^2
    
    # g.Xp_b = Xp_b^3 - 3 * Xp_b
    # JB1_p_b  = sum((apply(g.Xp_b, 1, sum))^2) / (6 * (n - order.p.b))
    # h.Xp_b = Xp_b^4 - 6 * Xp_b^2 + 3
    # JB2_p_b = sum(h.Xp_b) / sqrt(24 * (n - order.p.b) * d)
    # JB_p_b = JB1_p_b + JB2_p_b^2
    
    # JB_s = c(JB_s, JB_s_b)
    # JB_p = c(JB_p, JB_p_b)
    
    # KSD Test Statistic for bootstrap data
    Tp_b = KSD_new(x = t(Xp_b), score_function = sqx, width = -1)
    Tp_b = as.numeric(Tp_b)
    
    bootssample_p = c(bootssample_p, Tp_b)
  }
  
  # pval_s_sk = sum((skew_s_B - b1_s) > 0) / length(skew_s_B)
  # pval_p_sk = sum((skew_p_B - b1_p) > 0) / length(skew_p_B)
  # pval_s_ku = sum((kurt_s_B - b2_s) > 0) / length(kurt_s_B)
  # pval_p_ku = sum((kurt_p_B - b2_p) > 0) / length(kurt_p_B)
  # 
  # pval_s_skew = sum((skew_s - lambda_s.sk) > 0) / length(skew_s)
  # pval_s_kurt = sum((kurt_s - lambda_s.ku) > 0) / length(kurt_s)
  # pval_s_skku = sum((skku_s - lambda_s.sk - lambda_s.ku) > 0) / length(skku_s)
  # pval_p_skew = sum((skew_p - lambda_p.sk) > 0) / length(skew_p)
  # pval_p_kurt = sum((kurt_p - lambda_p.ku) > 0) / length(kurt_p)
  # pval_p_skku = sum((skku_p - lambda_p.sk - lambda_p.ku) > 0) / length(skku_p)
  # 
  # pval_s_JB = sum((JB_s - JB_s.obs) > 0) / length(JB_s)
  # pval_p_JB = sum((JB_p - JB_p.obs) > 0) / length(JB_p)
  
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
  
  ############################ Bai Test ################################
  #### for normal VAR
  # R2_1s = Rs[2,2] - Rs[1,2]^2 / Rs[1,1]
  # X1s = Ys.res[,1] / sqrt(Rs[1,1])
  # X2s = (Ys.res[,2] - Rs[1,2] / Rs[1,1] * Ys.res[,1]) / sqrt(R2_1s)
  # U1s = pnorm(X1s, df = nu)
  # U2s = pnorm(X2s, df = nu+1)
  #### for t VAR / GARCH
  Omega_hatp = Rp * (nu - 2) / nu
  a_p = (nu + Yp.res[,1]^2 / Omega_hatp[1,1]) / (nu + 1)
  Omega2_1p = a_p * (Omega_hatp[2,2] - Omega_hatp[2,1]^2 / Omega_hatp[1,1])
  X1p = Yp.res[,1] / sqrt(Omega_hatp[1,1])
  X2p = (Yp.res[,2] - Omega_hatp[2,1] / Omega_hatp[1,1] * Yp.res[,1]) / sqrt(Omega2_1p)
  U1p = pt(X1p, df = nu)
  U2p = pt(X2p, df = nu+1)
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
  dg = function(r){
    return(c(1,
             -(nu + 1) * qt(r, nu) / (nu + (qt(r, nu))^2),
             1 - (nu + 1) * (qt(r, nu))^2 / (nu + (qt(r, nu))^2),
             - (nu + 2) * qt(r, nu + 1) / (nu + 1 + (qt(r, nu + 1))^2),
             1 - (nu + 2) * (qt(r, nu + 1))^2 / (nu + 1 + (qt(r, nu + 1))^2)))
  }
  C = function(s){
    mat = diag(5)
    for (i in 1:5){
      for (j in 1:5){
        mat[i,j] = integral(function(r) (dg(r) %o% dg(r))[i,j], s, 1, weight = weight., alpha = alpha, beta = beta)
      }
    }
    return(mat)
  }
  int_dgdJ. = function(s, U1, U2){
    ind1 = (U1 >= s)
    ind2 = (U2 >= s)
    A1 = apply(as.array(U1), 1, dg) * repmat(t(ind1), 5, 1)
    A2 = apply(as.array(U2), 1, dg) * repmat(t(ind2), 5, 1)
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
  Sn.p = -DEoptim(W_J.p, lower = 0, upper = 1, control = ctrl)$optim$bestval
  #######################################
  return(c(
    # pval_s_energy, pval_p_energy,
    # pval_s_energy2, pval_p_energy2,
    # pval_s_micompr, pval_p_micompr,
    # pval_s_micompr2, pval_p_micompr2,
    # pval_s_skew, pval_p_skew,
    # pval_s_kurt, pval_p_kurt,
    # pval_s_skku, pval_p_skku,
    # pval_s_sk, pval_p_sk,
    # pval_s_ku, pval_p_ku,
    # pval_s_JB, pval_p_JB,
    pval_p_KSD,
    Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
  # pb$tick()
  # Sys.sleep(1 / 100)
}

set.seed(666)
# Repitition Number
nrep = 1000
# Bootstrap Number
nB = 1000
# Sample Size
n = 500
lambda = 1
result = pbmclapply(1:nrep, func_p, mc.cores = 30)
# typeof(result)
# bind_rows(as.data.frame(result))
df = data.frame(result)
cat("n = ",n, "lambda =", lambda)
cbind(t(t(apply(df[1,] < 0.01, 1, sum) / nrep)),
      t(t(apply(df[1,] < 0.05, 1, sum) / nrep)),
      t(t(apply(df[1,] < 0.10, 1, sum) / nrep)))
apply(df[2:4,], 1, sum) / nrep

