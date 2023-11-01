# Test for normal Copula
# for ARMA(1,1)-normal distribution
# dimension = 2

# library(mvnTest)
library(mvtnorm)
library(Matrix)
library(MTS) # msqrt, VAR
# library(progress) # progress bar
# library(KSD)
library(pryr)
library(compositions) # rlnorm.rplus
library(MVN)
# library(energy) # eqdist.test()
library(fMultivar) # rcauchy2d
# library(cramer) # cramer.test()
# library(micompr) # micompr, cmpoutput()
# library(Peacock.test) # peacock2()
library(statmod) # gauss.quad() used
library(DEoptim)
library(foreach)
library(doSNOW)

library(copula)
library(tseries)
library(QRM)
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
# Revised version of find_median_distance
find_median_distance <- function(Z){
  
  if(is.data.frame(Z)){
    Z = data.matrix(Z)
  }else{
    Z = as.array(Z)
  }
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]
  
  # if size of Z is greater than 100, randomly sample 100 points
  # if(size1 > 500){
  #   if(is.na(size2)){
  #     Zmed <- Z[sample(size1,500)]
  #   }else{
  #     Zmed <- Z[sample(size1,500),]
  #   }
  #   size1 = 500
  # }else{
  #   Zmed <- Z
  # }
  Zmed <- Z
  
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
#   score_function: "gaussian" or user-defined
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

# Function: KSD_Nres
# Input: data Y (n x d); KSD_new; nB; n
# Output: p-value or [test statistic (_obs and _b) for warp-speed method]
KSD_Nres = function(Y, KSD_new. = KSD_new, nB. = nB, n. = n){
  ################ Estimation #################
  est1 = arma(x = Y[,1], order = c(1, 1))
  est2 = arma(x = Y[,2], order = c(1, 1))
  X = rbind(est1$residuals[-1], est2$residuals[-1])
  sigma1_hat = sd(X[1,])
  sigma2_hat = sd(X[2,])
  X = rbind(est1$residuals[-1] / sigma1_hat,
            est2$residuals[-1] / sigma2_hat)
  n0 = n. - 1
  coef1 = c(est1$coef, sigma1_hat)
  coef2 = c(est2$coef, sigma2_hat)
  
  rho_hat = rho
  
  sqx = function(x){
    rho = rho_hat
    R = matrix(c(1, rho,
                 rho, 1), 2, 2)
    return(-(solve(R) ) %*% x)
  }
  
  ################### KSD Test ###################
  # Observed Test Statistics
  T_obs = KSD_new.(x = t(X), score_function = t(apply(t(X), 1, sqx)), width = -1)
  T_obs = as.numeric(T_obs)
  
  # Distribution of T
  # sample = NULL
  for (b in 1:nB.){
    # Generate residual under H0
    Y_b = DGP_b(coef1, coef2, rho_hat)
    
    est1_b = arma(x = Y_b[,1], order = c(1, 1))
    est2_b = arma(x = Y_b[,2], order = c(1, 1))
    X_b = rbind(est1_b$residuals[-1], est2_b$residuals[-1])
    sigma1_hat_b = sd(X_b[1,])
    sigma2_hat_b = sd(X_b[2,])
    X_b = rbind(est1_b$residuals[-1] / sigma1_hat_b,
              est2_b$residuals[-1] / sigma2_hat_b)
    
    rho_hat_b = rho_hat
    sqx_b = function(x){
      rho = rho_hat_b
      R = matrix(c(1, rho,
                   rho, 1), 2, 2)
      return(-(solve(R) ) %*% x)
    }
    # KSD Test Statistic for bootstrap data
    T_b = KSD_new.(x = t(X_b), score_function = t(apply(t(X_b), 1, sqx_b)), width = -1)
    T_b = as.numeric(T_b)
    # sample = c(sample, T_b)
  }
  
  # pval_KSD = sum((sample - T_obs) > 0) / nB.
  return(c(T_obs,T_b))
}


###############   Settings   ###############
# Dimension of data
d = 2
# Dimension of theta
# p = d + d^2
# Sample Size
n = 500
# Repitition Number
nrep = 10000
# Bootstrap Number
nB = 1

# Mean
# mu = rep(0,d)
# mu = c(0.5, 0.3)
# Sample Size & Block length
# Spec_Nsample = function(N, w = 0.3){
#   if (N == 100){
#     n <<- 100
#     m1 <<- n * w
#     m2 <<- n * (1 - w)
#   } else if (N == 200){
#     n <<- 200
#     m1 <<- n * w
#     m2 <<- n * (1 - w)
#   } else if (N == 500){
#     n <<- 500
#     m1 <<- n * w
#     m2 <<- n * (1 - w)
#   }
# }
# Spec_Nsample(100, 0.3)
# weight
# w = 0.3
# Variance-Covariance Matrix
# sigma = matrix(c(1,0.5,0.5,1), ncol = 2)

# d.f.
nu = 5
# nu_1 = 4
# nu_2 = 6

# ARMA parameters
delta0 = c(0.3, 0.4)
alpha0 = c(0.4, 0.8)
beta0  = c(0.5, 0.7)
sigma0 = c(0.2, 0.6)
ncut = 1000
tau = 0.3
# tau = 0.7
sigma = matrix(c(1, sin(tau * pi / 2),
                 sin(tau * pi / 2), 1), 2, 2)
rho = sin(tau * pi / 2)
DGP = function(cop){
  if (cop == 't'){
    cop = rcopula.t(n + ncut, df = nu, Sigma = sigma)
  } 
  else if (cop == 'N'){
    cop = rcopula.gauss(n + ncut, sigma)
  } else if (cop == 'Clayton'){
    cop = rcopula.clayton(n = n + ncut, theta = 2*tau/(1-tau), d = 2)
  } else if (cop == 'Gumbel'){
    cop = rcopula.gumbel(n = n + ncut, theta = 1 / (1-tau), d = 2)
  }
  # eps1 = sqrt(2) / sigma0[1] * rt(n = n + ncut, df = nu_1)
  # eps2 = sqrt(1.5) / sigma0[2] * rt(n = n + ncut, df = nu_2)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  for (t in 2:(n + ncut)){
    Y1[t] = delta0[1] + alpha0[1] * Y1[t-1] + 
      sigma0[1] * qnorm(p = cop[t,1]) +
      beta0[1] * sigma0[1] * qnorm(p = cop[t-1,1])
    Y2[t] = delta0[2] + alpha0[2] * Y2[t-1] + 
      sigma0[2] * qnorm(p = cop[t,2]) + 
      beta0[2] * sigma0[2] * qnorm(p = cop[t-1,2])
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  
  return(Y)
}
DGP_b = function(coef1, coef2, rho){
  sigma = matrix(c(1, rho, rho, 1), 2, 2)
  cop = rcopula.gauss(n + ncut, Sigma = sigma)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  delta = c(coef1[3], coef2[3])
  alpha = c(coef1[1], coef2[1])
  beta  = c(coef1[2], coef2[2])
  sigmahat = c(coef1[4], coef2[4])
  for (t in 2:(n + ncut)){
    Y1[t] = delta[1] + alpha[1] * Y1[t-1] + 
      sigmahat[1] * qnorm(p = cop[t,1]) + 
      beta[1] * sigmahat[1] * qnorm(p = cop[t-1,1])
    Y2[t] = delta[2] + alpha[2] * Y2[t-1] + 
      sigmahat[2] * qnorm(p = cop[t,2]) + 
      beta[2] * sigmahat[2] * qnorm(p = cop[t-1,2])
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  
  return(Y)
}
################### Execution ###################
func = function(cnt){
  ################# Generate Data #################
  # set.seed(count + 12345)
  # Generate Sample Data
  
  Ys = DGP('N')
  # colMeans(Ys)
  # cor(Ys)
  # cor(Ys, method = 'kendall')
  
  Yp = DGP('t')
  # cor(Yp)
  # cor(Yp, method = 'kendall')
  ############## Estimation#################
  # est1s = arma(x = Ys[,1], order = c(1, 1))
  # est2s = arma(x = Ys[,2], order = c(1, 1))
  # Xs = rbind(est1s$residuals[-1], est2s$residuals[-1])
  # sigma1_hats = sd(Xs[1,])
  # sigma2_hats = sd(Xs[2,])
  # n0 = n - 1
  # coef1s = c(est1s$coef, sigma1_hats)
  # coef2s = c(est2s$coef, sigma2_hats)
  # 
  # F1hats = function(x){
  #   return(1 / (n0 + 1) * sum(Xs[1,] <= x))
  # }
  # F2hats = function(x){
  #   return(1 / (n0 + 1) * sum(Xs[2,] <= x))
  # }
  # Fn1s = apply(as.array(Xs[1,]), 1, F1hats)
  # Fn2s = apply(as.array(Xs[2,]), 1, F2hats)
  # Fns = cbind(Fn1s, Fn2s)
  # 
  # # fitn = fitCopula(tCopula(dim = 2,dispstr = "ex", df = 5, df.fixed = T), Fn, method = 'ml')
  # # rho_hat = fitn@estimate[1]
  # objs = function(rho){
  #   sigma = matrix(c(1, rho, rho, 1), 2, 2)
  #   sum = 0
  #   for (t in 1:n0){
  #     x = c(qt(Fns[t,1],nu), qt(Fns[t,2],nu))
  #     sum = sum - 1 / 2 * log(det(sigma)) - (nu + d) / 2 *
  #       log(1 + t(x) %*% solve(sigma) %*% x / nu)
  #   }
  #   return(as.numeric(- sum / n0))
  # }
  # 
  # opts = optimize(f = objs, lower = -1 + 0.00001, upper = 1 - 0.00001, tol = 0.00001)
  # rho_hats = opts$minimum
  # 
  # sqxs = function(x){
  #   rho = rho_hats
  #   x1 = qt(pt(x[1], nu_1), nu)
  #   x2 = qt(pt(x[2], nu_2), nu)
  #   sqx1 = (- (nu + d) * (x1 - rho * x2) / 
  #             (nu * (1 - rho^2) + x1^2 - 2 * rho * x1 * x2 + x2^2) +
  #             (nu + 1) * x1 / (nu + x1^2))
  #   sqx2 = (- (nu + d) * (x2 - rho * x1) / 
  #             (nu * (1 - rho^2) + x1^2 - 2 * rho * x1 * x2 + x2^2) +
  #             (nu + 1) * x2 / (nu + x2^2))
  #   sqx1 = sqx1 * dt(x[1], nu_1) / dt(qt(pt(x[1], nu_1), nu), nu)
  #   sqx2 = sqx2 * dt(x[2], nu_2) / dt(qt(pt(x[2], nu_2), nu), nu)
  #   return(c(sqx1, sqx2))
  # }
  # 
  # 
  # est1p = arma(x = Yp[,1], order = c(1, 1))
  # est2p = arma(x = Yp[,2], order = c(1, 1))
  # Xp = rbind(est1p$residuals[-1], est2p$residuals[-1])
  # sigma1_hatp = sd(Xp[1,])
  # sigma2_hatp = sd(Xp[2,])
  # coef1p = c(est1p$coef, sigma1_hatp)
  # coef2p = c(est2p$coef, sigma2_hatp)
  # 
  # F1hatp = function(x){
  #   return(1 / (n0 + 1) * sum(Xp[1,] <= x))
  # }
  # F2hatp = function(x){
  #   return(1 / (n0 + 1) * sum(Xp[2,] <= x))
  # }
  # Fn1p = apply(as.array(Xp[1,]), 1, F1hatp)
  # Fn2p = apply(as.array(Xp[2,]), 1, F2hatp)
  # Fnp = cbind(Fn1p, Fn2p)
  # 
  # # fitn = fitCopula(tCopula(dim = 2,dispstr = "ex", df = 5, df.fixed = T), Fn, method = 'ml')
  # # rho_hat = fitn@estimate[1]
  # objp = function(rho){
  #   sigma = matrix(c(1, rho, rho, 1), 2, 2)
  #   sum = 0
  #   for (t in 1:n0){
  #     x = c(qt(Fnp[t,1],nu), qt(Fnp[t,2],nu))
  #     sum = sum - 1 / 2 * log(det(sigma)) - (nu + d) / 2 *
  #       log(1 + t(x) %*% solve(sigma) %*% x / nu)
  #   }
  #   return(as.numeric(- sum / n0))
  # }
  # 
  # optp = optimize(f = objp, lower = -1 + 0.00001, upper = 1 - 0.00001, tol = 0.00001)
  # rho_hatp = optp$minimum
  # 
  # sqxp = function(x){
  #   rho = rho_hatp
  #   x1 = qt(pt(x[1], nu_1), nu)
  #   x2 = qt(pt(x[2], nu_2), nu)
  #   sqx1 = (- (nu + d) * (x1 - rho * x2) / 
  #             (nu * (1 - rho^2) + x1^2 - 2 * rho * x1 * x2 + x2^2) +
  #             (nu + 1) * x1 / (nu + x1^2))
  #   sqx2 = (- (nu + d) * (x2 - rho * x1) / 
  #             (nu * (1 - rho^2) + x1^2 - 2 * rho * x1 * x2 + x2^2) +
  #             (nu + 1) * x2 / (nu + x2^2))
  #   sqx1 = sqx1 * dt(x[1], nu_1) / dt(qt(pt(x[1], nu_1), nu), nu)
  #   sqx2 = sqx2 * dt(x[2], nu_2) / dt(qt(pt(x[2], nu_2), nu), nu)
  #   return(c(sqx1, sqx2))
  # }
  # ############## KSD Block Test #############
  # # Block
  # X1s = Xs[,1:m1]
  # X2s = Xs[,(n0-m2+1):n0]
  # 
  # X1p = Xp[,1:m1]
  # X2p = Xp[,(n0-m2+1):n0]
  # 
  # # KSD Package Results
  # result1s <- KSD(t(X1s),score_function = t(apply(t(X1s), 1, sqxs)), 'rbf',-1.0)
  # result2s <- KSD(t(X2s),score_function = t(apply(t(X2s), 1, sqxs)), 'rbf',-1.0)
  # Ts1obs = result1s[["ksd"]]
  # Ts2obs = result2s[["ksd"]]
  # Tsobs = w * Ts1obs + (1-w) * Ts2obs
  # bootssample_s = w * result1s[["bootStrapSamples"]] + (1-w) * result2s[["bootStrapSamples"]]
  # pval_s_KSD_b = sum((bootssample_s - Tsobs) > 0) / 1000 # nboots=1000
  # 
  # result1p <- KSD(t(X1p),score_function = t(apply(t(X1p), 1, sqxp)),'rbf',-1.0)
  # result2p <- KSD(t(X2p),score_function = t(apply(t(X2p), 1, sqxp)),'rbf',-1.0)
  # Tp1obs = result1p[["ksd"]]
  # Tp2obs = result2p[["ksd"]]
  # Tpobs = w * Tp1obs + (1-w) * Tp2obs
  # bootssample_p = w * result1p[["bootStrapSamples"]] + (1-w) * result2p[["bootStrapSamples"]]
  # pval_p_KSD_b = sum((bootssample_p - Tpobs) > 0) / 1000
  ################### KSD Nres Test ###################
  KSD_s = KSD_Nres(Ys)
  Ts_obs = KSD_s[1]
  Ts_b = KSD_s[2]
  KSD_p = KSD_Nres(Yp) 
  Tp_obs = KSD_p[1]
  Tp_b = KSD_p[2]
  # pval_s_KSD. = KSD_Nres.(Ys)
  # pval_p_KSD. = KSD_Nres.(Yp)
  ################### return ###################
  return(c(
    # (pval_s_KSD < 0.005) | (pval_s_KSD > 0.995),
    # (pval_s_KSD < 0.025) | (pval_s_KSD > 0.975),
    # (pval_s_KSD < 0.05)  | (pval_s_KSD > 0.95),
    # (pval_p_KSD < 0.005) | (pval_p_KSD > 0.995),
    # (pval_p_KSD < 0.025) | (pval_p_KSD > 0.975),
    # (pval_p_KSD < 0.05)  | (pval_p_KSD > 0.95),
    # pval_s_KSD < 0.01,
    # pval_s_KSD < 0.05,
    # pval_s_KSD < 0.10,
    # pval_p_KSD < 0.01,
    # pval_p_KSD < 0.05,
    # pval_p_KSD < 0.10))
    Ts_obs, Tp_obs,
    Ts_b, Tp_b))
  
  # pval_s_KSD. < 0.01,
  # pval_s_KSD. < 0.05,
  # pval_s_KSD. < 0.10,
  # pval_p_KSD. < 0.01,
  # pval_p_KSD. < 0.05,
  # pval_p_KSD. < 0.10))
}


################### Size Evaluation ###################
library(pbmcapply)
# set.seed(666)
# Sample Size
n = 500
# Repitition Number
nrep = 10000
# Bootstrap Number
nB = 1
# lambda = 0.2
result = pbmclapply(1:nrep, func, mc.cores = 4)
# typeof(result)
# bind_rows(as.data.frame(result))
res = NULL
for (i in 1:nrep){
  res = rbind(res, result[[i]])
}
# cat("n =", n, "lambda =", lambda)
CV_s = quantile(res[,3], probs = c(0.99, 0.95, 0.90))
CV_p = quantile(res[,4], probs = c(0.99, 0.95, 0.90))
c(sum((res[,1] - CV_s[1]) > 0) / nrep,
  sum((res[,1] - CV_s[2]) > 0) / nrep,
  sum((res[,1] - CV_s[3]) > 0) / nrep)
c(sum((res[,2] - CV_p[1]) > 0) / nrep,
  sum((res[,2] - CV_p[2]) > 0) / nrep,
  sum((res[,2] - CV_p[3]) > 0) / nrep)

# # ---------- Use foreach and doSNOW package ---------- #
# n = 500
# # n = 50
# cl <- makeCluster(3)
# registerDoSNOW(cl)
# pb = txtProgressBar(max = nrep, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts = list(progress = progress)
# result = foreach(cnt = 1:nrep,
#                  .combine = 'rbind',
#                  .options.snow = opts,
#                  .packages = c("mvtnorm",
#                                "Matrix",
#                                "MTS",
#                                "progress",
#                                "pryr",
#                                "compositions",
#                                "MVN",
#                                # "energy",
#                                'mvnTest',
#                                'fMultivar',
#                                'DEoptim',
#                                'statmod',
#                                'copula',
#                                'tseries',
#                                'QRM')) %dopar% func(cnt)
# 
# 
# close(pb)
# stopCluster(cl)
# # result
# # cbind(apply(result < 0.01, 2, sum) / nrep,
# #       apply(result < 0.05, 2, sum) / nrep,
# #       apply(result < 0.10, 2, sum) / nrep)
# nNA = sum(is.na(result[,1]))
# res = na.omit(result)
# 
# apply(res, 2, sum) / length(res[,1])

