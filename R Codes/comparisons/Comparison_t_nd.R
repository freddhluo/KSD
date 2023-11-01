# Comparison of tests
# for Constant Mean Constant Variance t distribution
# dimension = d

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
library(foreach)
library(doSNOW)

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

# g(theta, Y)
# Input: theta (1 x p), Y: original data matrix (d x n)
# Output: value of g (d x n)
g = function(theta, Y){
  mu = repmat(theta[1:d], 1, n)
  # sigma_sqrtinv = matrix(theta[3:6],nrow = 2, ncol = 2)
  sigma = matrix(theta[(d+1):(d+d^2)], nrow = d, ncol = d)
  sigma_sqrtinv = msqrt(sigma)$invsqrt
  return(sigma_sqrtinv %*% (Y - mu))
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

# Function: KSD_Nres
# Input: data Y (n x d); estimation function g; KSD_new; nB; n
# Output: p-value
KSD_Nres = function(Y, g. = g, KSD_new. = KSD_new, nB. = nB, n. = n){
  ################ Estimation #################
  #Estimated Mean
  mu_hat = colMeans(Y)
  #Estimated Sigma
  sigma_hat = var(Y)
  #Standardize to get residual (d x T)
  theta_hat = c(mu_hat, c(sigma_hat))
  X = g.(theta_hat, t(Y))
  
  ################### KSD Test ###################
  # Observed Test Statistics
  T_obs = KSD_new.(x = t(X), score_function = sqx, width = -1)
  T_obs = as.numeric(T_obs)
  
  # Distribution of T
  sample = NULL
  for (b in 1:nB.){
    # Generate residual under H0
    # set.seed(b+1)
    ita.N = rmvt(n., delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)
    # Recover Ys_b (n x d)
    Y_b = t(mu_hat + msqrt(sigma_hat)$mtxsqrt %*% t(ita.N))
    # Estimation
    mu_hat_b = colMeans(Y_b)
    sigma_hat_b = var(Y_b)
    
    theta_hat_b = c(mu_hat_b, c(sigma_hat_b))
    X_b = g.(theta_hat_b, t(Y_b))
    # KSD Test Statistic for bootstrap data
    T_b = KSD_new.(x = t(X_b), score_function = sqx, width = -1)
    T_b = as.numeric(T_b)
    sample = c(sample, T_b)
  }
  
  pval_KSD = sum((sample - T_obs) > 0) / nB.
  return(pval_KSD)
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
# Dimension of theta
p = d + d^2
# Sample Size
n = 100
# Repitition Number
nrep = 1000
# Bootstrap Number
nB = 1000
# Mean
mu = rep(0,d)
# mu_log = rep(exp(1/2), 2)

# Variance-Covariance Matrix
sigma = matrix(c(1,0.5,0.5,1), ncol = 2)
# sigma = matrix(c(1, 0.5, 0.25,
#                  0.5, 1, 0.5,
#                  0.25, 0.5, 1), 3, 3)
# sigma = matrix(c(1, 0.5, 0.25, 0.125, 0.0625,
#                  0.5, 1, 0.5, 0.25, 0.125,
#                  0.25, 0.5, 1, 0.5, 0.25,
#                  0.125, 0.25, 0.5, 1, 0.5,
#                  0.0625, 0.125, 0.25, 0.5, 1), d, d)
# sigma_log = matrix(c(exp(2)-exp(1), exp(1.5)-exp(1),
#                      exp(1.5)-exp(1), exp(2)-exp(1)),
#                    ncol = 2)

# d.f.
nu = 5
# tolerance
tol = .Machine$double.eps^0.25
# Gaussian quadrature weight
# weight. = "chebyshev2"

lambda = 1
nu = 5
cosi = c(1, 1.3)

# gamma_hat = c(0, 0.2, -0.2, 0, -0.1)
gamma_hat = c(0, -0.6)
CP = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat)
DP = cp2dp(CP, "SN")


# d = 2
CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 17.22322)
DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# d = 5
# CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 70.55935)
# DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu

################### Execution ###################
func = function(cnt){
  # for (count in (1:nrep)){
  ################# Generate Data #################
  # set.seed(count + 12345)
  # Generate Sample Data
  # Ys = rmvt(n, delta = mu, sigma = sigma, df = nu)
  Ys = rmvt(n = n, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)
  Ys = t(mu + msqrt(sigma)$mtxsqrt %*% t(Ys))
  
  
  # Skew normal distibution
  Yp = rmsn(n, dp = DP)
  Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # Skew t distibution
  # Yp = rmst(n, dp = DP_ST)
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # Yp = rcauchy2d(n = n, rho = 0.5)
  # Yp = matrix(rlnorm.rplus(n, meanlog = mu, varlog = sigma),nrow = n)
  # Yp = rmvt(n, delta = mu, sigma = sigma, df = 7)
  # Yp = rmvt(n, delta = mu, sigma = sigma, df = 9)
  # Yp = rmvnorm(n, mean = mu, sigma = sigma)
  
  # normal
  # Yp = rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol")
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Ys))
  
  # Skew-student t
  # temp = rmvt(n = n, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)
  # Yp = NULL
  # for (i in 1:d){
  #   Wi = rbinom(n, 1, cosi[i]^2 / (1 + cosi[i]^2))
  #   tempi = Wi * abs(temp[,i]) * cosi[i] - 
  #     (1 - Wi) * abs(temp[,i]) / cosi[i]
  #   mi = gamma((nu - 1) / 2) * sqrt(nu - 2) / sqrt(pi) / gamma(nu / 2) * 
  #     (cosi[i] - 1 / cosi[i])
  #   si2 = (cosi[i]^2 + 1 / cosi[i]^2 - 1) - mi^2
  #   tempi = (tempi - mi) / sqrt(si2)
  #   Yp = rbind(Yp, tempi)
  # }
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% Yp)
  ################ Estimation#################
  #Estimated Mean
  mu_hats = colMeans(Ys)
  mu_hatp = colMeans(Yp)
  #Estimated Sigma
  sigma_hats = var(Ys)
  sigma_hatp = var(Yp)
  #Standardize to get residual (2 x T)
  theta_hats = c(mu_hats, c(sigma_hats))
  Xs = g(theta_hats, t(Ys))
  theta_hatp = c(mu_hatp, c(sigma_hatp))
  Xp = g(theta_hatp, t(Yp))
  ################### KSD Block Test ###################
  # # Block
  # X1s = Xs[,1:m1]
  # X2s = Xs[,(n-m2+1):n]
  # 
  # X1p = Xp[,1:m1]
  # X2p = Xp[,(n-m2+1):n]
  # 
  # # KSD Package Results
  # result1s <- KSD(t(X1s),score_function= sqx, 'rbf',-1.0)
  # result2s <- KSD(t(X2s),score_function= sqx, 'rbf',-1.0)
  # 
  # result1p <- KSD(t(X1p),score_function = sqx,'rbf',-1.0)
  # result2p <- KSD(t(X2p),score_function = sqx,'rbf',-1.0)
  # 
  # Ts1obs = result1s[["ksd"]]
  # Ts2obs = result2s[["ksd"]]
  # 
  # Tp1obs = result1p[["ksd"]]
  # Tp2obs = result2p[["ksd"]]
  # 
  # Tsobs = w * Ts1obs + (1-w) * Ts2obs
  # Tpobs = w * Tp1obs + (1-w) * Tp2obs
  # 
  # bootssample_s = w * result1s[["bootStrapSamples"]] + (1-w) * result2s[["bootStrapSamples"]]
  # bootssample_p = w * result1p[["bootStrapSamples"]] + (1-w) * result2p[["bootStrapSamples"]]
  # 
  # pval_s_KSD_b = sum((bootssample_s - Tsobs) > 0) / 1000 # nboots=1000
  # pval_p_KSD_b = sum((bootssample_p - Tpobs) > 0) / 1000
  # #
  # # plist_size_KSD = c(plist_size_KSD, pval_s_KSD)
  # # plist_power_KSD = c(plist_power_KSD, pval_p_KSD)
  
  ################### KSD Nres Test ###################
  pval_s_KSD = KSD_Nres(Ys)
  pval_p_KSD = KSD_Nres(Yp)
  ################### Bai Test ###################
  Omega_hats = sigma_hats * (nu - 2) / nu
  mu2_1s = rep(mu_hats[2], n) + Omega_hats[2,1] / Omega_hats[1,1] * (Ys[,1] - rep(mu_hats[1],n))
  a_s = (nu + (Ys[,1] - rep(mu_hats[1],n))^2 / Omega_hats[1,1] ) / (nu + 1)
  Omega2_1s = a_s * (Omega_hats[2,2] - Omega_hats[2,1] / Omega_hats[1,1] * Omega_hats[1,2])
  X1s = (Ys[,1] - mu_hats[1]) / sqrt(Omega_hats[1,1])
  X2s = (Ys[,2] - mu2_1s) / sqrt(Omega2_1s)
  U1s = pt(X1s, df = nu)
  U2s = pt(X2s, df = nu+1)
  Us = c(U1s, U2s)
  v1s = c(0, sort(U1s), 1)
  v2s = c(0, sort(U2s), 1)
  vs = c(0, sort(Us), 1) # i = 0 to 2n+1

  Omega_hatp = sigma_hatp * (nu - 2) / nu
  mu2_1p = rep(mu_hatp[2], n) + Omega_hatp[2,1] / Omega_hatp[1,1] * (Yp[,1] - rep(mu_hatp[1],n))
  a_p = (nu + (Yp[,1] - rep(mu_hatp[1],n))^2 / Omega_hatp[1,1] ) / (nu + 1)
  Omega2_1p = a_p * (Omega_hatp[2,2] - Omega_hatp[2,1] / Omega_hatp[1,1] * Omega_hatp[1,2])
  X1p = (Yp[,1] - mu_hatp[1]) / sqrt(Omega_hatp[1,1])
  X2p = (Yp[,2] - mu2_1p) / sqrt(Omega2_1p)
  U1p = pt(X1p, df = nu)
  U2p = pt(X2p, df = nu+1)
  Up = c(U1p, U2p)
  v1p = c(0, sort(U1p), 1)
  v2p = c(0, sort(U2p), 1)
  vp = c(0, sort(Up), 1) # i = 0 to 2n+1
  ################ Functions #################
  dg = function(r){
    return(c(1,
             -(nu + 1) * qt(r, nu) / (nu + (qt(r, nu))^2),
             1 - (nu + 1) * (qt(r, nu))^2 / (nu + (qt(r, nu))^2),
             - (nu + 2) * qt(r, nu + 1) / (nu + 1 + (qt(r, nu + 1))^2),
             1 - (nu + 2) * (qt(r, nu + 1))^2 / (nu + 1 + (qt(r, nu + 1))^2)))
  }
  ####### V method #####
  Ck = function(k, v){
    vi0 = v[(k+1):(2*n+1)]
    vi1 = v[(k+2):(2*n+2)]
    vi1_0 = vi1 - vi0
    # each column corresponds to a vectorized matrix of the result for v
    dgo = apply(array(vi0), 1, function(v) dg(v) %o% dg(v))
    vi1_0rep = repmat(t(vi1_0), 5^2, 1)
    vec = apply(dgo * vi1_0rep, 1, sum)
    mat = matrix(vec, 5, 5)
    return(mat)
  }
  Dk = function(k, v){
    vi0 = v[(k+1):(2*n+1)]
    dg1 = apply(array(vi0), 1, dg)
    result = apply(dg1, 1, sum)
    return(result)
  }
  W = function(j, v){
    A = j / (2 * n)
    sum = 0
    for (k in (1:j)){
      sum = sum + t(dg(v[k+1])) %*% solve(Ck(k, v)) %*% Dk(k, v) * (v[k+1] - v[k])
    }
    B = as.numeric(sum)
    return(sqrt(n*2) * abs(A - B / (n*2)))
  }
  
  # tic = Sys.time()
  # for (i in 1:(2*n)) print(c(i, W(i)))
  # toc = Sys.time()
  # cat("W time: ", toc - tic)
  # 
  # tic. = Sys.time()
  # for (i in 1:(2*n)) print(c(i, W.(i)))
  # toc. = Sys.time()
  # cat("W. time: ", toc. - tic.)
  
  # W.s = function(j){
  #   return(-W(j,vs))
  # }
  # W.p = function(j){
  #   return(-W(j,vp))
  # }
  # ctrl = DEoptim.control(itermax = 3, trace = F)
  # Sn.s = -DEoptim(W.s, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  # Sn.p = -DEoptim(W.p, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  
  maxlist_s = rep(0, 2*n-5) # i close to 2n, matrix inverse may not exist
  for (i in 1:(2*n-5)){
    maxlist_s[i] = W(i, vs)
  }
  Sn.s = max(maxlist_s)
  
  maxlist_p = rep(0, 2*n-5)
  for (i in 1:(2*n-5)){
    maxlist_p[i] = W(i, vp)
  }
  Sn.p = max(maxlist_p)
  
  
  ####### J method #######
  Ck.J = function(k, v){
    vi0 = v[(k+1):(n+1)]
    vi1 = v[(k+2):(n+2)]
    vi1_0 = vi1 - vi0
    # each column corresponds to a vectorized matrix of the result for v
    dgo = apply(array(vi0), 1, function(v) dg(v) %o% dg(v))
    vi1_0rep = repmat(t(vi1_0), 5^2, 1)
    vec = apply(dgo * vi1_0rep, 1, sum)
    mat = matrix(vec, 5, 5)
    return(mat)
  }
  Dk.J = function(k, v){
    vi0 = v[(k+1):(n+1)]
    dg1 = apply(array(vi0), 1, dg)
    result = apply(dg1, 1, sum)
    return(result)
  }
  W.J = function(j, v){
    A = j / n
    sum = 0
    for (k in (1:j)){
      sum = sum + t(dg(v[k+1])) %*% solve(Ck.J(k, v)) %*% Dk.J(k, v) * (v[k+1] - v[k])
    }
    B = as.numeric(sum)
    return(sqrt(n) * abs(A - B / n))
  }
  
  # W.J1s = function(j){
  #   return(-W.J(j,v1s))
  # }
  # W.J2s = function(j){
  #   return(-W.J(j,v2s))
  # }
  # W.J1p = function(j){
  #   return(-W.J(j,v1p))
  # }
  # W.J2p = function(j){
  #   return(-W.J(j,v2p))
  # }
  # ctrl = DEoptim.control(itermax = 3, trace = F)
  # W1.s = -DEoptim(W.J1s, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  # W2.s = -DEoptim(W.J2s, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  # W1.p = -DEoptim(W.J1p, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  # W2.p = -DEoptim(W.J2p, lower = 1, upper = 2*n-3, control = ctrl)$optim$bestval
  # Tn.s = max(W1.s, W2.s)
  # Tn.p = max(W1.p, W2.p)
  maxlist_s.J1 = rep(0, n-5)
  maxlist_s.J2 = rep(0, n-5)
  for (i in 1:(n-5)){
    maxlist_s.J1[i] = W.J(i, v1s)
    maxlist_s.J2[i] = W.J(i, v2s)
  }
  Tn.s = max(max(maxlist_s.J1), max(maxlist_s.J2))
  
  maxlist_p.J1 = rep(0, n-5)
  maxlist_p.J2 = rep(0, n-5)
  for (i in 1:(n-5)){
    maxlist_p.J1[i] = W.J(i, v1p)
    maxlist_p.J2[i] = W.J(i, v2p)
  }
  Tn.p = max(max(maxlist_p.J1), max(maxlist_p.J2))
  
  ####### K method
  # Kn.s = W1.s + W2.s
  # Kn.p = W1.p + W2.p
  Kn.s = max(maxlist_s.J1) + max(maxlist_s.J2)
  Kn.p = max(maxlist_p.J1) + max(maxlist_p.J2)    
  ################### return ###################
  return(c(
           # pval_s_skew, pval_p_skew,
           # pval_s_kurt, pval_p_kurt,
           # pval_s_sk, pval_p_sk,
           # pval_s_ku, pval_p_ku,
           # pval_s_skku, pval_p_skku,
           # pval_s_energy, pval_p_energy,
           # pval_s_micompr, pval_p_micompr,
           # pval_s_cramer, pval_p_cramer,
           # pval_s_JB, pval_p_JB,
           # pval_s_KSD_b, pval_p_KSD_b,
           pval_s_KSD, pval_p_KSD,
           Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           Tn.s > 2.993, Tn.s > 2.469, Tn.s > 2.211,
           Kn.s > 4.504, Kn.s > 3.792, Kn.s > 3.443,
           Sn.s, Sn.p, Tn.s, Tn.p, Kn.s, Kn.p))
           # Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
  # pb$tick()
  # Sys.sleep(1 / 100)
}

# ---------- Use foreach and doSNOW package ---------- #
n = 500
# n = 50
nrep = 100
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
                               # "energy",
                               'mvnTest',
                               'fMultivar',
                               'DEoptim',
                               'statmod',
                               'sn')) %dopar% {
                                 tryCatch(func(cnt),error = function(e) return(NA))
                               }


close(pb)
stopCluster(cl)
# result
res = result
nNA = sum(is.na(result[,1]))
result = na.omit(result)

cbind(apply(result[,1:2] < 0.01, 2, sum) / (nrep - nNA),
      apply(result[,1:2] < 0.05, 2, sum) / (nrep - nNA),
      apply(result[,1:2] < 0.10, 2, sum) / (nrep - nNA))



# apply(result[,3:11], 2, sum) / nrep
CV1_V = sort(result[,12])[(nrep - nNA)*0.9]
CV2_V = sort(result[,12])[(nrep - nNA)*0.95]
CV3_V = sort(result[,12])[(nrep - nNA)*0.99]
CV1_J = sort(result[,14])[(nrep - nNA)*0.9]
CV2_J = sort(result[,14])[(nrep - nNA)*0.95]
CV3_J = sort(result[,14])[(nrep - nNA)*0.99]
CV1_K = sort(result[,16])[(nrep - nNA)*0.9]
CV2_K = sort(result[,16])[(nrep - nNA)*0.95]
CV3_K = sort(result[,16])[(nrep - nNA)*0.99]

apply(result[,3:11], 2, sum) / (nrep - nNA)
c(sum(result[,13] > CV3_V) / (nrep - nNA),
  sum(result[,13] > CV2_V) / (nrep - nNA),
  sum(result[,13] > CV1_V) / (nrep - nNA))
c(sum(result[,15] > CV3_J) / (nrep - nNA),
  sum(result[,15] > CV2_J) / (nrep - nNA),
  sum(result[,15] > CV1_J) / (nrep - nNA))
c(sum(result[,17] > CV3_K) / (nrep - nNA),
  sum(result[,17] > CV2_K) / (nrep - nNA),
  sum(result[,17] > CV1_K) / (nrep - nNA))
