# Comparison of Various Normality Tests for GARCH
####################################################
######### remember to detach vars package ##########
####################################################
library(mvnTest) # Test for normality
library(MVN) # Test for normality
library(energy) # Test for normality
library(mvtnorm)
library(Matrix)
library(MTS)
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
d = 5
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

# GARCH = function(dist){
#   # Initialization
#   if (dist == 'N'){
#     ita = t(rmvnorm(n + ncut, mean = c(0, 0), sigma = diag(2), method = "chol"))
#   } else if (dist == 't'){
#     ita = t(rmvt(n = n + ncut, delta = c(0,0), sigma = (nu-2)/nu * diag(2),  df = nu))
#   }
#   y = array(rep(0,n*2),dim = c(2,n))
#   eps_t = c(0.5,0.5)
#   sigma_ii = c(0.3, 0.3)
#   # Iteration
#   for (t in 1:(n+ncut)){
#     sigma_ii = a + B %*% eps_t^2 + G %*% sigma_ii
#     sigma_21 = 0.5 * sqrt(sigma_ii[1]) * sqrt(sigma_ii[2])
#     sigma_t = matrix(c(sigma_ii[1], sigma_21,
#                        sigma_21, sigma_ii[2]),2,2)
#     eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
#     if (t > ncut){y[,t-ncut] = eps_t}
#   }
#   return(y)
# }
GARCH_d = function(dist){
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
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    if (t > ncut){y[,t-ncut] = eps_t}
  }
  return(y)
}
# GARCH. = function(dist, a., B., G., R.){
#   # Initialization
#   if (dist == 'N'){
#     ita = t(rmvnorm(n + ncut, mean = c(0, 0), sigma = diag(2), method = "chol"))
#   } else if (dist == 't'){
#     ita = t(rmvt(n = n + ncut, delta = c(0,0), sigma = (nu-2)/nu * diag(2),  df = nu))
#   }
#   y = array(rep(0,n*2),dim = c(2,n))
#   eps_t = c(0.5,0.5)
#   sigma_ii = c(0.3, 0.3)
#   r = R.[1,2]
#   # Iteration
#   for (t in 1:(n+ncut)){
#     sigma_ii = a. + B. %*% eps_t^2 + G. %*% sigma_ii
#     sigma_21 = r * sqrt(sigma_ii[1]) * sqrt(sigma_ii[2])
#     sigma_t = matrix(c(sigma_ii[1], sigma_21,
#                        sigma_21, sigma_ii[2]),2,2)
#     eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
#     if (t > ncut){y[,t-ncut] = eps_t}
#   }
#   return(y)
# }
GARCH_d. = function(dist, a., B., G., R.){
  # Initialization
  if (dist == 'N'){
    ita = t(rmvnorm(n + ncut, mean = rep(0, d), sigma = diag(d), method = "chol"))
  } else if (dist == 't'){
    ita = t(rmvt(n = n + ncut, delta = rep(0, d), sigma = (nu-2)/nu * diag(d),  df = nu))
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
    eps_t = msqrt(sigma_t)$mtxsqrt %*% ita[,t]
    if (t > ncut){y[,t-ncut] = eps_t}
  }
  return(y)
}

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
  } else if (N == 5.1){
    # Scenario 1 for d = 5
    a <<- rep(0.1, 5)
    B <<- matrix(c(0.3, 0.1, 0.1, 0.1, 0.1,
                   0.1, 0.2, 0.1, 0.1, 0.1,
                   0.1, 0.1, 0.25, 0.1, 0.1,
                   0.1, 0.1, 0.1, 0.15, 0.1,
                   0.1, 0.1, 0.1, 0.1, 0.1), 5, 5)
    G <<- matrix(c(0.2, 0.1, 0.01, 0.1, 0.01,
                   0.01, 0.3, 0.1, 0.1, 0.01,
                   0.01, 0.01, 0.1, 0.1, 0.01,
                   0.1, 0.01, 0.01, 0.15, 0.1,
                   0.01, 0.01, 0.1, 0.01, 0.2), 5, 5)
    R <<- matrix(0.5, d, d) + diag(1-0.5, d)
  }
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
nu = 5
Spec_GARCH(5.1)

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


n = 500
nB = 1000
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
# CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 17.22322)
# DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# d = 5
CP_ST = list(mean = rep(0, d), var.cov = diag(d), gamma1 = gamma_hat, gamma2 = 70.55935)
DP_ST = cp2dp(CP_ST, "ST")
# DP_ST$nu
# skew-t score function
# sqx = function(x){
#   h = 0.00001
#   df = rep(0, d)
#   for (i in 1:d){
#     vec01 = rep(0,d)
#     vec01[i] = 1
#     # five-point method
#     df[i] = (-dmst(x = x + 2 * vec01 * h, dp = DP_ST, log = T) +
#                8 * dmst(x = x + vec01 * h, dp = DP_ST, log = T) -
#                8 * dmst(x = x - vec01 * h, dp = DP_ST, log = T) +
#                dmst(x = x - 2 * vec01 * h, dp = DP_ST, log = T)) / 12 / h
#   }
#   
#   return(df)
# }


# pb <- progress_bar$new(
#   format = "  Percent [:bar] :percent Time :elapsed",
#   total = nrep, clear = FALSE, width= 60)
lambda = 1
# plist_size_KSD = NULL
# plist_power_KSD = NULL
# for (count in (1:nrep)){
func = function(cnt){
  ################# Generate Data #################
  # set.seed(count + 666)
  # Generate Sample Data
  # ccc_s = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
  #                  d.f = Inf, cut = ncut, model = "diagonal")
  # Ys = ccc_s$eps
  Ys = t(GARCH_d('N'))
  # ccc_p = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
  #                  d.f = nu, cut = ncut, model = "diagonal")
  # Yp = ccc_p$eps
  Yp = t(GARCH_d('N+logN'))
  ################# Estimation #################
  # est_s = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
  #                             R = R_init, dvar = Ys, model = "extended"),
  #             silent = TRUE)
  # if ('try-error' %in% class(est_s)) {return(c(NA, NA))}
  est_s = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                   R = R_init, dvar = Ys, model = "extended"),
                   error = function(e) {return(NA)},
                   warning = function(w) NA)
  if (is.na(est_s[1])){return(NA)}
  # if (is.na(est_s[1])){next}
  # est_p = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
  #                             R = R_init, dvar = Yp, model = "extended"),
  #             silent = TRUE)
  # if ('try-error' %in% class(est_p)) {return(c(NA, NA))}
  est_p = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                   R = R_init, dvar = Yp, model = "extended"),
                   error = function(e) {return(NA)},
                   warning = function(w) {return(NA)})
  if (is.na(est_p[1])) {return(NA)}
  # if (is.na(est_p[1])) {next}
  
  Ys.res = est_s$std.resid
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
  Rs = est_s$para.mat$R
  Rp = est_p$para.mat$R

  Rs.invsqrt = tryCatch(msqrt(Rs)$invsqrt,
                                warning = function(w) matrix(NA, d, d),
                                error = function(e) matrix(NA, d, d))
  if (is.na(Rs.invsqrt[1,1])) return(NA)
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                                warning = function(w) matrix(NA, d, d),
                                error = function(e) matrix(NA, d, d))
  if (is.na(Rp.invsqrt[1,1])) return(NA)

  Xs = Rs.invsqrt %*% t(Ys.res)
  Xp = Rp.invsqrt %*% t(Yp.res)

  ################### Skew / Kurt Test ###################
  # b1_hats = apply(Xs^3, 1, mean)
  # b1_hatp = apply(Xp^3, 1, mean)
  # b2_hats = apply(Xs^4, 1, mean)
  # b2_hatp = apply(Xp^4, 1, mean)
  # 
  # lambda_s.sk = n * sum(b1_hats^2) / 6 
  # lambda_p.sk = n * sum(b1_hatp^2) / 6
  # lambda_s.ku = n * sum((b2_hats - 3)^2) / 24
  # lambda_p.ku = n * sum((b2_hatp - 3)^2) / 24
  # 
  # pval_s_skew = 1 - pchisq(lambda_s.sk, d)
  # pval_s_kurt = 1 - pchisq(lambda_s.ku, d)
  # pval_p_skew = 1 - pchisq(lambda_p.sk, d)
  # pval_p_kurt = 1 - pchisq(lambda_p.ku, d)
  # pval_s_skku = 1 - pchisq(lambda_s.sk + lambda_s.ku, 2 * d)
  # pval_p_skku = 1 - pchisq(lambda_p.sk + lambda_p.ku, 2 * d)
  # 
  ################### DH Test ###################
  # result_s_DH = DH.test(data = t(Xs), qqplot = FALSE)
  # result_p_DH = DH.test(data = t(Xp), qqplot = FALSE)
  # 
  # # plist_size_DH = c(plist_size_DH, result_s_DH@p.value)
  # # plist_power_DH = c(plist_power_DH, result_p_DH@p.value)
  # 
  ################### AD Test ###################
  # result_s_AD = AD.test(data = t(Xs), qqplot = FALSE)
  # result_p_AD = AD.test(data = t(Xp), qqplot = FALSE)
  
  # plist_size_AD = c(plist_size_AD, result_s_AD@p.value)
  # plist_power_AD = c(plist_power_AD, result_p_AD@p.value)
  
  ################### CM Test ###################
  # result_s_CM = CM.test(data = t(Xs), qqplot = FALSE)
  # result_p_CM = CM.test(data = t(Xp), qqplot = FALSE)
  
  # plist_size_CM = c(plist_size_CM, result_s_CM@p.value)
  # plist_power_CM = c(plist_power_CM, result_p_CM@p.value)
  
  ################### HZ Test ###################
  # result_s_HZ = HZ.test(data = t(Xs), qqplot = FALSE)
  # result_p_HZ = HZ.test(data = t(Xp), qqplot = FALSE)
  # 
  # # plist_size_HZ = c(plist_size_HZ, result_s_HZ@p.value)
  # # plist_power_HZ = c(plist_power_HZ, result_p_HZ@p.value)
  # 
  ################### R Test ###################
  # result_s_R = R.test(data = t(Xs), qqplot = FALSE)
  # result_p_R = R.test(data = t(Xp), qqplot = FALSE)
  # 
  # # plist_size_R = c(plist_size_R, result_s_R@p.value)
  # # plist_power_R = c(plist_power_R, result_p_R@p.value)
  # 
  ################### Mardia Test ###################
  # result_s_mardia = mvn(data = t(Xs), mvnTest = "mardia")
  # result_p_mardia = mvn(data = t(Xp), mvnTest = "mardia")
  # 
  # # plist_size_mardia_sk = c(plist_size_mardia_sk, as.numeric(levels(result_s_mardia[["multivariateNormality"]][["p value"]]))[1])
  # # plist_power_mardia_sk = c(plist_power_mardia_sk, as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[1])
  # # plist_size_mardia_ku = c(plist_size_mardia_ku, as.numeric(levels(result_s_mardia[["multivariateNormality"]][["p value"]]))[2])
  # # plist_power_mardia_ku = c(plist_power_mardia_ku, as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[2])
  # 
  ################### Henze Test ###################
  # result_s_henze = mvn(data = t(Xs), mvnTest = "hz")
  # result_p_henze = mvn(data = t(Xp), mvnTest = "hz")
  # 
  # # plist_size_henze = c(plist_size_henze, result_s_henze[["multivariateNormality"]][["p value"]])
  # # plist_power_henze = c(plist_power_henze, result_p_henze[["multivariateNormality"]][["p value"]])
  # 
  ################### Roy Test ###################
  # result_s_roy = mvn(data = t(Xs), mvnTest = "royston")
  # result_p_roy = mvn(data = t(Xp), mvnTest = "royston")
  # 
  # # plist_size_roy = c(plist_size_roy, result_s_roy[["multivariateNormality"]][["p value"]])
  # # plist_power_roy = c(plist_power_roy, result_p_roy[["multivariateNormality"]][["p value"]])
  # 
  ################### doh Test ###################
  # result_s_doh = mvn(data = t(Xs), mvnTest = "dh")
  # result_p_doh = mvn(data = t(Xp), mvnTest = "dh")
  # 
  # # plist_size_doh = c(plist_size_doh, result_s_doh[["multivariateNormality"]][["p value"]])
  # # plist_power_doh = c(plist_power_doh, result_p_doh[["multivariateNormality"]][["p value"]])
  # 
  ################### Energy Test ###################
  # result_s_energy = mvn(data = t(Xs), mvnTest = "energy")
  # result_p_energy = mvn(data = t(Xp), mvnTest = "energy")
  # 
  # # plist_size_energy = c(plist_size_energy, result_s_energy[["multivariateNormality"]][["p value"]])
  # # plist_power_energy= c(plist_power_energy, result_p_energy[["multivariateNormality"]][["p value"]])
  # 
  # # result_s_energy = mvnorm.etest(t(Xs),1000)
  # # result_p_energy = mvnorm.etest(t(Xp),1000)
  # ################## JB Test ##################
  # # var(t(Xs))
  # g.Xs = Xs^3 - 3 * Xs
  # JB1_s = sum((apply(g.Xs, 1, sum))^2) / (6 * n)
  # h.Xs = Xs^4 - 6 * Xs^2 + 3
  # JB2_s = sum(h.Xs) / sqrt(24 * n * d)
  # pval_s_JB = 1 - pchisq(JB1_s + JB2_s^2, d + 1)
  # 
  # g.Xp = Xp^3 - 3 * Xp
  # JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * n)
  # h.Xp = Xp^4 - 6 * Xp^2 + 3
  # JB2_p = sum(h.Xp) / sqrt(24 * n * d)
  # pval_p_JB = 1 - pchisq(JB1_p + JB2_p^2, d + 1)
  # 
  ################### KSD Test ###################
  as = est_s$para.mat$a
  Bs = est_s$para.mat$A
  Gs = est_s$para.mat$B
  ap = est_p$para.mat$a
  Bp = est_p$para.mat$A
  Gp = est_p$para.mat$B
  # Rs = est_s$para.mat$R
  # Rp = est_p$para.mat$R
  # Observed Test Statistics
  Ts_obs = KSD_new(x = t(Xs), score_function = "gaussian", width = -1)
  Ts_obs = as.numeric(Ts_obs)
  Tp_obs = KSD_new(x = t(Xp), score_function = "gaussian", width = -1)
  Tp_obs = as.numeric(Tp_obs)

  # Bootstrap
  bootssample_s = NULL
  bootssample_p = NULL
  for (b in 1:nB){
    # Resample
    # set.seed(b+1)
    # Recover Ys_b (n x d)
    Ys_b = t(GARCH_d.('N', as, Bs, Gs, Rs))
    Yp_b = t(GARCH_d.('N', ap, Bp, Gp, Rp))
    # Estimation
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

    ################# M1 ################
    # #Estimated Mean
    # mu_hats_b = colMeans(Ys_b.res)
    # mu_hatp_b = colMeans(Yp_b.res)
    # #Estimated Sigma
    # sigma_hats_b = var(Ys_b.res)
    # sigma_hatp_b = var(Yp_b.res)
    #
    # sigma_hats_b.invsqrt = tryCatch(msqrt(sigma_hats_b)$invsqrt,
    #                               warning = function(w) matrix(NA, d, d),
    #                               error = function(e) matrix(NA, d, d))
    # # if (is.na(Rs.invsqrt)) {next}
    # if (is.na(sigma_hats_b.invsqrt[1,1])) next
    # sigma_hatp_b.invsqrt = tryCatch(msqrt(sigma_hatp_b)$invsqrt,
    #                               warning = function(w) matrix(NA, d, d),
    #                               error = function(e) matrix(NA, d, d))
    # # if (is.na(Rp.invsqrt)) {next}
    # if (is.na(sigma_hatp_b.invsqrt[1,1])) next
    #
    # Xs_b = sigma_hats_b.invsqrt %*% (t(Ys_b.res) - mu_hats_b)
    # Xp_b = sigma_hatp_b.invsqrt %*% (t(Yp_b.res) - mu_hatp_b)
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

    # KSD Test Statistic for bootstrap data
    Ts_b = KSD_new(x = t(Xs_b), score_function = "gaussian", width = -1)
    Ts_b = as.numeric(Ts_b)
    Tp_b = KSD_new(x = t(Xp_b), score_function = "gaussian", width = -1)
    Tp_b = as.numeric(Tp_b)
    bootssample_s = c(bootssample_s, Ts_b)
    bootssample_p = c(bootssample_p, Tp_b)
  }

  pval_s_KSD = sum((bootssample_s - Ts_obs) > 0) / length(bootssample_s)
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
  ################### KSD Block Test ###################
  # # Block
  # X1s = Xs[,1:m1]
  # X2s = Xs[,(n-m2+1):n]
  # 
  # X1p = Xp[,1:m1]
  # X2p = Xp[,(n-m2+1):n]
  # 
  # # KSD Package Results
  # result1s <- KSD(t(X1s),score_function = "gaussian", 'rbf',-1.0)
  # result2s <- KSD(t(X2s),score_function = "gaussian", 'rbf',-1.0)
  # Ts1obs = result1s[["ksd"]]
  # Ts2obs = result2s[["ksd"]]
  # Tsobs = w * Ts1obs + (1-w) * Ts2obs
  # bootssample_s = w * result1s[["bootStrapSamples"]] + (1-w) * result2s[["bootStrapSamples"]]
  # pval_s_KSD_b = sum((bootssample_s - Tsobs) > 0) / 1000 # nboots=1000
  # result1p <- KSD(t(X1p),score_function = "gaussian",'rbf',-1.0)
  # result2p <- KSD(t(X2p),score_function = "gaussian",'rbf',-1.0)
  # Tp1obs = result1p[["ksd"]]
  # Tp2obs = result2p[["ksd"]]
  # Tpobs = w * Tp1obs + (1-w) * Tp2obs
  # bootssample_p = w * result1p[["bootStrapSamples"]] + (1-w) * result2p[["bootStrapSamples"]]
  # pval_p_KSD_b = sum((bootssample_p - Tpobs) > 0) / 1000
  
  
  results = KSD(t(Xs), score_function = "gaussian", 'rbf',-1.0)
  pval_s_KSD = results$p
  resultp = KSD(t(Xp), score_function = "gaussian", 'rbf',-1.0)
  pval_p_KSD = resultp$p
  # ################### Bai Test ###################
  # ### for normal GARCH
  # # Normaldata = rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol")
  # 
  # sigma_ii.s = est_s[["h"]]
  # sigma_12.s = sqrt(sigma_ii.s[,1]) * sqrt(sigma_ii.s[,2]) * Rs[1,2]
  # sigma2_1.s = sqrt(sigma_ii.s[,2] - sigma_12.s^2 / sigma_ii.s[,1])
  # X1s = Ys[,1] / sqrt(sigma_ii.s[,1])
  # X2s = (Ys[,2] - sigma_12.s / sigma_ii.s[,1] * Ys[,1]) / sigma2_1.s
  # U1s = pnorm(X1s, 0, 1)
  # U2s = pnorm(X2s, 0, 1)
  # U1s[which(U1s == 1)] = 0.9999999
  # U2s[which(U2s == 1)] = 0.9999999
  # 
  # sigma_ii.p = est_p[["h"]]
  # sigma_12.p = sqrt(sigma_ii.p[,1]) * sqrt(sigma_ii.p[,2]) * Rp[1,2]
  # sigma2_1.p = sqrt(sigma_ii.p[,2] - sigma_12.p^2 / sigma_ii.p[,1])
  # X1p = Yp[,1] / sqrt(sigma_ii.p[,1])
  # X2p = (Yp[,2] - sigma_12.p / sigma_ii.p[,1] * Yp[,1]) / sigma2_1.p
  # U1p = pnorm(X1p, 0, 1)
  # U2p = pnorm(X2p, 0, 1)
  # 
  # U1p[which(U1p == 1)] = 0.9999999
  # U2p[which(U2p == 1)] = 0.9999999
  # 
  # # ds1 = density(X1s)
  # # ds2 = density(X2s)
  # # dp1 = density(X1p)
  # # dp2 = density(X2p)
  # # d0 = density(Normaldata)
  # # 
  # # # Comparison
  # # par(mfrow = c(1,1))
  # # plot(ds1, type="n", main= "Compare", ylim = c(0, 1), xlim = c(-5,5))
  # # polygon(ds1)
  # # polygon(dp1, border = 'red')
  # # polygon(ds2, border = 'green')
  # # polygon(dp2, border = 'blue')
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
  # 
  # 
  ###################################################
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
           pval_s_KSD, pval_p_KSD))
           # Sn.s > 2.787, Sn.s > 2.214, Sn.s > 1.940,
           # Sn.p > 2.787, Sn.p > 2.214, Sn.p > 1.940))
           # pval_s_KSD_b, pval_p_KSD_b))
  # pb$tick()
  # Sys.sleep(1 / 100)
}

#################### Size Evaluation ###################
# ---------- Use foreach and doSNOW package ---------- #
n = 500
lambda = 1

cl <- makeCluster(50)
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
                               "ccgarch",
                               "mvnTest",
                               "DEoptim",
                               "fMultivar",
                               "statmod",
                               "MVN",
                               "energy",
                               "sn")) %dopar% {
                                 tryCatch(func(cnt),error = function(e) return(NA))
                               }


close(pb)
stopCluster(cl)

# KSD Block Test Result
n.NA = sum(is.na(result[,1])) # NA number
index = which(is.na(result[,1])) # NA row index
if (n.NA > 0){
  temp = result[- index,] # remove NA rows
} else {
  temp = result
}
cat("n =",n, "lambda =",lambda)
cbind(apply(temp[,1:2]<0.01, 2, sum) / (nrep - n.NA),
      apply(temp[,1:2]<0.05, 2, sum) / (nrep - n.NA),
      apply(temp[,1:2]<0.10, 2, sum) / (nrep - n.NA))

# result
# cat("n =",n, "lambda =",lambda)
# cbind(apply(result[,1:28]<0.01, 2, sum) / nrep,
#       apply(result[,1:28]<0.05, 2, sum) / nrep,
#       apply(result[,1:28]<0.10, 2, sum) / nrep)
# apply(result[,29:34], 2, sum) / nrep


n.NA = sum(is.na(result[,1])) # NA number
index = which(is.na(result[,1])) # NA row index
if (n.NA > 0){
  temp = result[- index,] # remove NA rows
} else {
  temp = result
}

n.NA8 = sum(is.na(temp[,16]))
if (n.NA8 > 0){
  ind = which(is.na(temp[,16]))
  temp = temp[-ind,]
}
cat("n =",n, "lambda =",lambda)
cbind(apply(temp[,1:28]<0.01, 2, sum) / (nrep - n.NA - n.NA8),
      apply(temp[,1:28]<0.05, 2, sum) / (nrep - n.NA - n.NA8),
      apply(temp[,1:28]<0.10, 2, sum) / (nrep - n.NA - n.NA8))
apply(temp[,29:34], 2, sum) / (nrep - n.NA - n.NA8)





###################################################
###################### POWER ######################
###################################################
func_p = function(cnt){
  ################# Generate Data #################
  # set.seed(count + 666)
  # Generate Sample Data
  # ccc_p = eccc.sim(nobs = n, a = a, A = B, B = G, R = R,
  #                  d.f = nu, cut = ncut, model = "diagonal")
  # Yp = ccc_p$eps
  Yp = t(GARCH_d('N+logN'))
  ################# Estimation #################
  # est_p = try(eccc.estimation(a = a_init, A = B_init, B = G_init,
  #                             R = R_init, dvar = Yp, model = "extended"),
  #             silent = TRUE)
  # if ('try-error' %in% class(est_p)) {return(c(NA, NA))}
  est_p = tryCatch(eccc.estimation(a = a_init, A = B_init, B = G_init,
                                   R = R_init, dvar = Yp, model = "extended"),
                   error = function(e) {return(NA)},
                   warning = function(w) {return(NA)})
  if (is.na(est_p[1])) {return(NA)}
  # if (is.na(est_p[1])) {next}
  
  Yp.res = est_p$std.resid
  
  ########################## M1 ##########################
  # #Estimated Mean
  # mu_hatp = colMeans(Yp.res)
  # #Estimated Sigma
  # sigma_hatp = var(Yp.res)
  # 
  # sigma_hatp.invsqrt = tryCatch(msqrt(sigma_hatp)$invsqrt,
  #                       warning = function(w) matrix(NA, d, d),
  #                       error = function(e) matrix(NA, d, d))
  # if (is.na(sigma_hatp.invsqrt[1,1])) return(c(NA, NA))
  # 
  # Xp = sigma_hatp.invsqrt %*% (t(Yp.res) - mu_hatp)
  
  ########################## M2 #########################
  Rp = est_p$para.mat$R
  
  Rp.invsqrt = tryCatch(msqrt(Rp)$invsqrt,
                        warning = function(w) matrix(NA, d, d),
                        error = function(e) matrix(NA, d, d))
  if (is.na(Rp.invsqrt[1,1])) return(NA)
  
  Xp = Rp.invsqrt %*% t(Yp.res)
  
  ################### Skew / Kurt Test ###################
  b1_hatp = apply(Xp^3, 1, mean)
  b2_hatp = apply(Xp^4, 1, mean)
  
  lambda_p.sk = n * sum(b1_hatp^2) / 6
  lambda_p.ku = n * sum((b2_hatp - 3)^2) / 24
  
  pval_p_skew = 1 - pchisq(lambda_p.sk, d)
  pval_p_kurt = 1 - pchisq(lambda_p.ku, d)
  pval_p_skku = 1 - pchisq(lambda_p.sk + lambda_p.ku, 2 * d)
  
  ################### DH Test ###################
  result_p_DH = DH.test(data = t(Xp), qqplot = FALSE)
  ################### AD Test ###################
  # result_p_AD = AD.test(data = t(Xp), qqplot = FALSE)
  ################### CM Test ###################
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
  ################## JB Test ##################
  # var(t(Xs))
  g.Xp = Xp^3 - 3 * Xp
  JB1_p  = sum((apply(g.Xp, 1, sum))^2) / (6 * n)
  h.Xp = Xp^4 - 6 * Xp^2 + 3
  JB2_p = sum(h.Xp) / sqrt(24 * n * d)
  pval_p_JB = 1 - pchisq(JB1_p + JB2_p^2, d + 1)
  
  ################### KSD Test ###################
  ap = est_p$para.mat$a
  Bp = est_p$para.mat$A
  Gp = est_p$para.mat$B
  # Rp = est_p$para.mat$R
  # Observed Test Statistics
  Tp_obs = KSD_new(x = t(Xp), score_function = "gaussian", width = -1)
  Tp_obs = as.numeric(Tp_obs)
  
  # Bootstrap
  bootssample_p = NULL
  for (b in 1:nB){
    # Resample
    # set.seed(b+1)
    # Recover Ys_b (n x d)
    Yp_b = t(GARCH_d.('N', ap, Bp, Gp, Rp))
    # Estimation
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
    
    ################# M1 ################
    # #Estimated Mean
    # mu_hatp_b = colMeans(Yp_b.res)
    # #Estimated Sigma
    # sigma_hatp_b = var(Yp_b.res)
    # 
    # sigma_hatp_b.invsqrt = tryCatch(msqrt(sigma_hatp_b)$invsqrt,
    #                               warning = function(w) matrix(NA, d, d),
    #                               error = function(e) matrix(NA, d, d))
    # # if (is.na(Rp.invsqrt)) {next}
    # if (is.na(sigma_hatp_b.invsqrt[1,1])) next
    # 
    # Xp_b = sigma_hatp_b.invsqrt %*% (t(Yp_b.res) - mu_hatp_b)
    ################# M2 ################
    Rp_b = est_p.b$para.mat$R
    
    Rp_b.invsqrt = tryCatch(msqrt(Rp_b)$invsqrt,
                            warning = function(w) matrix(NA,d,d),
                            error = function(e) matrix(NA,d,d))
    if (is.na(Rp_b.invsqrt[1,1])) {next}
    
    Xp_b = Rp_b.invsqrt %*% t(Yp_b.res)
    
    # KSD Test Statistic for bootstrap data
    Tp_b = KSD_new(x = t(Xp_b), score_function = "gaussian", width = -1)
    Tp_b = as.numeric(Tp_b)
    bootssample_p = c(bootssample_p, Tp_b)
  }
  
  pval_p_KSD = sum((bootssample_p - Tp_obs) > 0) / length(bootssample_p)
  # ################### Bai Test ###################
  ### for normal GARCH
  sigma_ii.p = est_p[["h"]]
  sigma_12.p = sqrt(sigma_ii.p[,1]) * sqrt(sigma_ii.p[,2]) * Rp[1,2]
  sigma2_1.p = sqrt(sigma_ii.p[,2] - sigma_12.p^2 / sigma_ii.p[,1])
  X1p = Yp[,1] / sqrt(sigma_ii.p[,1])
  X2p = (Yp[,2] - sigma_12.p / sigma_ii.p[,1] * Yp[,1]) / sigma2_1.p
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
  
  
  ###################################################
  return(c(pval_p_skew,
           pval_p_kurt,
           pval_p_skku,
           result_p_DH@p.value,
           # result_p_AD@p.value,
           # result_p_CM@p.value,
           result_p_HZ@p.value,
           result_p_R@p.value,
           as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[1],
           as.numeric(levels(result_p_mardia[["multivariateNormality"]][["p value"]]))[2],
           result_p_henze[["multivariateNormality"]][["p value"]],
           result_p_roy[["multivariateNormality"]][["p value"]],
           result_p_doh[["multivariateNormality"]][["p value"]],
           result_p_energy[["multivariateNormality"]][["p value"]],
           pval_p_JB,
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
                               "ccgarch")) %dopar% {
                                 tryCatch(func_p(cnt),error = function(e) return(NA))
                               }

close(pb)
stopCluster(cl)
# result
# cat("n =",n, "lambda =",lambda)
# cbind(apply(result[,1:14]<0.01, 2, sum) / nrep,
#       apply(result[,1:14]<0.05, 2, sum) / nrep,
#       apply(result[,1:14]<0.10, 2, sum) / nrep)
# apply(result[,15:17], 2, sum) / nrep

n.NA = sum(is.na(result[,1])) # NA number
index = which(is.na(result[,1])) # NA row index
if (n.NA > 0){
  temp = result[- index,] # remove NA rows
} else {
  temp = result
}

n.NA8 = sum(is.na(temp[,8]))
if (n.NA8 > 0){
  ind = which(is.na(temp[,8]))
  temp = temp[-ind,]
}
cat("n =",n, "lambda =",lambda)
cbind(apply(temp[,1:14]<0.01, 2, sum) / (nrep - n.NA - n.NA8),
      apply(temp[,1:14]<0.05, 2, sum) / (nrep - n.NA - n.NA8),
      apply(temp[,1:14]<0.10, 2, sum) / (nrep - n.NA - n.NA8))
apply(temp[,15:17], 2, sum) / (nrep - n.NA - n.NA8)

########################################################
#################### not in use ########################
########################################################
plist_size_skew = result[,1]
plist_power_skew = result[,2]
plist_size_kurt = result[,3]
plist_power_kurt = result[,4]
plist_size_skku = result[,5]
plist_power_skku = result[,6]
plist_size_DH = result[,7]
plist_power_DH = result[,8]
# plist_size_AD = result[,]
# plist_power_AD = result[,]
# plist_size_CM = result[,]
# plist_power_CM = result[,]
plist_size_HZ = result[,9]
plist_power_HZ = result[,10]
plist_size_R = result[,11]
plist_power_R = result[,12]
plist_size_mardia_sk = result[,13]
plist_power_mardia_sk = result[,14]
plist_size_mardia_ku = result[,15]
plist_power_mardia_ku = result[,16]
plist_size_henze = result[,17]
plist_power_henze = result[,18]
plist_size_roy = result[,19]
plist_power_roy = result[,20]
plist_size_doh = result[,21]
plist_power_doh = result[,22]
plist_size_energy = result[,23]
plist_power_energy = result[,24]
plist_size_JB = result[,25]
plist_power_JB = result[,26]
plist_size_JB_m = result[,27]
plist_power_JB_m = result[,28]
plist_size_KSD = result[,29]
plist_power_KSD = result[,30]
eval_test = function(plist_size, plist_power){
  #Remove NAs
  ns = sum(is.na(plist_size))
  np = sum(is.na(plist_power))
  plist_size[is.na(plist_size)] = 1
  plist_power[is.na(plist_power)] = 1
  size1 = sum(plist_size - 0.01 < 0) / (length(plist_size) - ns)
  size5 = sum(plist_size - 0.05 < 0) / (length(plist_size) - ns)
  size10 = sum(plist_size - 0.10 < 0) /(length(plist_size) - ns)
  power1 = sum(plist_power - 0.01 < 0) / (length(plist_power) - np)
  power5 = sum(plist_power - 0.05 < 0) / (length(plist_power) - np)
  power10 = sum(plist_power - 0.10 < 0) / (length(plist_power) - np)
  
  # size1 = sum(plist_size - 0.01 < 0) / nrep
  # size5 = sum(plist_size - 0.05 < 0) / nrep
  # size10 = sum(plist_size - 0.10 < 0) / nrep
  # power1 = sum(plist_power - 0.01 < 0) / nrep
  # power5 = sum(plist_power - 0.05 < 0) / nrep
  # power10 = sum(plist_power - 0.10 < 0) / nrep
  print(c(size1, size5, size10, power1, power5, power10))
}

eval_test(plist_size_skew, plist_power_skew)
eval_test(plist_size_kurt, plist_power_kurt)
eval_test(plist_size_skku, plist_power_skku)
eval_test(plist_size_DH, plist_power_DH)
# eval_test(plist_size_AD, plist_power_AD)
# eval_test(plist_size_CM, plist_power_CM)
eval_test(plist_size_HZ, plist_power_HZ)
eval_test(plist_size_R, plist_power_R)
eval_test(plist_size_mardia_sk, plist_power_mardia_sk)
eval_test(plist_size_mardia_ku, plist_power_mardia_ku)
eval_test(plist_size_henze, plist_power_henze)
eval_test(plist_size_roy, plist_power_roy)
eval_test(plist_size_doh, plist_power_doh)
eval_test(plist_size_energy, plist_power_energy)
eval_test(plist_size_JB, plist_power_JB)
eval_test(plist_size_JB_m, plist_power_JB_m)
eval_test(plist_size_KSD, plist_power_KSD)

