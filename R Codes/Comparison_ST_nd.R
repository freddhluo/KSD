# Comparison of tests
# for Constant Mean Constant Variance 
# dimension = d

library(mvnTest)
library(mvtnorm)
library(Matrix)
library(MTS)
library(progress)
# library(KSD)
library(pryr)
library(compositions)
library(MVN)
library(energy)
library(DEoptim)
library(fMultivar) # rcauchy2d
library(statmod) # gauss.quad() used

library(foreach)
library(doSNOW)

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
  T_obs = KSD_new.(x = t(X), score_function = t(apply(t(X), 1, sqx)), width = -1)
  T_obs = as.numeric(T_obs)
  
  # Distribution of T
  sample = NULL
  for (b in 1:nB.){
    # Generate residual under H0
    # set.seed(b+1)
    ita.N = rmst(n, dp = DP_ST)
    # Recover Ys_b (n x d)
    Y_b = t(mu_hat + msqrt(sigma_hat)$mtxsqrt %*% t(ita.N))
    # Estimation
    mu_hat_b = colMeans(Y_b)
    sigma_hat_b = var(Y_b)
    
    theta_hat_b = c(mu_hat_b, c(sigma_hat_b))
    X_b = g.(theta_hat_b, t(Y_b))
    # KSD Test Statistic for bootstrap data
    T_b = KSD_new.(x = t(X_b), score_function = t(apply(t(X_b), 1, sqx)), width = -1)
    T_b = as.numeric(T_b)
    sample = c(sample, T_b)
  }
  
  pval_KSD = sum((sample - T_obs) > 0) / nB.
  return(pval_KSD)
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



# Sample Size & Block length
Spec_Nsample = function(N, w = 0.3){
  if (N == 1){
    n <<- 100
    m1 <<- n * w
    m2 <<- n * (1 - w)
  } else if (N == 2){
    n <<- 200
    m1 <<- n * w
    m2 <<- n * (1 - w)
  } else if (N == 3){
    n <<- 500
    m1 <<- n * w
    m2 <<- n * (1 - w)
  }
}
#Mean
mu = rep(0,d)
# mu_log = rep(exp(1/2), 2)

#Variance-Covariance Matrix
sigma = matrix(c(1,0.5,0.5,1), ncol = 2)
# sigma = matrix(c(1, 0.5, 0.25, 0.125, 0.0625,
#                  0.5, 1, 0.5, 0.25, 0.125,
#                  0.25, 0.5, 1, 0.5, 0.25,
#                  0.125, 0.25, 0.5, 1, 0.5,
#                  0.0625, 0.125, 0.25, 0.5, 1), d, d)
# sigma_log = matrix(c(exp(2)-exp(1), exp(1.5)-exp(1),
#                      exp(1.5)-exp(1), exp(2)-exp(1)),
#                    ncol = 2)

# Weights
# w = 0.3

#Settings for H0 model
# mu0 = matrix(c(0,0),ncol =1)
# sigma0 = array(diag(2),c(2,2,1))
# model = gmm(nComp = 1, mu = mu0, sigma = sigma0, d = 2)
# score_function = partial(scorefunctiongmm, model=model)

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
nu = 5
cosi = c(1, 1.3)
# Spec_Nsample(1, 0.3)
# Spec_Nsample(2, 0.3)
# Spec_Nsample(3, 0.3)

# skew-normal parameters
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

################### Execution ###################
func = function(cnt){
  # for (count in (1:nrep)){
  ################# Generate Data #################
  # set.seed(count + 12345)
  # Generate Sample Data
  Ys = rmst(n, dp = DP_ST) # Skew t distibution
  # Ys = rmsn(n, dp = DP)
  # Ys = rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol")
  Ys = t(mu + msqrt(sigma)$mtxsqrt %*% t(Ys))
  
  
  # Skew t distibution
  # Yp = rmst(n, dp = DP_ST)
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # Skew normal distibution
  Yp = rmsn(n, dp = DP)
  Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # Yp = logN(n, mu_log, sigma_log)
  # Yp = rmvt(n, df = 5)
  # Yp = rmvt(n = n, delta = rep(0,d), sigma = sigma,  df = 5)
  
  # Normal as H1 distribution
  # Yp = rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol")
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # [(1 - lambda) x N + lambda x t(5)]
  # Yp = ((1 - lambda) * rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol") +
  #   lambda * rmvt(n = n, delta = rep(0,d), sigma = (nu - 2) / nu * diag(d), df = nu)) /
  #   sqrt(lambda^2 + (1 - lambda)^2)
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
  # [(1 - lambda) x N + lambda x logN]
  # temp = t(matrix(rlnorm.rplus(n, meanlog = log(rep(0.3, d)), varlog = diag(d)), nrow = n))
  # temp = t(msqrt(var(t(temp)))$invsqrt %*% (temp - rowMeans(temp)))
  # Yp = ((1 - lambda) * rmvnorm(n, mean = rep(0,d), sigma = diag(d), method = "chol") +
  #         lambda * temp) / 
  #   sqrt(lambda^2 + (1 - lambda)^2)
  # Yp = t(mu + msqrt(sigma)$mtxsqrt %*% t(Yp))
  
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
  
  # var(Yp)
  # ds = density(Ys[,1])
  # dp = density(Yp[,1])
  # par(mfrow = c(1,2))
  # plot(ds, type="n", main="Ys")
  # polygon(ds)
  # plot(dp, type = "n", main = "Yp")
  # polygon(dp)
  ################ Visualize Ys #################
  # library(MASS)
  # # Calculate kernel density estimate
  # Ys.kde <- kde2d(Ys[,1], Ys[,2], n = 50)   # from MASS package
  # 
  # # Contour plot overlayed on heat map image of results
  # # image(Ys.kde)       # from base graphics package
  # # contour(Ys.kde, add = TRUE)     # from base graphics package
  # 
  # # Classic Bivariate Normal Diagram
  # library(ellipse)
  # # rho <- cor(Ys)
  # # y_on_x <- lm(Ys[,2] ~ Ys[,1])    # Regressiion Y ~ X
  # # x_on_y <- lm(Ys[,1] ~ Ys[,2])    # Regression X ~ Y
  # # plot_legend <- c("99% CI green", "95% CI red","90% CI blue",
  # #                  "Y on X black", "X on Y brown")
  # # 
  # # plot(Ys, xlab = "X", ylab = "Y",
  # #      col = "dark blue",
  # #      main = "Bivariate Normal with Confidence Intervals")
  # # lines(ellipse(rho), col="red")       # ellipse() from ellipse package
  # # lines(ellipse(rho, level = .99), col="green")
  # # lines(ellipse(rho, level = .90), col="blue")
  # # abline(y_on_x)
  # # abline(x_on_y, col="brown")
  # # legend(3,1,legend=plot_legend,cex = .5, bty = "n")
  # 
  # # Three dimensional surface
  # # Basic perspective plot
  # persp(Ys.kde, phi = 45, theta = 30, shade = .1, border = NA) # from base graphics package
  # 
  # # RGL interactive plot
  # library(rgl)
  # col2 <- heat.colors(length(Ys.kde$z))[rank(Ys.kde$z)]
  # persp3d(x=Ys.kde, col = col2)
  # 
  # # threejs Javascript plot
  # library(threejs)
  # # Unpack data from kde grid format
  # x <- Ys.kde$x; y <- Ys.kde$y; z <- Ys.kde$z
  # # Construct x,y,z coordinates
  # xx <- rep(x,times=length(y))
  # yy <- rep(y,each=length(x))
  # zz <- z; dim(zz) <- NULL
  # # Set up color range
  # ra <- ceiling(16 * zz/max(zz))
  # col <- rainbow(16, 2/3)
  # # 3D interactive scatter plot
  # scatterplot3js(x=xx,y=yy,z=zz,size=0.4,color = col[ra],bg="black")
  # ################ Visualize Yp #################
  # # Calculate kernel density estimate
  # Yp.kde <- kde2d(Yp[,1], Yp[,2], n = 50)   # from MASS package
  # 
  # # Contour plot overlayed on heat map image of results
  # # image(Yp.kde)       # from base graphics package
  # # contour(Yp.kde, add = TRUE)     # from base graphics package
  # 
  # # Classic Bivariate Normal Diagram
  # library(ellipse)
  # # rho <- cor(Yp)
  # # y_on_x <- lm(Yp[,2] ~ Yp[,1])    # Regressiion Y ~ X
  # # x_on_y <- lm(Yp[,1] ~ Yp[,2])    # Regression X ~ Y
  # # plot_legend <- c("99% CI green", "95% CI red","90% CI blue",
  # #                  "Y on X black", "X on Y brown")
  # # 
  # # plot(Yp, xlab = "X", ylab = "Y",
  # #      col = "dark blue",
  # #      main = "Bivariate Normal with Confidence Intervals")
  # # lines(ellipse(rho), col="red")       # ellipse() from ellipse package
  # # lines(ellipse(rho, level = .99), col="green")
  # # lines(ellipse(rho, level = .90), col="blue")
  # # abline(y_on_x)
  # # abline(x_on_y, col="brown")
  # # legend(3,1,legend=plot_legend,cex = .5, bty = "n")
  # 
  # # Three dimensional surface
  # # Basic perspective plot
  # # persp(Yp.kde, phi = 45, theta = 30, shade = .1, border = NA) # from base graphics package
  # 
  # # RGL interactive plot
  # library(rgl)
  # col2 <- heat.colors(length(Yp.kde$z))[rank(Yp.kde$z)]
  # persp3d(x=Yp.kde, col = col2)
  # 
  # # threejs Javascript plot
  # library(threejs)
  # # Unpack data from kde grid format
  # x <- Yp.kde$x; y <- Yp.kde$y; z <- Yp.kde$z
  # # Construct x,y,z coordinates
  # xx <- rep(x,times=length(y))
  # yy <- rep(y,each=length(x))
  # zz <- z; dim(zz) <- NULL
  # # Set up color range
  # ra <- ceiling(16 * zz/max(zz))
  # col <- rainbow(16, 2/3)
  # # 3D interactive scatter plot
  # scatterplot3js(x=xx,y=yy,z=zz,size=0.4,color = col[ra],bg="black")
  # 
  # # Draw from multi-t distribution without truncation
  # Yp.kde <- kde2d(Yp[,1], Yp[,2], n = 50)   # from MASS package
  # col2 <- heat.colors(length(Yp.kde$z))[rank(Yp.kde$z)]
  # persp3d(x=Yp.kde, col = col2)
  ################ Estimation#################
  # #Estimated Mean
  # mu_hats = colMeans(Ys)
  # mu_hatp = colMeans(Yp)
  # #Estimated Sigma
  # sigma_hats = var(Ys)
  # sigma_hatp = var(Yp)
  # #Standardize to get residual (2 x T)
  # theta_hats = c(mu_hats, c(sigma_hats))
  # Xs = g(theta_hats, t(Ys))
  # theta_hatp = c(mu_hatp, c(sigma_hatp))
  # Xp = g(theta_hatp, t(Yp))
  ################### KSD Nres Test ###################
  pval_s_KSD = KSD_Nres(Ys)
  pval_p_KSD = KSD_Nres(Yp)
  
  # plist_size_KSD = c(plist_size_KSD, pval_s_KSD)
  # plist_power_KSD = c(plist_power_KSD, pval_p_KSD)
  ################### return ###################
  return(c(pval_s_KSD, pval_p_KSD))
  # pb$tick()
  # Sys.sleep(1 / 100)
}

################### Size & Power Evaluation ###################
# ---------- Use foreach and doSNOW package ---------- #
# n = 500
# for (n in c(100, 200, 500)){
# for (lambda in seq(0.2,1,0.2)){
n = 500
lambda = 1
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
                               'sn')) %dopar% func(cnt)
# {
#                                  tryCatch(func(cnt),error = function(e) return(NA))
#                                }

close(pb)
stopCluster(cl)
# result
cat("n =",n, "lambda =",lambda)
cbind(apply(result[,1:2]<0.01, 2, sum) / nrep,
      apply(result[,1:2]<0.05, 2, sum) / nrep,
      apply(result[,1:2]<0.10, 2, sum) / nrep)


# apply(result[,29:34], 2, sum) / nrep
# which(is.na(result))
