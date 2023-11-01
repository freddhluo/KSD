# Simulations for copula-EWMA model

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

################### Load KSD Dependencies ################### 
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

# g(theta, Y)
# Input: theta (1 x p), Y: original data matrix (d x n)
# Output: value of g (d x n)
g = function(theta, Y){
  mu = repmat(theta[1:d], 1, n)
  # sigma_sqrtinv = matrix(theta[3:6],nrow = 2, ncol = 2)
  sigma = diag(c(theta[d+1], theta[d+d^2]))
  sigma_sqrtinv = msqrt(sigma)$invsqrt
  return(sigma_sqrtinv %*% (Y - mu))
}

# score function for t(0, (nu-2)/nu*I, nu) distribution
# Input: data x (n x d); d.f. nu; dimension d
# Output: Sq(x) (n x d)
# sqx = function(x){
#   return(t(apply(x, 1, function(x) -(nu + d) / (nu - 2) * 1 / (1 + 1 / (nu - 2) * sum(x^2)) * x)))
# }

# check
# temp = matrix(1:6, 3, 2, byrow = T)
# sqx(temp)
# -(nu + d) * temp[1,] / (nu - 2 + sum((temp[1,])^2))

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
KSD_Nres = function(Y, cop, dist, g. = g, KSD_new. = KSD_new, nB. = nB, n. = n){
  ################ Estimation #################
  est = eccc.estimation(a = c(0,0), A = A_init, B = B_init,
                          R = diag(2), dvar = Y, model = "EWMA-cop")
  X = t(est$std.resid)
  A_hat = est$para.mat$A
  B_hat = est$para.mat$B
  
  n0 = n.
  rho_hat = rho
  
  # Set score function
  if (cop == 'N'){
    if (dist == 'N'){
      # normal copula with normal margins
      sqx = function(x){
          rho = rho_hat
          R = matrix(c(1, rho,
                       rho, 1), 2, 2)
          return(-(solve(R) ) %*% x)
        } 
    } else if (dist == 't'){
      # normal copula with t margins
      sqx = function(x){
        rho = rho_hat
        # nu_1 = nu_hat_GBP
        # nu_2 = nu_hat_EUR
        Ax1 = qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1))
        Bx2 = qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2))
        dAx1 = sqrt(nu_1 / (nu_1 - 2)) * dt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1) /
          dnorm(qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1)))
        dBx2 = sqrt(nu_2 / (nu_2 - 2)) * dt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2) /
          dnorm(qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2)))
        sqx1 = - rho^2 / (1 - rho^2) * Ax1 * dAx1 +
          rho / (1 - rho^2) * dAx1 * Bx2 -
          (nu_1 + 1) * nu_1 / (nu_1 - 2) * x[1] / (nu_1 + nu_1 / (nu_1 - 2) * x[1]^2)
        sqx2 = - rho^2 / (1 - rho^2) * Bx2 * dBx2 +
          rho / (1 - rho^2) * Ax1 * dBx2 -
          (nu_2 + 1) * nu_2 / (nu_2 - 2) * x[2] / (nu_2 + nu_2 / (nu_2 - 2) * x[2]^2)
        return(c(sqx1, sqx2))
      }
    }
  } else if (cop == 't'){
    if (dist == 'N'){
      # t copula with normal margins
      sqx = function(x){
        if (pnorm(x[1]) == 1){
          x[1] = 8
        }
        if (pnorm(x[2]) == 1){
          x[2] = 8
        }
        rho = rho_hat
        # nu = nu_hat
        A1 = qt(pnorm(x[1]), nu)
        A2 = qt(pnorm(x[2]), nu)
        dA1 = dnorm(x[1]) / dt(qt(pnorm(x[1]), nu), nu)
        dA2 = dnorm(x[2]) / dt(qt(pnorm(x[2]), nu), nu)
        sqx1 = - (nu + d) / 2 * (2 * A1 * dA1 - 2 * rho * dA1 * A2) /
          (A1^2 - 2 * rho * A1 * A2 + A2^2 + nu * (1 - rho^2)) +
          (nu + 1) * (A1 * dA1) / (nu + A1^2) - x[1]
        sqx2 = - (nu + d) / 2 * (2 * A2 * dA2 - 2 * rho * dA2 * A1) /
          (A2^2 - 2 * rho * A2 * A1 + A1^2 + nu * (1 - rho^2)) +
          (nu + 1) * (A2 * dA2) / (nu + A2^2) - x[2]
        return(c(sqx1, sqx2))
      }
    } else if (dist == 't'){
      # t copula with t margins
      sqx = function(x){
        rho = rho_hat
        # nu = nu_hat
        # nu_1 = nu_hat_GBP
        # nu_2 = nu_hat_EUR
        b1 = sqrt(nu_1 / (nu_1 - 2))
        b2 = sqrt(nu_2 / (nu_2 - 2))
        A1 = qt(pt(b1 * x[1], nu_1), nu)
        A2 = qt(pt(b2 * x[2], nu_2), nu)
        dA1 = b1 * dt(b1 * x[1], nu_1) /
          dt(qt(pt(b1 * x[1], nu_1), nu), nu)
        dA2 = b2 * dt(b2 * x[2], nu_2) /
          dt(qt(pt(b2 * x[2], nu_2), nu), nu)
        sqx1 = - (nu + d) / 2 * (2 * A1 * dA1 - 2 * rho * dA1 * A2) /
          (A1^2 - 2 * rho * A1 * A2 + A2^2 + nu * (1 - rho^2)) +
          (nu + 1) / 2 * (2 * A1 * dA1) / (nu + A1^2) - 
          (nu_1 + 1) * b1^2 * x[1] / (nu_1 + b1^2 * x[1]^2)
        sqx2 = - (nu + d) / 2 * (2 * A2 * dA2 - 2 * rho * dA2 * A1) /
          (A2^2 - 2 * rho * A2 * A1 + A1^2 + nu * (1 - rho^2)) +
          (nu + 1) / 2 * (2 * A2 * dA2) / (nu + A2^2) - 
          (nu_2 + 1) * b2^2 * x[2] / (nu_2 + b2^2 * x[2]^2)
        return(c(sqx1, sqx2))
      }
    }
  }
  
  # require(optimx)
  # opt = optimx(par = 0, fn = obj, lower = -1 + 0.00001, upper = 1 - 0.00001, method = 'L-BFGS-B')
  # temp = NULL
  # for (rho in seq(-0.99, 0.99, 0.01)){
  #   temp = rbind(temp, c(rho, obj(rho)))
  # }
  # plot(temp[,1],temp[,2])
  ################### KSD Test ###################
  # Observed Test Statistics
  T_obs = KSD_new.(x = t(X), score_function = t(apply(t(X), 1, sqx)), width = -1)
  T_obs = as.numeric(T_obs)
  
  # Distribution of T
  # sample = NULL
  for (b in 1:nB.){
    # Generate residual under H0
    Y_b = DGP_b(cop, dist, A_hat, B_hat, rho_hat)
    
    est_b = eccc.estimation(a = c(0,0), A = A_init, B = B_init,
                          R = diag(2), dvar = Y_b, model = "EWMA-cop")
    X_b = t(est_b$std.resid)
    # A_hat = est$para.mat$A
    # B_hat = est$para.mat$B
    
    rho_hat_b = rho_hat
    
    # Set score function
    if (cop == 'N'){
      if (dist == 'N'){
        # normal copula with normal margins
        sqx_b = function(x){
          rho = rho_hat_b
          R = matrix(c(1, rho,
                       rho, 1), 2, 2)
          return(-(solve(R) ) %*% x)
        } 
      } else if (dist == 't'){
        # normal copula with t margins
        sqx_b = function(x){
          rho = rho_hat_b
          # nu_1 = nu_hat_GBP
          # nu_2 = nu_hat_EUR
          Ax1 = qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1))
          Bx2 = qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2))
          dAx1 = sqrt(nu_1 / (nu_1 - 2)) * dt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1) /
            dnorm(qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1)))
          dBx2 = sqrt(nu_2 / (nu_2 - 2)) * dt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2) /
            dnorm(qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2)))
          sqx1 = - rho^2 / (1 - rho^2) * Ax1 * dAx1 +
            rho / (1 - rho^2) * dAx1 * Bx2 -
            (nu_1 + 1) * nu_1 / (nu_1 - 2) * x[1] / (nu_1 + nu_1 / (nu_1 - 2) * x[1]^2)
          sqx2 = - rho^2 / (1 - rho^2) * Bx2 * dBx2 +
            rho / (1 - rho^2) * Ax1 * dBx2 -
            (nu_2 + 1) * nu_2 / (nu_2 - 2) * x[2] / (nu_2 + nu_2 / (nu_2 - 2) * x[2]^2)
          return(c(sqx1, sqx2))
        }
      }
    } else if (cop == 't'){
      if (dist == 'N'){
        # t copula with normal margins
        sqx_b = function(x){
          if (pnorm(x[1]) == 1){
            x[1] = 8
          }
          if (pnorm(x[2]) == 1){
            x[2] = 8
          }
          rho = rho_hat_b
          # nu = nu_hat
          A1 = qt(pnorm(x[1]), nu)
          A2 = qt(pnorm(x[2]), nu)
          dA1 = dnorm(x[1]) / dt(qt(pnorm(x[1]), nu), nu)
          dA2 = dnorm(x[2]) / dt(qt(pnorm(x[2]), nu), nu)
          sqx1 = - (nu + d) / 2 * (2 * A1 * dA1 - 2 * rho * dA1 * A2) /
            (A1^2 - 2 * rho * A1 * A2 + A2^2 + nu * (1 - rho^2)) +
            (nu + 1) * (A1 * dA1) / (nu + A1^2) - x[1]
          sqx2 = - (nu + d) / 2 * (2 * A2 * dA2 - 2 * rho * dA2 * A1) /
            (A2^2 - 2 * rho * A2 * A1 + A1^2 + nu * (1 - rho^2)) +
            (nu + 1) * (A2 * dA2) / (nu + A2^2) - x[2]
          return(c(sqx1, sqx2))
        }
      } else if (dist == 't'){
        # t copula with t margins
        sqx_b = function(x){
          rho = rho_hat_b
          # nu = nu_hat
          # nu_1 = nu_hat_GBP
          # nu_2 = nu_hat_EUR
          b1 = sqrt(nu_1 / (nu_1 - 2))
          b2 = sqrt(nu_2 / (nu_2 - 2))
          A1 = qt(pt(b1 * x[1], nu_1), nu)
          A2 = qt(pt(b2 * x[2], nu_2), nu)
          dA1 = b1 * dt(b1 * x[1], nu_1) /
            dt(qt(pt(b1 * x[1], nu_1), nu), nu)
          dA2 = b2 * dt(b2 * x[2], nu_2) /
            dt(qt(pt(b2 * x[2], nu_2), nu), nu)
          sqx1 = - (nu + d) / 2 * (2 * A1 * dA1 - 2 * rho * dA1 * A2) /
            (A1^2 - 2 * rho * A1 * A2 + A2^2 + nu * (1 - rho^2)) +
            (nu + 1) / 2 * (2 * A1 * dA1) / (nu + A1^2) - 
            (nu_1 + 1) * b1^2 * x[1] / (nu_1 + b1^2 * x[1]^2)
          sqx2 = - (nu + d) / 2 * (2 * A2 * dA2 - 2 * rho * dA2 * A1) /
            (A2^2 - 2 * rho * A2 * A1 + A1^2 + nu * (1 - rho^2)) +
            (nu + 1) / 2 * (2 * A2 * dA2) / (nu + A2^2) - 
            (nu_2 + 1) * b2^2 * x[2] / (nu_2 + b2^2 * x[2]^2)
          return(c(sqx1, sqx2))
        }
      }
    }
    
    # KSD Test Statistic for bootstrap data
    T_b = KSD_new.(x = t(X_b), score_function = t(apply(t(X_b), 1, sqx_b)), width = -1)
    T_b = as.numeric(T_b)
    # sample = c(sample, T_b)
  }
  
  # pval_KSD = sum((sample - T_obs) > 0) / nB.
  return(c(T_obs,T_b))
}

########## load dependencies for eccc-estimation ##########
library(ccgarch)
# esitimating an (E)CCC-GARCH(1,1) model
eccc.estimation <- function(a, A, B, R, dvar, model, method="BFGS"){
  dvar <- as.matrix(dvar)
  nobs <- dim(dvar)[1]
  ndim <- dim(dvar)[2]
  if(!is.matrix(A)||!is.matrix(B)){
    stop("A or B or both must be matrices")
  }
  if(model=="diagonal"){
    init <- c(sqrt(a), diag(sqrt(A)), diag(sqrt(B)), R[lower.tri(R)])
  } else if (model == 'EWMA-cop'){
    init <- c(as.vector(sqrt(A)), sqrt(B[1,1]), sqrt(B[2,2]))
  }
  results <- optim(par=init, fn=loglik.eccc, method=method, dvar=dvar, model=model, control=list(maxit=10^5, reltol=1e-15))
  
  if(results$convergence != 0){
    cat("***********************************************************\n")
    cat("* Optimization is FAILED.                                 *\n")
    stop("* Fine tuning is required.                                *\n")
    cat("***********************************************************\n")
  }
  
  estimates <- p.mat(results$par, model=model, ndim=ndim)
  esta <- c(0,0)
  estA <- estimates$A
  estB <- estimates$B
  estR <- diag(2)
  
  h <- vector.garch(dvar, esta, estA, estB)    # estimated volatility
  std.resid <- dvar/sqrt(h)                          # std. residuals
  
  # grad <- analytical.grad(a=esta, A=estA, B=estB, R=estR, u=dvar, model=model)
  # grad <- grad%*%t(grad)/nobs
  # H    <- analytical.Hessian(a=esta, A=estA, B=estB, R = estR, u=dvar, model=model)/nobs
  # invH <- solve(H)
  # Avar <- (invH%*%grad%*%invH)/nobs
  # rob.se <- sqrt(diag(Avar))                # robust standard errors
  # out.se <- sqrt(diag(solve(grad)/nobs))                # standard errors based on the outer product of the gradient
  # H.se   <- sqrt(diag(invH/nobs))                # standard errors based on the inverted Hessian
  
  # std.error <- rbind(H.se, out.se, rob.se); rownames(std.error) <- c("Inv. Hessian", "Outer Prod.", "Robust")
  
  name.a <- numeric(0)
  name.A <- numeric(0)
  name.B <- numeric(0)
  name.R <- numeric(0)
  for(i in 1:ndim){          # index for column
    name.a <- c(name.a, paste("a", paste(i), sep=""))
    for(j in 1:ndim){    # index for row
      name.A <- c(name.A, paste("A",paste(paste(j), i, sep=""), sep=""))
      name.B <- c(name.B, paste("B",paste(paste(j), i, sep=""), sep=""))
      name.R <- c(name.R, paste("R",paste(paste(j), i, sep=""), sep=""))
    }
  }
  names(esta) <- name.a
  if(model=="diagonal"){
    ind <- as.vector(diag(ndim))
    vecA <- diag(estA);        name.A <- name.A[ind==1]
    vecB <- diag(estB);        name.B <- name.B[ind==1]
  } else if (model == 'EWMA-cop'){
    vecA <- as.vector(estA)
    vecB <- as.vector(estB)
  }
  names(vecA) <- name.A
  names(vecB) <- name.B
  vecR <- estR[lower.tri(estR)];         name.R <- matrix(name.R, ndim, ndim);  name.R <- name.R[lower.tri(name.R)]
  names(vecR) <- name.R
  # colnames(std.error) <- c(name.a, name.A, name.B, name.R)
  para.estimates <- c(esta, vecA, vecB, vecR)
  
  list(out=para.estimates, h=h, std.resid=std.resid, opt=results, para.mat=estimates)
}

# computing a likelihood value for the (E)CCC-GARCH(1,1) model
loglik.eccc <- function(param, dvar, model){
  nobs <- dim(dvar)[1]
  ndim <- dim(dvar)[2]
  para.mat <- p.mat(param, model, ndim)
  a <- para.mat$a
  A <- para.mat$A
  B <- para.mat$B
  R <- para.mat$R
  
  # check if R is positive definite
  eigenR <- eigen(R)$values
  if(max(abs(R[lower.tri(R)]))>1.0||min(eigenR)<0||!is.double(eigenR)){
    R <- diag(ndim)
  }
  h <- vector.garch(dvar, a, A, B)
  z <- dvar/sqrt(h)
  lndetR <- log(det(R))
  invR <- solve(R)
  lf <- -0.5*nobs*ndim*log(2*pi) - 0.5*sum(log(h)) - 0.5*nobs*lndetR - 0.5*sum((z%*%invR)*z)
  -lf
}

# constructing parameter vector and matrices. positivity constraints are imposed in "lleccc.c"
p.mat <- function(para, model, ndim){
  npara <- length(para)
  if(model=="diagonal"){                         # for the diagonal vector GARCH equation
    a <- para[1:ndim]^2                          # constant in variance
    A <- diag(para[(ndim+1):(2*ndim)]^2)               # ARCH parameter
    B <- diag(para[(2*ndim+1):(3*ndim)]^2)                   # GARCH parameter
    R <- diag(ndim)                              # Constant Conditional Correlation Matrix
    R[lower.tri(R)] <- para[(3*ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
  } else if(model=="EWMA-cop"){   # for the EWMA-copula model equation
    a <- c(0,0)
    A <- matrix(para[(1):(4)]^2, ndim, ndim)
    B <- diag(para[5:6]^2)
    R <- diag(ndim)
  } else if(model=="ECCC.neg"){   # the extended model with negative interaction
    a <- para[1:ndim]^2
    A <- matrix(para[(ndim+1):(ndim^2+ndim)]^2, ndim, ndim)
    B <- matrix(para[(ndim^2+ndim+1):(2*ndim^2+ndim)], ndim, ndim)
    R <- diag(ndim)
    R[lower.tri(R)] <- para[(2*ndim^2+ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
  }
  list(a=a, A=A, B=B, R=R)
}

#computing vector GARCH volatilities: valid for CCC, ECCC, DCC, EDCC models
vector.garch <- function(dvar, a, A, B){
  dvar <- dvar^2           # dvar = eps
  .Call("vector_garch", dvar, a, A, B)
}

###############   Settings   ###############
# Dimension of data
d = 2
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
nu_1 = 4
nu_2 = 6
# tolerance
# tol = .Machine$double.eps^0.25
# EWMA parameters
A = matrix(c(0.1, 0.1,
             0.05, 0.15),2,2)
B = matrix(c(0.85, 0,
             0, 0.8),2,2)
# Specify initial parameters for estimation
A_init <<- diag(rep(0.2 - 0.01, d)) + matrix(0.01, d, d)
B_init <<- diag(rep(0.6, d))
# B_init <<- diag(rep(0.6 - 0.01, d)) + matrix(0.01, d, d)

# copula parameter
tau = 0.3
# tau = 0.7
sigma = matrix(c(1, sin(tau * pi / 2),
                 sin(tau * pi / 2), 1), 2, 2)
rho = sin(tau * pi / 2)

# Simulate EWMA-copula Samples
DGP = function(cop, dist){
  if (cop == 't'){
    cop = rcopula.t(n + ncut, df = nu, Sigma = sigma)
  } 
  else if (cop == 'N'){
    cop = rcopula.gauss(n + ncut, sigma)
  } 
  else if (cop == 'Clayton'){
    cop = rcopula.clayton(n = n + ncut, theta = 2*tau/(1-tau), d = 2)
  } 
  else if (cop == 'Gumbel'){
    cop = rcopula.gumbel(n = n + ncut, theta = 1 / (1-tau), d = 2)
  }
  ita1 = rep(0, n + ncut)
  ita2 = rep(0, n + ncut)
  # ita = array(rep(0,d*n),dim = c(d,n))
  for (t in 1:(n+ncut)){
    if (dist == 'N'){
      ita1[t] = qnorm(p = cop[t,1])
      ita2[t] = qnorm(p = cop[t,2])
    } else if (dist == 't'){
      ita1[t] = sqrt((nu_1 - 2) / nu_1) * 
        qt(p = cop[t,1], df = nu_1)
      ita2[t] = sqrt((nu_2 - 2) / nu_2) *
        qt(p = cop[t,2], df = nu_2)
    }
  }

  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0.1, n + ncut)
  sigma2_2 = rep(0.1, n + ncut)
  
  for (t in 2:(n + ncut)){
    sigma2_1[t] = B[1,1] * sigma2_1[t-1] +
      A[1,1] * sigma2_1[t-1] * ita1[t-1]^2 +
      A[1,2] * sigma2_2[t-1] * ita2[t-1]^2
    sigma2_2[t] = B[2,2] * sigma2_2[t-1] +
      A[2,1] * sigma2_1[t-1] * ita1[t-1]^2 +
      A[2,2] * sigma2_2[t-1] * ita2[t-1]^2
    Y1[t] = sqrt(sigma2_1[t]) * ita1[t]
    Y2[t] = sqrt(sigma2_2[t]) * ita2[t]
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}

DGP_b = function(cop, dist, A., B., rho.){
  sigma = matrix(c(1, rho., rho., 1), 2, 2)
  if (cop == 't'){
    cop = rcopula.t(n + ncut, df = nu, Sigma = sigma)
  } 
  else if (cop == 'N'){
    cop = rcopula.gauss(n + ncut, sigma)
  } 
  ita1 = rep(0, n + ncut)
  ita2 = rep(0, n + ncut)
  # ita = array(rep(0,d*n),dim = c(d,n))
  for (t in 1:(n+ncut)){
    if (dist == 'N'){
      ita1[t] = qnorm(p = cop[t,1])
      ita2[t] = qnorm(p = cop[t,2])
    } else if (dist == 't'){
      ita1[t] = sqrt((nu_1 - 2) / nu_1) * 
        qt(p = cop[t,1], df = nu_1)
      ita2[t] = sqrt((nu_2 - 2) / nu_2) *
        qt(p = cop[t,2], df = nu_2)
    }
  }
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0.1, n + ncut)
  sigma2_2 = rep(0.1, n + ncut)
  
  for (t in 2:(n + ncut)){
    sigma2_1[t] = B.[1,1] * sigma2_1[t-1] +
      A.[1,1] * sigma2_1[t-1] * ita1[t-1]^2 +
      A.[1,2] * sigma2_2[t-1] * ita2[t-1]^2
    sigma2_2[t] = B.[2,2] * sigma2_2[t-1] +
      A.[2,1] * sigma2_1[t-1] * ita1[t-1]^2 +
      A.[2,2] * sigma2_2[t-1] * ita2[t-1]^2
    Y1[t] = sqrt(sigma2_1[t]) * ita1[t]
    Y2[t] = sqrt(sigma2_2[t]) * ita2[t]
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}

# temp = Ys
# par(mfcol = c(2,1))
# plot(1:n, temp[,1], 'l')
# plot(1:n, temp[,2], 'l')
# colMeans(temp)
# var(temp)
# est_temp = eccc.estimation(a = c(0,0), A = A_init, B = B_init,
#                 R = diag(2), dvar = temp, model = "EWMA-cop")
# temp_X = est_temp$std.resid
# est_temp$para.mat$A
# est_temp$para.mat$B
# plot(1:n, temp_X[,1], 'l')
# plot(1:n, temp_X[,2], 'l')
# colMeans(temp_X)
# var(temp_X)

cop_KSD_s = 't'
cop_KSD_p1 = 'N'
cop_KSD_p2 = 'Clayton'
cop_KSD_p3 = 'Gumbel'
dist_KSD = 't'
################### Execution ###################
func = function(cnt){
  sign = 0
  while (sign == 0){
    ################# Generate Data #################
    # set.seed(count + 12345)
    # Generate Sample Data
    Ys = DGP(cop_KSD_s, dist_KSD)
    Yp1 = DGP(cop_KSD_p1, dist_KSD)
    Yp2 = DGP(cop_KSD_p2, dist_KSD)
    Yp3 = DGP(cop_KSD_p3, dist_KSD)
    
    ################### KSD Nres Test ###################
    KSD_s = tryCatch(KSD_Nres(Ys, cop_KSD_s, dist_KSD),
                     error = function(e) {return(NA)},
                     warning = function(w) NA)
    if (is.na(KSD_s[1])){next}
    # KSD_s = KSD_Nres(Ys, cop_KSD_s, dist_KSD)
    Ts_obs = KSD_s[1]
    Ts_b = KSD_s[2]
    
    KSD_p1 = tryCatch(KSD_Nres(Yp1, cop_KSD_s, dist_KSD),
                     error = function(e) {return(NA)},
                     warning = function(w) NA)
    if (is.na(KSD_p1[1])){next}
    # KSD_p1 = KSD_Nres(Yp1, cop_KSD_s, dist_KSD) 
    Tp1_obs = KSD_p1[1]
    Tp1_b = KSD_p1[2]
    
    KSD_p2 = tryCatch(KSD_Nres(Yp2, cop_KSD_s, dist_KSD),
                      error = function(e) {return(NA)},
                      warning = function(w) NA)
    if (is.na(KSD_p2[1])){next}
    # KSD_p2 = KSD_Nres(Yp2, cop_KSD_s, dist_KSD) 
    Tp2_obs = KSD_p2[1]
    Tp2_b = KSD_p2[2]
    
    KSD_p3 = tryCatch(KSD_Nres(Yp3, cop_KSD_s, dist_KSD),
                      error = function(e) {return(NA)},
                      warning = function(w) NA)
    if (is.na(KSD_p3[1])){next}
    # KSD_p3 = KSD_Nres(Yp3, cop_KSD_s, dist_KSD) 
    Tp3_obs = KSD_p3[1]
    Tp3_b = KSD_p3[2]
    
  sign = 1
  }
  ################### return ###################
  return(c(Ts_obs, Tp1_obs, Tp2_obs, Tp3_obs,
           Ts_b, Tp1_b, Tp2_b, Tp3_b))
}

################# Size & Power Evaluation #################
# ----------- Use foreach and doSNOW package ----------- #
nrep = 10000
cl <- makeCluster(4)
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
                               'copula',
                               'tseries',
                               'QRM',
                               'ccgarch')) %dopar% func(cnt)


close(pb)
stopCluster(cl)

res = result
CV_s = quantile(res[,5], probs = c(0.99, 0.95, 0.90))
CV_p1 = quantile(res[,6], probs = c(0.99, 0.95, 0.90))
CV_p2 = quantile(res[,7], probs = c(0.99, 0.95, 0.90))
CV_p3 = quantile(res[,8], probs = c(0.99, 0.95, 0.90))
cat("Marginal:", dist_KSD, 
    "| H0 Copula:", cop_KSD_s,
    "| H1 Copula:", cop_KSD_p1,
    "/", cop_KSD_p2,
    "/", cop_KSD_p3)
rbind(c(sum((res[,1] - CV_s[1]) > 0) / nrep,
        sum((res[,1] - CV_s[2]) > 0) / nrep,
        sum((res[,1] - CV_s[3]) > 0) / nrep),
      c(sum((res[,2] - CV_p1[1]) > 0) / nrep,
        sum((res[,2] - CV_p1[2]) > 0) / nrep,
        sum((res[,2] - CV_p1[3]) > 0) / nrep),
      c(sum((res[,3] - CV_p2[1]) > 0) / nrep,
        sum((res[,3] - CV_p2[2]) > 0) / nrep,
        sum((res[,3] - CV_p2[3]) > 0) / nrep),
      c(sum((res[,4] - CV_p3[1]) > 0) / nrep,
        sum((res[,4] - CV_p3[2]) > 0) / nrep,
        sum((res[,4] - CV_p3[3]) > 0) / nrep))



