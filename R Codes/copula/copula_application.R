
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
find_median_distance <- function(Z){
  
  if(is.data.frame(Z)){
    Z = data.matrix(Z)
  }else{
    Z = as.array(Z)
  }
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]
  
  # if size of Z is greater than 100, randomly sample 100 points
  # if(size1 > 100){
  #   if(is.na(size2)){
  #     Zmed <- Z[sample(size1,100)]
  #   }else{
  #     Zmed <- Z[sample(size1,100),]
  #   }
  #   size1 = 100
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


################### Data ########################
# Import data
data_cop_total = read.csv("C:/Users/luodh/Desktop/Thesis/data/data_cop_total.csv", header = TRUE)
data_cop1 = read.csv("C:/Users/luodh/Desktop/Thesis/data/data_cop1.csv", header = TRUE)
data_cop2 = read.csv("C:/Users/luodh/Desktop/Thesis/data/data_cop2.csv", header = TRUE)
data_cop3 = read.csv("C:/Users/luodh/Desktop/Thesis/data/data_cop3.csv", header = TRUE)

data0 = data_cop3

GBP = data0$GBPUSD
EUR = data0$EURUSD
N0 = length(GBP)
d = 2
# Daily log differences in percentages
lr_GBP = (log(GBP[2:N0]) - log(GBP[1:(N0-1)])) * 100
lr_EUR = (log(EUR[2:N0]) - log(EUR[1:(N0-1)])) * 100

data = rbind(lr_GBP, lr_EUR)

################### Plotting #################
# library(lubridate)
# library(scales)
# library(ggplot2)
# 
# date0 = as.character(data0$Date)
# date0 = as.Date(date0, "%Y/%m/%d")
# plot.data0 = data.frame(date0, data0$GBPUSD, data0$EURUSD)
# p = ggplot(data = plot.data0, aes(x = date0)) +
#   geom_line(aes(y = data0.GBPUSD), color = 'red', size = 0.5) +
#   geom_line(aes(y = data0.EURUSD), color = 'blue', size = 0.5) +
#   scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
#   xlab('Year') + ylab('Exchange Rate(/USD)') + ggtitle('Daily Exchange Rate') +
#   theme(plot.title = element_text(hjust = 0.5)) 
# p
# ggsave('plot_cop1.eps',p)
# 
# date = as.character(data0$Date[2:N0])
# date = as.Date(date, "%Y/%m/%d")
# plot.data = data.frame(date, lr_GBP, lr_EUR)
# 
# par(mfcol = c(2,1))
# p1 = ggplot(data = plot.data, aes(x = date, y = lr_GBP)) +
#   geom_line(color = "black", size = 0.5) +
#   scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
#   xlab('Year') + ylab('log differences') + ggtitle('(A) GBP/USD')
# 
# p2 = ggplot(data = plot.data, aes(x = date, y = lr_EUR)) +
#   geom_line(color = "black", size = 0.5) +
#   scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
#   xlab('Year') + ylab('log differences') + ggtitle('(B) EUR/USD')
# 
# p1
# p2
# 
# library(cowplot)
# p3 <- cowplot::plot_grid(p1, p2, nrow = 2)#将p1-p4四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）
# p3
# ggsave('plot_cop2.eps',p3)
################### Summary Statistics #################
# mu_hat = rowMeans(data)
# sigma_hat = c(sd(lr_GBP),sd(lr_EUR))
# # rho_hat = cor(t(data))
# c(skewness(lr_GBP),skewness(lr_EUR))
# c(kurtosis(lr_GBP,method = 'moment'),kurtosis(lr_EUR,method = 'moment'))
# c(max(lr_GBP), max(lr_EUR))
# c(min(lr_GBP), min(lr_EUR))
# # data0[which(lr_EUR == max(lr_EUR)),]
# # ljung.box.test(lr_SP5)

################## Preliminary testing ###################
# Augmented Dickey-Fuller test
library(tseries)
adf.test(lr_GBP)
adf.test(lr_EUR)

# Engle's ARCH test
# library(aTSA)
# arch.test(lr_GBP)
library(FinTS)
ArchTest(lr_GBP,lags = 20)
ArchTest(lr_EUR,lags = 20)

# Ljung-Box test
Box.test(lr_GBP, lag = 20, type = "Ljung-Box")
Box.test(lr_EUR, lag = 20, type = "Ljung-Box")
################## model fitting ##################
library(rugarch)
library(fGarch)
fit_GBP = garchFit(formula = ~arma(1,0) + garch(1,1), data = lr_GBP,
         cond.dist = 'QMLE', 
         algorithm = 'lbfgsb',
         trace = FALSE)

fit_EUR = garchFit(formula = ~garch(1,1), data = lr_EUR,
                   cond.dist = 'QMLE', 
                   algorithm = 'lbfgsb',
                   trace = FALSE)
X = rbind(residuals(fit_GBP, standardize = T), residuals(fit_EUR, standardize = T))
n = N0 - 1
coef_GBP = coef(fit_GBP)
coef_EUR = coef(fit_EUR)


######### estimate marginal t dist. d.f. parameter nu #########

nu_hat_GBP = fitdist(x = X[1,], distribution = "std")$pars[["shape"]]
nu_hat_EUR = fitdist(x = X[2,], distribution = "std")$pars[["shape"]]

# test for normal marginal
shapiro.test(X[1,])
shapiro.test(X[2,])

# test for t marginal
library(goftest)
dist.t_GBP = function(q){
  return(pt(q, df = nu_hat_GBP))
}
dist.t_EUR = function(q){
  return(pt(q, df = nu_hat_EUR))
}
ad.test(X[1,],dist.t_GBP)
ad.test(X[2,],dist.t_EUR)
cvm.test(X[1,],dist.t_GBP)
cvm.test(X[2,],dist.t_EUR)

########### estimate rho_hat & nu_hat ###########
F1hat = function(x){
  return(1 / (n + 1) * sum(X[1,] <= x))
}
F2hat = function(x){
  return(1 / (n + 1) * sum(X[2,] <= x))
}
Fn1 = apply(as.array(X[1,]), 1, F1hat)
Fn2 = apply(as.array(X[2,]), 1, F2hat)
Fn = cbind(Fn1, Fn2)
# t copula
fitn = fitCopula(tCopula(dim = 2,dispstr = "un", df.fixed = F), Fn, method = 'ml')
rho_hat = fitn@estimate[1]
nu_hat = fitn@estimate[2]
# normal copula
fit_normal = fitCopula(normalCopula(dim = 2, dispstr = 'un'), Fn, method = 'ml')
rho_hat = fit_normal@estimate

############## score functions ################
# normal copula with normal margins
sqx = function(x){
  rho = rho_hat
  R = matrix(c(1, rho,
               rho, 1), 2, 2)
  return(-(solve(R) ) %*% x)
}

# normal copula with t margins
sqx = function(x){
  rho = rho_hat
  nu_1 = nu_hat_GBP
  nu_2 = nu_hat_EUR
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

# t copula with normal margins
sqx = function(x){
  if (pnorm(x[1]) == 1){
    x[1] = 8
  }
  if (pnorm(x[2]) == 1){
    x[2] = 8
  }
  rho = rho_hat
  nu = nu_hat
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

# t copula with t margins
sqx = function(x){
  rho = rho_hat
  nu = nu_hat
  nu_1 = nu_hat_GBP
  nu_2 = nu_hat_EUR
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

################ KSD test ####################
T_obs = KSD_new(x = t(X), score_function = t(apply(t(X), 1, sqx)), width = -1)
T_obs = as.numeric(T_obs)

# Distribution of T
ncut = 1000

# normal copula + AR(1)-GARCH(1,1) + normal margins
DGP_b = function(coef1, coef2, rho){
  n = N0 - 1
  sigma0 = matrix(c(1, rho, rho, 1), 2, 2)
  cop = rcopula.gauss(n + ncut, Sigma = sigma0)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0, n + ncut)
  sigma2_2 = rep(0, n + ncut)

  # mu = c(coef1[1], coef2[1])
  # phi = c(0, 0)
  # omega = c(coef1[2], coef2[2])
  # alpha = c(coef1[3], coef2[3])
  # beta  = c(coef1[4], coef2[4])

  # mu = c(coef1[1], coef2[1])
  # phi = c(0, coef2[2])
  # omega = c(coef1[2], coef2[3])
  # alpha = c(coef1[3], coef2[4])
  # beta  = c(coef1[4], coef2[5])
  
  mu = c(coef1[1], coef2[1])
  phi = c(coef1[2], 0)
  omega = c(coef1[3], coef2[2])
  alpha = c(coef1[4], coef2[3])
  beta  = c(coef1[5], coef2[4])

  for (t in 2:(n + ncut)){
    sigma2_1[t] = omega[1] + beta[1] * sigma2_1[t-1] +
      alpha[1] * sigma2_1[t-1] * qnorm(p = cop[t-1,1])^2
    sigma2_2[t] = omega[2] + beta[2] * sigma2_2[t-1] +
      alpha[2] * sigma2_2[t-1] * qnorm(p = cop[t-1,2])^2
    Y1[t] = mu[1] + phi[1] * Y1[t-1] +
      sqrt(sigma2_1[t]) * qnorm(p = cop[t,1])
    Y2[t] = mu[2] + phi[2] * Y2[t-1] +
      sqrt(sigma2_2[t]) * qnorm(p = cop[t,2])
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}

# normal copula + AR(1)-GARCH(1,1) + t margins
DGP_b = function(coef1, coef2, rho){
  n = N0 - 1
  sigma = matrix(c(1, rho, rho, 1), 2, 2)
  # cop = rcopula.t(n + ncut, df = nu, Sigma = sigma)
  cop = rcopula.gauss(n + ncut, sigma)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0, n + ncut)
  sigma2_2 = rep(0, n + ncut)
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, 0)
  # omega = c(coef1[2], coef2[2])
  # alpha = c(coef1[3], coef2[3])
  # beta  = c(coef1[4], coef2[4])
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, coef2[2])
  # omega = c(coef1[2], coef2[3])
  # alpha = c(coef1[3], coef2[4])
  # beta  = c(coef1[4], coef2[5])
  
  mu = c(coef1[1], coef2[1])
  phi = c(coef1[2], 0)
  omega = c(coef1[3], coef2[2])
  alpha = c(coef1[4], coef2[3])
  beta  = c(coef1[5], coef2[4])
  
  for (t in 2:(n + ncut)){
    sigma2_1[t] = omega[1] + beta[1] * sigma2_1[t-1] +
      alpha[1] * sigma2_1[t-1] * (nu_hat_GBP - 2) / nu_hat_GBP * 
      qt(p = cop[t-1,1], df = nu_hat_GBP)^2
    sigma2_2[t] = omega[2] + beta[2] * sigma2_2[t-1] +
      alpha[2] * sigma2_2[t-1] * (nu_hat_EUR - 2) / nu_hat_EUR * 
      qt(p = cop[t-1,2], df = nu_hat_EUR)^2
    Y1[t] = mu[1] + phi[1] * Y1[t-1] + 
      sqrt(sigma2_1[t]) * sqrt((nu_hat_GBP - 2) / nu_hat_GBP) *
      qt(p = cop[t,1], df = nu_hat_GBP)
    Y2[t] = mu[2] + phi[2] * Y2[t-1] + 
      sqrt(sigma2_2[t]) * sqrt((nu_hat_EUR - 2) / nu_hat_EUR) * 
      qt(p = cop[t,2], df = nu_hat_EUR)
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}

# t copula + AR(1)-GARCH(1,1) + normal margins
DGP_b = function(coef1, coef2, rho){
  n = N0 - 1
  sigma = matrix(c(1, rho, rho, 1), 2, 2)
  # cop = rcopula.gauss(n + ncut, Sigma = sigma)
  cop = rcopula.t(n + ncut, df = nu_hat, Sigma = sigma)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0, n + ncut)
  sigma2_2 = rep(0, n + ncut)
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, 0)
  # omega = c(coef1[2], coef2[2])
  # alpha = c(coef1[3], coef2[3])
  # beta  = c(coef1[4], coef2[4])
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, coef2[2])
  # omega = c(coef1[2], coef2[3])
  # alpha = c(coef1[3], coef2[4])
  # beta  = c(coef1[4], coef2[5])
  
  mu = c(coef1[1], coef2[1])
  phi = c(coef1[2], 0)
  omega = c(coef1[3], coef2[2])
  alpha = c(coef1[4], coef2[3])
  beta  = c(coef1[5], coef2[4])
  
  for (t in 2:(n + ncut)){
    sigma2_1[t] = omega[1] + beta[1] * sigma2_1[t-1] +
      alpha[1] * sigma2_1[t-1] * qnorm(p = cop[t-1,1])^2
    sigma2_2[t] = omega[2] + beta[2] * sigma2_2[t-1] +
      alpha[2] * sigma2_2[t-1] * qnorm(p = cop[t-1,2])^2
    Y1[t] = mu[1] + phi[1] * Y1[t-1] +
      sqrt(sigma2_1[t]) * qnorm(p = cop[t,1])
    Y2[t] = mu[2] + phi[2] * Y2[t-1] +
      sqrt(sigma2_2[t]) * qnorm(p = cop[t,2])
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}

# t copula + AR(1)-GARCH(1,1) + t margins
DGP_b = function(coef1, coef2, rho){
  n = N0 - 1
  sigma = matrix(c(1, rho, rho, 1), 2, 2)
  cop = rcopula.t(n + ncut, df = nu_hat, Sigma = sigma)
  # cop = rcopula.gauss(n + ncut, sigma)
  Y1 = rep(0, n + ncut)
  Y2 = rep(0, n + ncut)
  sigma2_1 = rep(0, n + ncut)
  sigma2_2 = rep(0, n + ncut)
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, 0)
  # omega = c(coef1[2], coef2[2])
  # alpha = c(coef1[3], coef2[3])
  # beta  = c(coef1[4], coef2[4])
  
  # mu = c(coef1[1], coef2[1])
  # phi = c(0, coef2[2])
  # omega = c(coef1[2], coef2[3])
  # alpha = c(coef1[3], coef2[4])
  # beta  = c(coef1[4], coef2[5])
  
  mu = c(coef1[1], coef2[1])
  phi = c(coef1[2], 0)
  omega = c(coef1[3], coef2[2])
  alpha = c(coef1[4], coef2[3])
  beta  = c(coef1[5], coef2[4])
  
  for (t in 2:(n + ncut)){
    sigma2_1[t] = omega[1] + beta[1] * sigma2_1[t-1] +
      alpha[1] * sigma2_1[t-1] * (nu_hat_GBP - 2) / nu_hat_GBP * 
      qt(p = cop[t-1,1], df = nu_hat_GBP)^2
    sigma2_2[t] = omega[2] + beta[2] * sigma2_2[t-1] +
      alpha[2] * sigma2_2[t-1] * (nu_hat_EUR - 2) / nu_hat_EUR * 
      qt(p = cop[t-1,2], df = nu_hat_EUR)^2
    Y1[t] = mu[1] + phi[1] * Y1[t-1] + 
      sqrt(sigma2_1[t]) * sqrt((nu_hat_GBP - 2) / nu_hat_GBP) *
      qt(p = cop[t,1], df = nu_hat_GBP)
    Y2[t] = mu[2] + phi[2] * Y2[t-1] + 
      sqrt(sigma2_2[t]) * sqrt((nu_hat_EUR - 2) / nu_hat_EUR) * 
      qt(p = cop[t,2], df = nu_hat_EUR)
  }
  Y = cbind(Y1[(ncut+1):(n+ncut)], Y2[(ncut+1):(n+ncut)])
  return(Y)
}


sample = NULL
nB = 1000
for (b in 1:nB){
  # Generate residual under H0
  Y_b = DGP_b(coef_GBP, coef_EUR, rho_hat)
  
  fit_GBP_b = garchFit(formula = ~arma(1,0) + garch(1,1), data = Y_b[,1],
                     cond.dist = 'QMLE', 
                     algorithm = 'lbfgsb',
                     trace = FALSE)
  
  fit_EUR_b = garchFit(formula = ~garch(1,1), data = Y_b[,2],
                     cond.dist = 'QMLE', 
                     algorithm = 'lbfgsb',
                     trace = FALSE)
  X_b = rbind(residuals(fit_GBP_b, standardize = T), residuals(fit_EUR_b, standardize = T))
  
  rho_hat_b = rho_hat
  
  # F1hat_b = function(x){
  #   return(1 / (n + 1) * sum(X_b[1,] <= x))
  # }
  # F2hat_b = function(x){
  #   return(1 / (n + 1) * sum(X_b[2,] <= x))
  # }
  # Fn1_b = apply(as.array(X_b[1,]), 1, F1hat_b)
  # Fn2_b = apply(as.array(X_b[2,]), 1, F2hat_b)
  # Fn_b = cbind(Fn1_b, Fn2_b)
  # t copula
  # fitn = fitCopula(tCopula(dim = 2,dispstr = "un", df.fixed = F), Fn, method = 'ml')
  # rho_hat = fitn@estimate[1]
  # nu_hat = fitn@estimate[2]
  # normal copula
  # fit_normal_b = fitCopula(normalCopula(dim = 2, dispstr = 'un'), Fn_b, method = 'ml')
  # rho_hat_b = fit_normal_b@estimate
  
  # normal copula with normal margins
  # sqx_b = function(x){
  #   rho = rho_hat_b
  #   R = matrix(c(1, rho,
  #                rho, 1), 2, 2)
  #   return(-(solve(R) ) %*% x)
  # }
  
  # normal copula with t margins
  # sqx_b = function(x){
  #   rho = rho_hat_b
  #   nu_1 = nu_hat_GBP
  #   nu_2 = nu_hat_EUR
  #   Ax1 = qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1))
  #   Bx2 = qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2))
  #   dAx1 = sqrt(nu_1 / (nu_1 - 2)) * dt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1) /
  #     dnorm(qnorm(pt(sqrt(nu_1 / (nu_1 - 2)) * x[1], nu_1)))
  #   dBx2 = sqrt(nu_2 / (nu_2 - 2)) * dt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2) /
  #     dnorm(qnorm(pt(sqrt(nu_2 / (nu_2 - 2)) * x[2], nu_2)))
  #   sqx1 = - rho^2 / (1 - rho^2) * Ax1 * dAx1 +
  #     rho / (1 - rho^2) * dAx1 * Bx2 -
  #     (nu_1 + 1) * nu_1 / (nu_1 - 2) * x[1] / (nu_1 + nu_1 / (nu_1 - 2) * x[1]^2)
  #   sqx2 = - rho^2 / (1 - rho^2) * Bx2 * dBx2 +
  #     rho / (1 - rho^2) * Ax1 * dBx2 -
  #     (nu_2 + 1) * nu_2 / (nu_2 - 2) * x[2] / (nu_2 + nu_2 / (nu_2 - 2) * x[2]^2)
  #   return(c(sqx1, sqx2))
  # }
  
  # t copula with normal margins
  # sqx_b = function(x){
  #   if (pnorm(x[1]) == 1){
  #     x[1] = 8
  #   }
  #   if (pnorm(x[2]) == 1){
  #     x[2] = 8
  #   }
  #   rho = rho_hat_b
  #   nu = nu_hat
  #   A1 = qt(pnorm(x[1]), nu)
  #   A2 = qt(pnorm(x[2]), nu)
  #   dA1 = dnorm(x[1]) / dt(qt(pnorm(x[1]), nu), nu)
  #   dA2 = dnorm(x[2]) / dt(qt(pnorm(x[2]), nu), nu)
  #   sqx1 = - (nu + d) / 2 * (2 * A1 * dA1 - 2 * rho * dA1 * A2) /
  #     (A1^2 - 2 * rho * A1 * A2 + A2^2 + nu * (1 - rho^2)) +
  #     (nu + 1) * (A1 * dA1) / (nu + A1^2) - x[1]
  #   sqx2 = - (nu + d) / 2 * (2 * A2 * dA2 - 2 * rho * dA2 * A1) /
  #     (A2^2 - 2 * rho * A2 * A1 + A1^2 + nu * (1 - rho^2)) +
  #     (nu + 1) * (A2 * dA2) / (nu + A2^2) - x[2]
  #   return(c(sqx1, sqx2))
  # }
  
  # t copula with t margins
  sqx_b = function(x){
    rho = rho_hat_b
    nu = nu_hat
    nu_1 = nu_hat_GBP
    nu_2 = nu_hat_EUR
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

  # KSD Test Statistic for bootstrap data
  T_b = KSD_new(x = t(X_b), score_function = t(apply(t(X_b), 1, sqx_b)), width = -1)
  T_b = as.numeric(T_b)
  sample = c(sample, T_b)
}
pval_KSD = sum((sample - T_obs) > 0) / nB
pval_KSD

