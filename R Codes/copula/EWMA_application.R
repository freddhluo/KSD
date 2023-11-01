# load dependencies from EWMA_simulation.R

######################################################
# Import data
# data_AAPL = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/AAPL.csv", header = TRUE)
# data_CSCO = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/CSCO.csv", header = TRUE)
# data_INTC = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/INTC.csv", header = TRUE)
# data_MSFT = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/MSFT.csv", header = TRUE)
# 
# data_AA = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/AA.csv", header = TRUE)
# data_CAT = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/CAT.csv", header = TRUE)
# data_DIS = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/DIS.csv", header = TRUE)
# 
# data_ESXB = read.csv("C:/Users/luodh/Desktop/Thesis/data/ESXB Historical Data.csv", header = TRUE)
# data_MCBF = read.csv("C:/Users/luodh/Desktop/Thesis/data/MCBF Historical Data.csv", header = TRUE)

# ESXB = data_ESXB$Price
# MCBF = data_MCBF$Price
# N0 = length(ESXB)
# d = 2
# data_SP500 = read.csv("C:/Users/luodh/OneDrive - connect.hku.hk/Time Series Research/KSD paper/KSD Real Data Application/SP5001.csv", header = TRUE)
# SP500 = data_SP500$Adj.Close
# N0 = length(SP500)
# lr_SP500 = (log(SP500[2:N0]) - log(SP500[1:(N0-1)])) * 100

# mean(lr_SP500)
# var(lr_SP500)
# plot(1:(N0-1), lr_SP500, 'l')
# plot(1:551, lr_SP500[1750:2300], 'l')
data_cop_total = read.csv("C:/Users/luodh/Desktop/Thesis/data/data_cop_total.csv", header = TRUE)
data0 = data_cop_total[911:1406,]

GBP = data0$GBPUSD
EUR = data0$EURUSD
N0 = length(GBP)
d = 2
# Daily log differences in percentages
lr_GBP = (log(GBP[2:N0]) - log(GBP[1:(N0-1)])) * 100
lr_EUR = (log(EUR[2:N0]) - log(EUR[1:(N0-1)])) * 100
# lr_ESXB = (log(ESXB[2:N0]) - log(ESXB[1:(N0-1)])) * 100
# plot(1:(N0-1), lr_ESXB, 'l')

data = rbind(lr_GBP, lr_EUR)
# plot(1:500, lr_GBP[1:500], 'l')
mu_hat = rowMeans(data)
data = data - mu_hat
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
#   scale_x_date(labels = date_format("%Y/%m"), breaks = date_breaks("4 months")) +
#   xlab('Date') + ylab('Exchange Rate(/USD)') + ggtitle('Daily Exchange Rate') +
#   theme(plot.title = element_text(hjust = 0.5))
# p
# ggsave('plot_nst1.eps',p)
# 
# date = as.character(data0$Date[2:N0])
# date = as.Date(date, "%Y/%m/%d")
# plot.data = data.frame(date, lr_GBP, lr_EUR)
# 
# par(mfcol = c(2,1))
# p1 = ggplot(data = plot.data, aes(x = date, y = lr_GBP)) +
#   geom_line(color = "black", size = 0.5) +
#   scale_x_date(labels = date_format("%Y/%m"), breaks = date_breaks("4 months")) +
#   xlab('Date') + ylab('log differences') + ggtitle('(A) GBP/USD')
# 
# p2 = ggplot(data = plot.data, aes(x = date, y = lr_EUR)) +
#   geom_line(color = "black", size = 0.5) +
#   scale_x_date(labels = date_format("%Y/%m"), breaks = date_breaks("4 months")) +
#   xlab('Date') + ylab('log differences') + ggtitle('(B) EUR/USD')
# 
# p1
# p2
# 
# library(cowplot)
# p3 <- cowplot::plot_grid(p1, p2, nrow = 2)#将p1-p4四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D。（LETTERS[1:4] 意为提取26个大写英文字母的前四个:A、B、C、D）
# p3
# ggsave('plot_nst2.eps',p3)
################### Summary Statistics #################
mu_hat = rowMeans(data)
sigma_hat = c(sd(lr_GBP),sd(lr_EUR))
# rho_hat = cor(t(data))
c(skewness(lr_GBP),skewness(lr_EUR))
c(kurtosis(lr_GBP,method = 'moment'),kurtosis(lr_EUR,method = 'moment'))
c(max(lr_GBP), max(lr_EUR))
c(min(lr_GBP), min(lr_EUR))
################## Preliminary testing ###################
# Augmented Dickey-Fuller test
library(tseries)
# adf.test(lr_ESXB)
adf.test(lr_GBP)
adf.test(lr_EUR)

kpss.test(lr_GBP,null = 'Trend')
kpss.test(lr_EUR,null = 'Trend')
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
est = eccc.estimation(a = c(0,0), A = A_init, B = B_init,
                      R = diag(2), dvar = t(data), model = "EWMA-cop")
X = t(est$std.resid)
A_hat = est$para.mat$A
B_hat = est$para.mat$B
n = N0 - 1


######### estimate marginal t dist. d.f. parameter nu #########
library(rugarch)
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

cop = 'N'
dist = 't'
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
  } else if (dist == 't'){
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
  }
}
################ KSD test ####################
T_obs = KSD_new(x = t(X), score_function = t(apply(t(X), 1, sqx)), width = -1)
T_obs = as.numeric(T_obs)

# Distribution of T
ncut = 1000

DGP_b = function(cop, dist, A., B., rho.){
  nu = nu_hat
  nu_1 = nu_hat_GBP
  nu_2 = nu_hat_EUR
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



sample = NULL
nB = 1000
for (b in 1:nB){
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
    } else if (dist == 't'){
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
    }
  }
  
  # KSD Test Statistic for bootstrap data
  T_b = KSD_new(x = t(X_b), score_function = t(apply(t(X_b), 1, sqx_b)), width = -1)
  T_b = as.numeric(T_b)
  sample = c(sample, T_b)
}

pval_KSD = sum((sample - T_obs) > 0) / nB
cat("Marginal:", dist, 
    "| H0 Copula:", cop)
pval_KSD





