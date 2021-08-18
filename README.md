# Semi-parametric Bayesian Additive Regression trees
This repository houses R scripts and data sets that can be used to reproduce the results presented in the paper [_Semi-parametric Bayesian Additive Regression Trees_. arXiv (2021)](https://arxiv.org/abs/2108.07636).

In addition, it provides an implementation of SP-BART in the format of an R package named ```spbart```.

## Installation
``` r
library(devtools)
install_github("ebprado/SP-BART/spbart",ref='main')
```
## Example
``` r
library(spbart)
rm(list = ls())

# ---------------------------------
# SP-BART for a continuous response
# ---------------------------------

# Simulate from Friedman equation
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = as.data.frame(x)))
}

n = 200
ncov = 5
var = 1
data = friedman_data(n, ncov, sqrt(var))

y = data$y # response variable
X1 = as.data.frame(cbind(y,data$x)) # all covariates + the response
X2 = data$x # all covariates

# Run the semi-parametric BART (WITHOUT intercept)
spbart.fit = spbart::semibart(formula = y ~ 0 + V4 + V5, x1 = X1, x2 = X2, ntrees = 10, nburn = 2000, npost = 1000)

# Run the semi-parametric BART (WITH intercept)
# spbart.fit = spbart::semibart(formula = y ~ V4 + V5, x1 = X1, x2 = X2, ntrees = 10, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat)
yhat = apply(spbart.fit$y_hat,2,mean)
betahat = apply(spbart.fit$beta_hat,2,mean)

# Predict on a new dataset
yhat_pred = spbart::predict_semibart(spbart.fit, newdata_x1 = X1, newdata_x2 = X2, type = 'mean')
cor(yhat,yhat_pred) == 1

# Plot 
plot(y, yhat);abline(0,1)
plot(1:2, c(10,5), main = 'True versus estimates', ylim=c(3,12))
points(1:2, betahat, col=2, pch=2)

# -----------------------------
# SP-BART for a binary response
# -----------------------------
n = 200
ncov = 5
var = 1
data = friedman_data(n, ncov, sqrt(var))

aux = data$y
y = ifelse(aux > median(aux), 1, 0)
X1 = as.data.frame(cbind(y,data$x)) # all covariates + the response
X2 = data$x # all covariates

# Run the semi-parametric BART (WITH intercept)
spbart.fit = spbart::cl_semibart(formula = y ~ V4 + V5, x1 = X1, x2 = X2, ntrees = 1, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat)
yhat = pnorm(apply(spbart.fit$y_hat,2,mean))
betahat = apply(spbart.fit$beta_hat,2,mean)

# Predict on a new dataset
yhat_pred = spbart::cl_predict_semibart(spbart.fit, newdata_x1 = X1, newdata_x2 = X2, type = 'mean')
cor(yhat,yhat_pred) == 1
```
