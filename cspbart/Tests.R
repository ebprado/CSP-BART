library(devtools)
library(dbarts)
remove.packages('cspbart')
load_all()
document()
check()
build()
install_github("ebprado/CSP-BART/cspbart", ref='main')
library(cspbart)
citation("cspbart")
?cspbart

# Simulate from Friedman equation -------------------------
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

# Run the semi-parametric BART (WITHOUT intercept)--------------

cspbart.fit = cspbart(formula = y ~ 0 + V4 + V5, x1 = X1, x2 = X2, ntrees = 1, nburn = 2000, npost = 1000)

# Run the semi-parametric BART (WITH intercept)--------------
# cspbart.fit = cspbart::cspbart(formula = y ~ V4 + V5, sparse = TRUE, x1 = X1, x2 = X2, ntrees = 1, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat) ------
yhat = colMeans(cspbart.fit$y_hat)
betahat = colMeans(cspbart.fit$beta_hat)

# Predict on a new dataset
yhat_pred = predict(cspbart.fit, newdata_x1 = X1, newdata_x2 = X2, type = 'mean')
cor(yhat,yhat_pred) == 1

# Plot --------------
plot(y, yhat);abline(0,1)
plot(1:2, c(10,5), main = 'True versus estimates', ylim=c(3,12))
points(1:2, betahat, col=2, pch=2)

## CSP-BART for a binary response --------------
n = 200
ncov = 5
var = 1
data = friedman_data(n, ncov, sqrt(var))

aux = data$y
y = ifelse(aux > median(aux), 1, 0)
X1 = as.data.frame(cbind(y,data$x)) # all covariates + the response
X2 = data$x # all covariates

# Run the semi-parametric BART (WITH intercept)--------------
cspbart.fit = cspbart::cl_cspbart(formula = y ~ V4 + V5, x1 = X1, x2 = X2, ntrees = 1, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat) ------
yhat = pnorm(colMeans(cspbart.fit$y_hat))
betahat = colMeans(cspbart.fit$beta_hat)