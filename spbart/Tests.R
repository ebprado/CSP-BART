library(devtools)
library(dbarts)
remove.packages('spbart')
load_all()
document()
check()
build()
install_github("ebprado/SP-BART/spbart", ref='main')
library(spbart)

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
x = data$x # all covariates
friedman_data = as.data.frame(cbind(y,x)) # all covariates + the response

# Run the semi-parametric BART --------------
spbart.fit = spbart::semibart(formula = y ~ 0 + V4 + V5, sparse = TRUE, x1 = friedman_data, x2 = x, ntrees = 1, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat) ------
yhat = apply(spbart.fit$y_hat,2,mean)
betahat = apply(spbart.fit$beta_hat,2,mean)

# Predict on a new dataset
yhat_pred = spbart::predict_semibart(spbart.fit, newdata_x1 = friedman_data, newdata_x2 = x, type = 'mean')

# Plot --------------
plot(y, yhat);abline(0,1)
plot(1:2, c(10,5), main = 'True versus estimates')
points(1:2, betahat, col=2, pch=2)
