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
X1 = data.frame(y=data$y, data$x)
X1$letters = factor(rep(letters[1:25], 4))
X1$LET = rep(LETTERS[1:4], 25)
X2 = data.frame(data$x)
head(X2)
# Run the semi-parametric BART (WITHOUT intercept)--------------
# set.seed(002)
cspbart.fit = cspbart(formula = y ~ V4 + V5 + LET, x1 = X1, x2 = X2,
                      ntrees = 10, nburn = 100, npost = 10,
                      alpha = 0.99, beta = 0.01)

colMeans(cspbart.fit$beta_hat)

aa = var_used_trees(cspbart.fit, raw = FALSE)

df <- data.frame(x = factor(rep(c("a", "b", "c"), times = 3)),
                 y = as.character(rep(c("d", "e", "f"), times = 3)),
                 z = 1:9)

cspbart.fit = cspbart(formula = z ~ x + y, x1 = df, x2 = df[,-3],
                      ntrees = 10, nburn = 100, npost = 10,
                      alpha = 0.99, beta = 0.01)
colMeans(cspbart.fit$beta_hat)

## CSP-BART for a binary response --------------
n = 200
ncov = 5
var = 1
data = friedman_data(n, ncov, sqrt(var))

aux = data$y
y = ifelse(aux > median(aux), 1, 0)
X1 = as.data.frame(cbind(y,data$x)) # all covariates + the response
X2 = data$x # all covariates

# Run the semi-parametric BART (WITHOUT intercept)--------------
cspbart.fit = cspbart::cl_cspbart(formula = y ~ 0 + V4 + V5, x1 = X1, x2 = X2, ntrees = 1, nburn = 2000, npost = 1000)

# Calculate the predicted values (yhat) and parameter estimates (betahat) ------
yhat = pnorm(colMeans(cspbart.fit$y_hat))
betahat = colMeans(cspbart.fit$beta_hat)