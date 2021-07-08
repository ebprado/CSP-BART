library(devtools)
library(dbarts)
remove.packages('semiBART')
load_all()
document()
check()
build()
install_github("ebprado/ExtensionsBART/semiBART")

apply(ambarti$g_hat, 2, mean) - ambarti$y_mean






# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  # y = 20 + 10*x[,4]-5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = as.data.frame(x)))
}
# Training data
data = friedman_data(200, 10, 1)
y = data$y
x = data$x

lm(y ~ V4 +V5, data=x)
lm(y ~ -1 + V4 +V5, data=x)

x1 = data.frame(y = data$y, data$x)
formula = y ~ V1 + V2
x2 = x1[,-1]
oi = semibart(formula, x1, x2, ntrees = 50, nburn = 200, npost = 100)
tab1 = var_used_trees(oi)
