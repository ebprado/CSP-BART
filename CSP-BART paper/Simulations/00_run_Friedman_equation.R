library(devtools)
install_github("ebprado/SP-BART/spbart", ref = 'main')
install_github("ebprado/semibart") # a modified version of the implementation of the semi-parametric BART (Zeldow et al, 2019)

library(spbart)
library(semibart)
library(mgcv)
library(VCBART)

list_n = c(1000) # number of observations
list_var = c(1, 10) # error variance
list_ncov = c(10, 50) # number of covariates
all.comb = expand.grid(n=list_n, var=list_var, ncov = list_ncov) # create all scenarios given n, var, and ncov
seed = 0
save.path = "~/R/semiBART/results/"

for (i in 1:10){ # Number of Monte Carlo repetitions
  for (j in 1:nrow(all.comb)){
    
    seed = seed + 1
    n = all.comb[j,'n']
    var = all.comb[j,'var']
    ncov = all.comb[j,'ncov']
    
    set.seed(seed)
    
    # Simulate from Friedman equation -------------------------
    friedman_data = function(n, num_cov, sd_error){
      x = matrix(runif(n*num_cov),n,num_cov)
      y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
      return(list(y = y,
                  x = as.data.frame(x)))
    }
    
    data = friedman_data(n, ncov, sqrt(var))
    y = data$y # response variable
    x = data$x # all covariates
    aux_data = as.data.frame(cbind(y,x))
    
    x2 = x[,-c(4,5)]
    
    xx = as.matrix(x)
    xx1 = xx[,4:5] # covariates of primary interest (x4 and x5)
    xx2 = xx[,-c(4:5)]
    colnames(xx1) = c('x4','x5')
    cutpoints <- lapply(as.data.frame(xx), unique)
    n_train = length(y)
    
    # Run models -------------------------------------------
    
    # Semi-parametric BART of Zeldow et al (2019)
    old_semibart2 = semibart::semibart(x.train = xx2, a.train = xx1, y.train = y, ndpost = 4000, ntree = 50)
    
    # Proposed method
    new_semibart = spbart::semibart(formula = y ~ 0 + V4 + V5, sparse = TRUE, aux_data, x2 = x, ntrees = 50, nburn = 2000, npost = 1000)
    
    # Generalised Additive Models
    gam_mod <- gam(y ~ x[,4] + x[,5] + s(x[,1], x[,2]) + s(x[,3]))
    
    # Varying Coefficient BART
    vcbart <- VCBART::VCBART(y, X_train = xx1, Z_train = xx, X_test = xx1, Z_test = xx, n_train = n_train, n_test = n_train, cutpoints = cutpoints, intercept = FALSE, burn = 2000, nd = 1000, error_structure = "ind", split_probs_type = 'fixed')
    
    id = paste('Fried','n', n, 'sd', var, 'ncov', ncov, 'rep', i,  sep='')
    
    file_name_old2   = paste(save.path, id, 'old_semi_bart2', '.RData', sep='')
    file_name_new    = paste(save.path, id, 'new_semi_bart',  '.RData', sep='')
    file_name_vcbart = paste(save.path, id, 'vc_bart',        '.RData', sep='')
    file_name_gam    = paste(save.path, id, 'gam',            '.RData', sep='')
    file_name_data   = paste(save.path, id, 'data',           '.RData', sep='')
    
    # Save results -------------------------------------------
  
    save(old_semibart2, file = file_name_old2)
    save(new_semibart,  file = file_name_new)
    save(vcbart,        file = file_name_vcbart)
    save(gam_mod,       file = file_name_gam)
    save(data,          file = file_name_data)
  } 
}
