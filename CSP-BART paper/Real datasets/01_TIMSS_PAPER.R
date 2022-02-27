library(haven)
library(tidyverse)
library(dbarts)
library(reshape2)
library(dplyr)
library(devtools)
install_github("ebprado/CSP-BART/cspbart", ref = 'main')
# install.packages('bcf')
library(cspbart)
library(bcf)
library(VCBART)
library(xtable)
source('00_aux.R')

RMAE = function(true, predicted){
  sqrt(mean(abs(true - predicted)))
}
RMedAE = function(true, predicted){
  sqrt(median(abs(true - predicted)))
}
# The datasets below can be downloaded in the section "Download TIMSS 2019 Data Files by Country" at https://timss2019.org/international-database/.
# To reproduce our results, tick "Ireland" in the "Grade 8" and the dowload the files (SPSS Data).


# BSBG05B B eighth grade, S student, B rasch score, G general question, 
# BSMMAT01 B eighth grade, S student, M Mathematics, 

# Page 78 explains the variable naming convention
# https://timss2019.org/international-database/downloads/TIMSS-2019-User-Guide-for-the-International-Database.pdf

bcg = read_sav('bcgirlm7.sav') # Grade 8 school context data files ()
# bsa = read_sav('bsairlm7.sav') # Grade 8 student achievement data files
bsg = read_sav('bsgirlm7.sav') # Grade 8 student context data files
# bsr = read_sav('bsrirlm7.sav') # Grade 8 within-country scoring reliability data files
bst = read_sav('bstirlm7.sav') # Grade 8 student-teacher linkage files
btm = read_sav('btmirlm7.sav') # Grade 8 mathematics teacher context data files
# bts = read_sav('btsirlm7.sav') # Grade 8 science teacher context data files

bcg = bcg[,keep_bcg]
bsg = bsg[,keep_bsg]
bst = bst[,keep_bst]
btm = btm[,keep_btm]

# Merge data sets ------------------ 

aa = left_join(bst, btm, by='IDTEALIN') # a lot of missing from btm

bb = left_join(bsg, bcg, by='IDSCHOOL')

cc = left_join(aa, bb, by='IDSTUD')

cc = cc[,keep_20ish]

# Find repeated columns
aux = as.data.frame(names(cc))
aux$rep = grepl('\\.',aux$`names(cc)`)
non_rep_vars = aux %>% filter(rep == FALSE) %>% dplyr::select('names(cc)')
non_rep_vars = non_rep_vars[,1]

# Remove repeated columns
dd = cc[,non_rep_vars]
aux_na = as.data.frame(sapply(dd, function(x) sum(is.na(x))))
colnames(aux_na) = 'count_na'
aux_na = aux_na/nrow(dd)
aux_na$cov = rownames(aux_na)
ee = aux_na %>% filter(count_na < 0.1) %>% dplyr::select(cov)

set.seed(001)
ff = na.omit(dd[, ee[,1]])
# ff$ITSEX = as.factor(ff$ITSEX)
# ff$BSDGEDUP = as.factor(ff$BSDGEDUP)

s <- sample(nrow(ff), round(.7*nrow(ff)))
y = ff$BSMMAT01[s]
y_test = ff$BSMMAT01[-s]
x = as.data.frame(ff[s,-which(names(ff) %in% c('BSMMAT01', "IDTEALIN", "IDSTUD", "IDBOOK", "IDCLASS"))])
x_test = as.data.frame(ff[-s,-which(names(ff) %in% c('BSMMAT01', "IDTEALIN", "IDSTUD", "IDBOOK", "IDCLASS"))])

# -------------------------------------------------------------------------
# Causal type analysis 
# -------------------------------------------------------------------------

data = as.data.frame(cbind(x, BSMMAT01 = ff$BSMMAT01[s]))
data$BCDGDAS = ifelse(data$BCDGDAS %in% c(1,2), 0, 1) # 1 = severe discipline issues
data$BSMMAT01 = as.numeric(data$BSMMAT01)

x2 = data[,-which(colnames(data) == 'BSMMAT01')]

# Proposed method ---------------------------------- 
cspbart = cspbart::cspbart(BSMMAT01 ~ 1 + BCDGDAS, x1 = data, x2, ntrees = 50, nburn = 10000, npost=2000)
apply(cspbart$beta_hat,2,quantile, probs=0.05)
apply(cspbart$beta_hat,2,quantile, probs=0.95)
apply(cspbart$beta_hat,2,quantile, probs=0.5)
colMeans(cspbart$beta_hat)
hist(cspbart$beta_hat[,2])
save(cspbart,file = 'results_causal_double_BART_TIMSS2019.RData')
load('results_causal_double_BART_TIMSS2019.RData')
auxtrees = var_used_trees(cspbart)

# BCF ---------------------------------- 
estimate_z = dbarts::bart2(BCDGDAS ~ ., data=data)
# pihat = apply(pnorm(estimate_z$yhat.train),3,mean)
pihat = rep(0.5, nrow(data))
data2 = as.matrix(data[,-which(colnames(data) %in% c('BSMMAT01', 'BCDGDAS'))])
# data2 = as.matrix(data[,-which(colnames(data) %in% c('BSMMAT01'))])
bcf = bcf(y = data$BSMMAT01, z = data$BCDGDAS, x_control = data2, x_moderate = data2,  pihat = pihat, nburn = 5000, nsim = 2000)
save(bcf,file = 'results_bcf_TIMSS2019.RData')
load('results_bcf_TIMSS2019.RData')

# mean(colMeans(bcf$tau)[data$BCDGDAS==1])
# mean(colMeans(bcf$tau)[data$BCDGDAS==0])
mean(colMeans(bcf$tau))
mean(apply(bcf$tau,2,quantile, probs=c(0.025)))
mean(apply(bcf$tau,2,quantile, probs=c(0.5)))
mean(apply(bcf$tau,2,quantile, probs=c(0.975)))
hist(colMeans(bcf$tau))

### plot -----------------------------------------

a1 = data.frame(id = 'BCF', estimate = apply(bcf$tau,2,mean))
a2 = data.frame(id = 'cspbart', estimate = cspbart$beta_hat[,2])
a3 = rbind(a1,a2)

a3 %>% ggplot(aes(x=id, y=estimate, colour=id)) +
  geom_boxplot() +
  labs(title='Treatment effect estimates',
       x = 'Methods',
       y = 'Parameter estimates') + 
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ------------------------------------------------------------ 
# VCBART and SP-BART
# ------------------------------------------------------------ 
# VCBART
seed = 001
set.seed(seed)
data = as.data.frame(cbind(x, BSMMAT01 = ff$BSMMAT01[s]))
#data$BCDGDAS = ifelse(data$BCDGDAS %in% c(1,2), 0, 1) # 1 = severe discipline issues
data$BSMMAT01 = as.numeric(data$BSMMAT01)
s_train = sample(nrow(data), round(0.8*nrow(data)))
Y_train = data[s_train,'BSMMAT01']
Z_train = as.matrix(data[s_train, which(colnames(data) != 'BSMMAT01')])
X_train = data[s_train,c('BSDGEDUP', 'BSBM42BA', 'BCDGDAS')]
X_train$BSDGEDUP = as.factor(X_train$BSDGEDUP)
X_train$BSBM42BA = as.factor(X_train$BSBM42BA)
X_train$BCDGDAS = as.factor(X_train$BCDGDAS)
X_train = dbarts::makeModelMatrixFromDataFrame(X_train)
n_train = nrow(Z_train)
cutpoints <- lapply(as.data.frame(Z_train), unique)

Y_test =  data[-s_train,'BSMMAT01']
X_test = data[-s_train,c('BSDGEDUP', 'BSBM42BA', 'BCDGDAS')]
X_test$BSDGEDUP = as.factor(X_test$BSDGEDUP)
X_test$BSBM42BA = as.factor(X_test$BSBM42BA)
X_test$BCDGDAS = as.factor(X_test$BCDGDAS)
X_test = dbarts::makeModelMatrixFromDataFrame(X_test)
Z_test = as.matrix(data[-s_train, which(colnames(data) != 'BSMMAT01')])
n_test = nrow(Z_test)

chain1 <- VCBART(Y_train, X_train, Z_train, n_train, cutpoints = cutpoints,
                 X_test, Z_test, n_test,
                 intercept = TRUE, M = 50, error_structure = "ind",
                 split_probs_type = "adaptive", ht_sigma_y = TRUE,
                 nd = 2000, burn = 5000, verbose = TRUE, print_every = 50)
save(chain1,file = 'results_VCBART_no_vary_by_sex_TIMSS2019.RData')
load('results_VCBART_no_vary_by_sex_TIMSS2019.RData')

beta_summary <- summarize_beta(chain1, chain1, burn = 5000)
ystar_summary <- summarize_posterior_predictive(chain1, chain1, burn = 5000)
RMSE(Y_test, ystar_summary$test[,'MEAN'])
apply(beta_summary$train[,'MEAN',],2,mean) # intercept

# CSP-BART ----------------------

dat1 = as.data.frame(cbind(x, BSMMAT01 = ff$BSMMAT01[s]))
dat1$BSBG10 = as.factor(dat1$BSBG10)
dat1$BSBG11B = as.factor(dat1$BSBG11B)
dat1$BCDGDAS = as.factor(dat1$BCDGDAS)
dat1$BSDGEDUP = as.factor(dat1$BSDGEDUP)
dat1$BSBM42BA = as.factor(dat1$BSBM42BA)
dat1$ITSEX = as.factor(dat1$ITSEX)
dat1$BSMMAT01 = as.numeric(dat1$BSMMAT01)
dat1_train = dat1[s_train,]
dat1_test = dat1[-s_train,]

datx2_train = dat1[s_train, -which(colnames(dat1) == 'BSMMAT01')]
datx2_test = dat1[-s_train, -which(colnames(dat1) == 'BSMMAT01')]

cspbart = cspbart::cspbart(BSMMAT01 ~ 1 + BSDGEDUP + BSBM42BA + BCDGDAS, x1 = dat1_train,
                            datx2_train, ntrees = 50, nburn = 5000, npost=2000)
head(cspbart$beta_hat)
save(cspbart,file = 'results_CSP_BART_TIMSS2019.RData')
load('results_CSP_BART_TIMSS2019.RData')
y_hat_test = predict_cspbart(cspbart, dat1_test, datx2_test, type = 'mean')
colMeans(cspbart$beta_hat)

auxtrees = var_used_trees(cspbart, raw = FALSE)
auxtrees = auxtrees %>%
  mutate(prop = round((Freq/sum(Freq))*100,4)) %>% 
  arrange(desc(prop)) %>% 
  mutate(propcum = cumsum(prop))
auxtrees$absenteeism = ifelse(grepl('BSBG10', auxtrees$vars_trees, fixed=TRUE) == TRUE, 1, 0)
auxtrees$disciplineissues = ifelse(grepl('BCDGDAS', auxtrees$vars_trees, fixed=TRUE) == TRUE, 1, 0)
auxtrees$parentseducation = ifelse(grepl('BSDGEDUP', auxtrees$vars_trees, fixed=TRUE) == TRUE, 1, 0)
auxtrees$homework = ifelse(grepl('BSBM42BA', auxtrees$vars_trees, fixed=TRUE) == TRUE, 1, 0)
auxtrees$hungry = ifelse(grepl('BSBG11B', auxtrees$vars_trees, fixed=TRUE) == TRUE, 1, 0)