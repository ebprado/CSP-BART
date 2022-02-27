library(devtools)
install_github("ebprado/CSP-BART/cspbart", ref='main')
library(cspbart)
library(HRW)
library(mlbench)
library(mgcv)
library(caret)
library(MKmisc)

# Calculate missclassification error -------- 
missclass_rate <- function(y_train, y_pred){
  cut_off = seq(0,1,by=0.01)
  save = NULL
  for (i in 1:length(cut_off)){
    y_pred_aux = ifelse (y_pred > cut_off[i], 1, 0)
    miss_rate = sum(y_train != y_pred_aux)/length(y_train)
    save = rbind(save, cbind(cutoff = cut_off[i], miss_rate))
  }
  rate = save[which(save[,2] == min(save[,2])),]
  if (class(rate) == 'matrix') {return(rate[1,])}
  if (class(rate) == 'numeric') {return(rate)}
}
# Pima dataset ----------
data(PimaIndiansDiabetes2)
d <- na.omit(PimaIndiansDiabetes2)
data = as.data.frame(dbarts::makeModelMatrixFromDataFrame(d))

for(i in seq_len((ncol(d) - 1))){
  plot(d[,i], d[,'diabetes'], main=colnames(d)[i])  
}

test2 <- data %>% 
  mutate(pregnant_bin = cut_number(pregnant, 5),
         glucose_bin = cut_number(glucose,10),
         pressure_bin = cut_number(pressure,10),
         triceps_bin = cut_number(triceps,10),
         insulin_bin = cut_number(insulin,10),
         mass_bin = cut_number(mass,10),
         pedigree_bin = cut_number(pedigree,10),
         age_bin = cut_number(age,10),
         )
# glucose_bin is somehow linear
# age_bin is somehow linear

tab = test2 %>%
  # group_by(age_bin) %>% 
  group_by(glucose_bin) %>% 
  summarise(n=n(),
            positive.rate=(sum(diabetes.pos)/n())*100)
  
# Which effects are approximately linear? ------------------ 
# plot(tab$age_bin, tab$positive.rate)
plot(tab$glucose_bin, tab$positive.rate)
save_results = NULL
set.seed(42)
names_x1 = c('glucose','age')
index_train = sample(1:nrow(data), ceiling(nrow(data)*0.8))
x_train = as.data.frame(data[index_train, -ncol(data)])
x_test = as.data.frame(data[-index_train, -ncol(data)])
y_train = data[index_train,'diabetes.pos']
y_test = data[-index_train,'diabetes.pos']
x1_train = x_train[,names_x1]
aux_data_train = as.data.frame(cbind(y = y_train, x_train))
aux_data_test = as.data.frame(cbind(y = y_test, x_test))

# Proposed method ---------------------------------- 
new_semi_bart = cspbart::cl_cspbart(y ~ glucose + age, x1 = aux_data_train, x2 = x_train, ntrees = 50, nburn = 5000, npost = 2000)
save(new_semi_bart,file = 'results_CSP_BART_pima.RData')
load('results_CSP_BART_pima.RData')
y_hat_new = cl_predict_cspbart(new_semi_bart,aux_data_train,x_train,'mean')
yhat_test = cl_predict_cspbart(new_semi_bart,aux_data_test,x_test,'mean')
tabA = var_used_trees(new_semi_bart, FALSE)

cspbart_miss_rate = missclass_rate(y_train, y_hat_new) # missclassification rate (training)
cspbart_miss_rate_test = missclass_rate(y_test, yhat_test) # missclassification rate (test)

# Our implementation of the HYBRID under CSP-BART priors ---------------------------------- 
aux_data_train_hybrid1 = aux_data_train
test_hybrid1 = x_test
aux_data_test_hybrid1 = aux_data_test
auxx1 = which(colnames(aux_data_train_hybrid1) %in% c('age','glucose'))
auxx2 = which(colnames(test_hybrid1) %in% c('age','glucose'))
auxx3 = which(colnames(aux_data_test_hybrid1) %in% c('age','glucose'))
colnames(aux_data_train_hybrid1)[auxx1] = c('age1', 'glucose1')
colnames(test_hybrid1)[auxx2] = c('age1', 'glucose1')
colnames(aux_data_test_hybrid1)[auxx3] = c('age1', 'glucose1')

hybrid_new_semi_bart = cspbart::cl_cspbart(y ~ glucose1 + age1, x1 = aux_data_train_hybrid1, x2 = x_train, ntrees = 50, nburn = 5000, npost = 2000)
save(hybrid_new_semi_bart,file = 'results_hybrid_under_CSP_BART_pima.RData')
load('results_hybrid_under_CSP_BART_pima.RData')
hybridy_hat_new = predict(hybrid_new_semi_bart,aux_data_train_hybrid1,x_train,'mean')
hybridyhat_test = predict(hybrid_new_semi_bart,aux_data_test_hybrid1,test_hybrid1,'mean')

hybridcspbart_miss_rate = missclass_rate(y_train, hybridy_hat_new) # missclassification rate (training)
hybridcspbart_miss_rate_test = missclass_rate(y_test, hybridyhat_test) # missclassification rate (test)

# Our implementation of the HYBRID under SSP-BART ---------------------------------- 
our_ssp_bart = cspbart::cl_sspbart(y ~ glucose + age, x1 = aux_data_train, x2 = x_train, ntrees = 50, nburn = 5000, npost = 2000)
save(our_ssp_bart,file = 'results_hybrid_under_SSP_BART_pima.RData')
load('results_hybrid_under_SSP_BART_pima.RData')
y_hat_new_our_ssp = cl_predict_cspbart(our_ssp_bart,aux_data_train,x_train,'mean')
yhat_test_our_ssp = cl_predict_cspbart(our_ssp_bart,aux_data_test,x_test,'mean')

hybridsspbart_miss_rate = missclass_rate(y_train, y_hat_new_our_ssp) # missclassification rate (training)
hybridsspbart_miss_rate_test = missclass_rate(y_test, yhat_test_our_ssp) # missclassification rate (test)

# Our implementation of the SSP-BART (WITHOUT AGE AND GLUCOSE IN X2) ---------------------------------- 
x_train22 = x_train[,-which(colnames(x_train) %in% c('age', 'glucose'))]
x_test22 = x_test[,-which(colnames(x_test) %in% c('age', 'glucose'))]

our_ssp_bart22 = cspbart::cl_sspbart(y ~ glucose + age, x1 = aux_data_train, x2 = x_train22, ntrees = 50, nburn = 5000, npost = 2000)
save(our_ssp_bart22,file = 'results_our_SPP_BART_noAGEandGLUCOSEinX2_pima.RData')
load('results_our_SPP_BART_noAGEandGLUCOSEinX2_pima.RData')
y_hat_new_our_ssp22 = predict(our_ssp_bart22,aux_data_train,x_train22,'mean')
yhat_test_our_ssp22 = predict(our_ssp_bart22,aux_data_test,x_test22,'mean')

our_sspbart_miss_rate = missclass_rate(y_train, y_hat_new_our_ssp22) # missclassification rate (training)
our_sspbart_miss_rate_test = missclass_rate(y_test, yhat_test_our_ssp22) # missclassification rate (test)

# Store misclassification rates --- 
save_results = rbind(save_results,
                     cbind(model = 'CSP-BART', set='training', rate = t(cspbart_miss_rate)),
                     cbind(model = 'CSP-BART', set='test',     rate = t(cspbart_miss_rate_test)),
                     
                     cbind(model = 'Hybrid CSP-BART', set='training', rate = t(hybridcspbart_miss_rate)),
                     cbind(model = 'Hybrid CSP-BART', set='test',     rate = t(hybridcspbart_miss_rate_test)),
                     
                     cbind(model = 'Hybrid SSP-BART', set='training', rate = t(hybridsspbart_miss_rate)),
                     cbind(model = 'Hybrid SSP-BART', set='test',     rate = t(hybridsspbart_miss_rate_test)),
                     
                     cbind(model = 'Our SSP-BART', set='training', rate = t(our_sspbart_miss_rate)),
                     cbind(model = 'Our SSP-BART', set='test',     rate = t(our_sspbart_miss_rate_test))
)

save(save_results, file='00_misclassification_rates.RData')

# Parameter estimates --------------------------- 

# CSP-BART
beta_hat_new = colMeans(new_semi_bart$beta_hat)
q5_beta_new = apply(new_semi_bart$beta_hat, 2, quantile, prob=0.05)
q95_beta_new = apply(new_semi_bart$beta_hat, 2, quantile, prob=0.95)

# SSP-BART
beta_hat_new_ssp22 = colMeans(our_ssp_bart22$beta_hat)
q5_beta_new_ssp22 = apply(our_ssp_bart22$beta_hat, 2, quantile, prob=0.05)
q95_beta_new_ssp22 = apply(our_ssp_bart22$beta_hat, 2, quantile, prob=0.95)
