library(devtools)
install_github("ebprado/SP-BART/spbart", ref = 'main')
install_github("ebprado/semibart")

library(spbart)
library(semibart)
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
  return(rate)
}

# Pima dataset ----------
data(PimaIndiansDiabetes2)
d <- na.omit(PimaIndiansDiabetes2)
data = as.data.frame(dbarts::makeModelMatrixFromDataFrame(d))

for(i in 1:(ncol(d) - 1)){
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
new_semi_bart = spbart::cl_semibart(y ~ glucose + age, x1 = aux_data_train, x2 = x_train, ntrees = 50, nburn = 5000, npost = 2000)
save(new_semi_bart,file = 'results_double_BART_pima.RData')
load('results_double_BART_pima.RData')
y_hat_new = cl_predict_semibart(new_semi_bart,aux_data_train,x_train,'mean')
yhat_test = cl_predict_semibart(new_semi_bart,aux_data_test,x_test,'mean')
beta_hat_new = apply(new_semi_bart$beta_hat, 2, mean)
cor(y_train, y_hat_new)
cor(y_test, yhat_test)
tabA = var_used_trees(new_semi_bart, FALSE)

spbart_miss_rate = missclass_rate(y_train, y_hat_new) # missclassification rate (training)
spbart_miss_rate_test = missclass_rate(y_test, yhat_test) # missclassification rate (test)

# Semi-parametric BART of Zeldow et al (2019) (X1 \cap X2 = empty) ------------------- 
L1 = as.matrix(x_train[,!(colnames(x_train) %in% names_x1)]) # nuisance covariates
L2 = as.matrix(cbind(intercept = 1, x1_train)) # variables of interest
cur_semi_bart = semibart::semibart(x.train = L1, a.train = L2, y.train = y_train, binarylink='probit', ntree = 50, ndpost = 7000)
save(cur_semi_bart,file = 'results_semi_BART_pima_intersection_empty.RData')
load('results_semi_BART_pima_intersection_empty.RData')
beta_hat = apply(cur_semi_bart$beta[-c(1:5000),], 2, mean) # remove burn-in
y_hat_sbart1 = pnorm(L2%*%beta_hat + apply(cur_semi_bart$bartfit, 2, mean))
cor(y_train, y_hat_sbart1)
cur_semi_bart1_miss_rate = missclass_rate(y_train, y_hat_sbart1) # missclassification rate (training)

# Semi-parametric BART of Zeldow et al (2019) (X1 \cap X2 not equal to empty) ---------------------
cur_semi_bart2 = semibart::semibart(x.train = as.matrix(x_train), a.train = L2, y.train = y_train, binarylink='probit',ntree = 50, ndpost = 7000)
save(cur_semi_bart2,file = 'results_semi_BART_pima_intersection_NO_empty.RData')
load('results_semi_BART_pima_intersection_NO_empty.RData')
beta_hat2 = apply(cur_semi_bart2$beta[-c(1:5000),], 2, mean) # remove burn-in
y_hat_sbart2 = pnorm(L2%*%beta_hat2 + apply(cur_semi_bart2$bartfit, 2, mean))
cor(y_train, y_hat_sbart2)
cur_semi_bart2_miss_rate = missclass_rate(y_train, y_hat_sbart2) # missclassification rate (training)

# Credible intervals -------------------

q5_beta_new = apply(new_semi_bart$beta_hat, 2, quantile, prob=0.05)
q95_beta_new = apply(new_semi_bart$beta_hat, 2, quantile, prob=0.95)


q5_beta_old = apply(cur_semi_bart$beta, 2, quantile, prob=0.05)
q95_beta_old = apply(cur_semi_bart$beta, 2, quantile, prob=0.95)


new_sbart = cbind(estimate = beta_hat_new, lower = q5_beta_new, upper = q95_beta_new)
old_sbart = cbind(estimate = beta_hat, lower = q5_beta_old, upper = q95_beta_old)
rbind(new_sbart,old_sbart)

# Misclassification error ----------------------

spbart_miss_rate

cur_semi_bart1_miss_rate

cur_semi_bart2_miss_rate

bart_miss_rate

# Most used variables ---------------------------

saveTrees = bart$fit$getTrees() 
names_x <- colnames(x_train)
auxTrees = saveTrees %>% filter(var!=-1)

# auxTrees2 = dcast(data = auxTrees, chain + sample + tree ~ var, fun.aggregate = function(x) sum(!is.na(x)), value.var='value')
auxTrees2 = dcast(data = auxTrees, chain + sample + tree ~ var, fun.aggregate = function(x) sum(!is.na(x)), value.var='value')
namesauxTree2 = colnames(auxTrees2)[-c(1,2,3)]
test = auxTrees2[, namesauxTree2]
test$var_used_tree2 = apply(test, 1, function(x) paste(sort(names_x[as.numeric(namesauxTree2[x > 0])]), collapse = ','))

tab = test %>% 
  group_by(var_used_tree2) %>% 
  summarise(n=n()) %>% 
  mutate(pc=round(n/sum(n),3)) %>% 
  arrange(desc(pc)) %>% 
  head(20)

tab$var_used_tree2 = factor(tab$var_used_tree2, levels = tab$var_used_tree2, labels=tab$var_used_tree2)

ggplot(data=tab, aes(x=var_used_tree2, y=pc)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_bw(base_size = 12) +
  labs(title='Most used variables by the trees',
       x = 'Variables',
       y = 'Proportion of times') +
  theme(text = element_text(size=20),
        plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(aes(label=pc), vjust=1.6, color="white", size=5)
