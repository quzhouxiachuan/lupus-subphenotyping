library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr)
setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
x = read.csv('cld_full_slicc.csv')
x = x[,-1] 
x$dermatitis = x$slicc_alopecia + x$slicc_discoid + x$slicc_acute
x$hemo = x$slicc_leuk + x$slicc_throm + x$slicc_hemanemia
x$immu = x$slicc_dsdna +x$slicc_compliment + x$slicc_coombs + x$slicc_apa + x$slicc_sm
x$dermatitis[x$dermatitis>1] =1
x$hemo[x$hemo>1] =1
x$immu[x$immu>1] =1
x= x[,c('dermatitis','hemo','immu','slicc_ulcer','slicc_arthritis','slicc_serositis','slicc_renal','slicc_neuro','slicc_ana')]

#x = x[,3:12]
#x= sapply(x, function(x) as.factor(x))
# split x into training and testing dataset 
res = data.frame() 
for (i in 1112:1115)
{
  set.seed(i)
 # smp_size <- floor(0.9 * nrow(x))
  smp_size <- floor(1 * nrow(x))
  train_ind <- sample(seq_len(nrow(x)), size = smp_size)
  train <- x[train_ind, ]
  test <- x[-train_ind, ]
  mydata = train + 1 
  mydata = sapply(mydata, function(x) as.factor(x))
  mydata = as.data.frame(mydata)
  #f<-with(mydata, cbind(slicc_acute, slicc_discoid, slicc_ulcer,    
              #          slicc_alopecia, slicc_arthritis, slicc_serositis, slicc_renal,   
               #         slicc_neuro, slicc_hemanemia, slicc_leuk, slicc_throm, slicc_ana, slicc_dsdna, 
                #        slicc_sm, slicc_apa)~1) 

f<-with(mydata, cbind(dermatitis,hemo,immu,slicc_ulcer,slicc_arthritis,slicc_serositis,slicc_renal,slicc_neuro,slicc_ana)~1) 
  ####################################################
  # using the whole dataset for clustering 
  #mydata = x + 1 
  #mydata = sapply(mydata, function(x) as.factor(x))
  #mydata = as.data.frame(mydata)
  
  #set.seed(01125)
  #lc_whole<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)
  lc_test<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000) #,graph=TRUE)
  #lc_train<-poLCA(f, data=train+1, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)
 # new = cbind(lc_test$probs$slicc_acute[4:6],lc_test$probs$slicc_discoid[4:6],lc_test$probs$slicc_ulcer[4:6], lc_test$probs$slicc_alopecia[4:6]
 #             ,lc_test$probs$slicc_arthritis[4:6],lc_test$probs$slicc_serositis[4:6],lc_test$probs$slicc_renal[4:6]
 #             ,lc_test$probs$slicc_neuro[4:6],lc_test$probs$slicc_hemanemia[4:6],lc_test$probs$slicc_leuk[4:6],lc_test$probs$slicc_throm[4:6],
 #             lc_test$probs$slicc_ana[4:6],lc_test$probs$slicc_dsdna[4:6],lc_test$probs$slicc_sm[4:6],lc_test$probs$slicc_apa[4:6],
 #             lc_test$P)
  
new = cbind(lc_test$probs$dermatitis[4:6],lc_test$probs$hemo[4:6],lc_test$probs$immu[4:6], lc_test$probs$slicc_ulcer[4:6]
                         ,lc_test$probs$slicc_arthritis[4:6],lc_test$probs$slicc_serositis[4:6],lc_test$probs$slicc_renal[4:6]
                         ,lc_test$probs$slicc_neuro[4:6],lc_test$probs$slicc_ana[4:6],
                         lc_test$P)
  res = rbind(new, res)
  
}
colnames(res) = c('dermatitis','hemo','immu','slicc_ulcer','slicc_arthritis','slicc_serositis','slicc_renal','slicc_neuro','slicc_ana',
                  'population_share')
write.csv(res, 'LCA_lupus_bootstrapping_slicc_collapse.csv')








###additional analysis 
x = read.csv('cld_full_acr.csv')
x$age = as.Date(as.character(x$DXDATE), format="%m/%d/%Y")-as.Date(as.character(x$BIRTHDAY), format="%m/%d/%Y")
x$age = x$age/365
x$age = as.numeric(age)
kk = cbind(lc_whole$predclass, x)
kk1 = kk[kk$`lc_whole$predclass`==1,]
kk2 = kk[kk$`lc_whole$predclass`==2,]
kk3 = kk[kk$`lc_whole$predclass`==3,]
summary(kk1)
summary(kk2)
summary(kk3)







