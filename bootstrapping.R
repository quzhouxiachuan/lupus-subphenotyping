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
x = read.csv('cld_full_acr.csv')
x = x[,6:16] 
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
x = x[,3:12]
#x= sapply(x, function(x) as.factor(x))
# split x into training and testing dataset 
res = data.frame() 
for (i in 1112:1120)
{
  set.seed(i)
smp_size <- floor(0.9 * nrow(x))
train_ind <- sample(seq_len(nrow(x)), size = smp_size)
train <- x[train_ind, ]
test <- x[-train_ind, ]
mydata = train + 1 
mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)
f<-with(mydata, cbind(rash, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER)~1) 
####################################################
# using the whole dataset for clustering 
#mydata = x + 1 
#mydata = sapply(mydata, function(x) as.factor(x))
#mydata = as.data.frame(mydata)

#set.seed(01125)
#lc_whole<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)
lc_test<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)
#lc_train<-poLCA(f, data=train+1, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)
new = cbind(lc_test$probs$rash[4:6],lc_test$probs$ACRPHOTO[4:6],lc_test$probs$ACRJOINT[4:6], lc_test$probs$ACRSERO[4:6]
            ,lc_test$probs$ACRRENAL[4:6],lc_test$probs$ACRNEURO[4:6],lc_test$probs$ACRHEME[4:6]
            ,lc_test$probs$ACRANA[4:6],lc_test$probs$ACRIMMUN[4:6],lc_test$probs$ACRULCER[4:6],lc_test$P)

res = rbind(new, res)
write.csv(res, 'LCA_lupus_bootstrapping.csv')

}








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







