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

x$DXDATE <- as.Date(strptime(x$DXDATE, "%m/%d/%Y"))
x$BIRTHDAY <- as.Date(strptime(x$BIRTHDAY, "%m/%d/%Y")) 
x$DXAGE = x$DXDATE - x$BIRTHDAY
x$DXAGE = x$DXAGE/365
x$AGE = as.numeric(x$DXAGE)
x = x[,c(3,6:16,18)]
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
x = x[,c(1,4:14)]

mydata = x + 1 
mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)
f<-with(mydata, cbind(rash, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER)~1)



max_II <- -100000
min_bic <- 100000
for(i in 2:4){
  lc <- poLCA(f, mydata, nclass=i, maxiter=3000, 
              tol=1e-5, na.rm=FALSE,  
              nrep=10, verbose=TRUE, calc.se=TRUE)
  if(lc$bic < min_bic){
    min_bic <- lc$bic
    LCA_best_model<-lc
  }
}       

lc_test = LCA_best_model

kk = cbind(LCA_best_model$predclass, x)
kk$`LCA_best_model$predclass` = as.factor(kk$`LCA_best_model$predclass`)
colnames(kk)[1] = 'class'
summary(aov(AGE ~ class, data = kk))

boxplot(AGE~class,data=kk, main="boxplot", 
        xlab="class", ylab="AGE")


table(kk$class,kk$SEX)

res = cbind(lc_test$probs$rash[4:6],lc_test$probs$ACRPHOTO[4:6],lc_test$probs$ACRJOINT[4:6],lc_test$probs$ACRSERO[4:6], lc_test$probs$ACRRENAL[4:6]
            ,lc_test$probs$ACRNEURO[4:6],lc_test$probs$ACRHEME[4:6],lc_test$probs$ACRANA[4:6]
            ,lc_test$probs$ACRIMMUN[4:6],lc_test$probs$ACRULCER[4:6],
            lc_test$P)

colnames(res) = c('rash','photo','ACRJOINT','ACRSERO','ACRRENAL','ACRNEURO','ACRHEME','ACRANA','ACRIMMUN','ACRULCER',
                  'population_share')

