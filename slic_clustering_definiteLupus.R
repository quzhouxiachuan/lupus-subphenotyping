setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/YU/medoinfo Conference')
library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr)
x = read.table('LUPUS_DATA1 1Dec2017.CSV', sep = ',')
#delete duplicated patients; 08506490 and 07923366
x$V1 = as.character(x$V1)
x= x[x$V1 != '08506490',]
x = x[x$V1 != '07923366',]
x = x[x$V1 != 'MRN',]
colns = colnames (read.csv('LUPUS_DATA1 1Dec2017.CSV'))
colnames(x) = colns 
x$ACRSCORE = as.numeric(as.character(x$ACRSCORE))
x = x[x$ACRSCORE> 3,]
slic = read.table('../cld_slicc_final.csv',sep = ',')
col_slic = colnames (read.csv('../cld_slicc_final.csv'))
colnames(slic) = col_slic 
slic = slic[slic$MRN %in% x$MRN,] #839 patients left after deleting the duplicated patients 
#and the probable lupus patients 
kk = merge(slic, x, by = 'MRN')


###############################################
#get patient cohort 1). having dx before 2014-03-21  2)no missing in the death data  

#enc = read.table('./death_data_oldedw.txt',sep = '\t')
enc= read.table('./death_data_oldedw.txt', head=F,
                colClasses=c('character','character','character'), sep = '\t')
enc$V1 = as.character(enc$V1)
enc[1,1] = '00001347'
kk$MRN = as.character(kk$MRN)
kk =merge(enc, kk, by.x= 'V1',by.y = 'MRN')
colnames(kk)[1:4] = c('MRN','death_flag','alive','death_date')
#00317490 07439482 07679455  07727243 08260731 13524359 00311497 do not have death date 
kk = kk[kk$MRN != '00311497',]
kk = kk[kk$MRN != '00317490',]
kk = kk[kk$MRN != '07439482',]
kk = kk[kk$MRN != '07679455',]
kk = kk[kk$MRN != '07727243',]
kk = kk[kk$MRN != '08260731',]
kk = kk[kk$MRN != '13524359',]


## only keep individual who has DXDate before 2014-03-21
kk$DXDATE <- anydate(kk$DXDATE)
kk$DX_timediff = as.Date('2014-03-21') - kk$DXDATE
#994 patients 
kk = kk[kk$DX_timediff>0,]
#982 patients left after filter out thos whose dx date is after 2014-03-21
kk = kk[!is.na(kk$DXDATE),]
kk = kk[kk$RACE!='', ]



#####################################################perform LCA on the remaining cohort #########################

mydata = kk 
mydata$slicc_acute = as.numeric(as.character(mydata$slicc_acute)) + 1 
mydata$slicc_discoid = as.numeric(as.character(mydata$slicc_discoid)) + 1 
mydata$slicc_ulcer = as.numeric(as.character(mydata$slicc_ulcer)) + 1 
mydata$slicc_alopecia = as.numeric(as.character(mydata$slicc_alopecia)) + 1 
mydata$slicc_arthritis = as.numeric(as.character(mydata$slicc_arthritis)) + 1 
mydata$slicc_serositis = as.numeric(as.character(mydata$slicc_serositis)) + 1 
mydata$slicc_renal = as.numeric(as.character(mydata$slicc_renal)) + 1 
mydata$slicc_neuro = as.numeric(as.character(mydata$slicc_neuro)) + 1 
mydata$slicc_hemanemia = as.numeric(as.character(mydata$slicc_hemanemia)) + 1 
mydata$slicc_leuk = as.numeric(as.character(mydata$slicc_leuk)) + 1 
mydata$slicc_throm = as.numeric(as.character(mydata$slicc_throm)) + 1 
mydata$slicc_ana = as.numeric(as.character(mydata$slicc_ana)) + 1 
mydata$slicc_dsdna = as.numeric(as.character(mydata$slicc_dsdna)) + 1 
mydata$slicc_sm = as.numeric(as.character(mydata$slicc_sm)) + 1 
mydata$slicc_apa = as.numeric(as.character(mydata$slicc_apa)) + 1 

mydata$slicc_compliment = as.numeric(as.character(mydata$slicc_compliment)) + 1 

mydata$slicc_coombs = as.numeric(as.character(mydata$slicc_coombs)) + 1 

mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)



f<-with(mydata, cbind(slicc_acute,slicc_discoid,slicc_ulcer,slicc_alopecia, slicc_arthritis, slicc_serositis, slicc_renal, slicc_neuro, slicc_hemanemia,
                      slicc_leuk, slicc_throm, slicc_ana, slicc_dsdna, slicc_sm, slicc_apa, slicc_compliment, slicc_coombs)~1) 


max_II <- -100000
min_bic <- 100000
for(i in 2:10){
  lc <- poLCA(f, mydata, nclass=i, maxiter=3000, 
              tol=1e-5, na.rm=FALSE,  
              nrep=10, verbose=TRUE, calc.se=TRUE)
  if(lc$bic < min_bic){
    min_bic <- lc$bic
    LCA_best_model<-lc
  }
}       
LCA_best_model

kkk = cbind(LCA_best_model$predclass, kk)
kkk$`LCA_best_model$predclass` = as.factor(kkk$`LCA_best_model$predclass`)
colnames(kkk)[1] = 'class'

kkk$ttoevent = 'empty'
kkk$death_flag = as.character(kkk$death_flag)
for (i in (1:dim(kkk)[1]))
{
  if (kkk[i,'death_flag'] == 'NULL'){
    kkk[i,'ttoevent'] = as.Date('2014-03-21') - as.Date(kkk$DXDATE[i])
  }
  else
  {
    kkk[i,'ttoevent'] = as.Date(kkk$death_date[i]) -  as.Date(kkk$DXDATE[i])
    
  }
  
}
kkk$ttoevent = as.numeric(kkk$ttoevent)
kkk$death_flag[kkk$death_flag=='NULL'] = 0 





library(survival)
library(survminer)
#colnames(big)[2] = 'class'

fit <- survfit(Surv(ttoevent, as.numeric(death_flag)) ~ class, data = kkk)
ggsurvplot(fit, data = kkk, risk.table = TRUE, pval = TRUE, conf.int = TRUE )



test = cbind(LCA_best_model$predclass,mydata)
test = test[order(test$`LCA_best_model$predclass`),]
test = test[,6:22]
test = sapply(test, function(x) as.numeric(as.character(x)))
matrix = data.matrix(test)
#matrix = as.matrix(test[,4:13])
heatmap(as.matrix(matrix), scale = "none", Rowv = NA, Colv = NA, col = cm.colors(2), main = "HeatMap Example")




#downstream analysis 
#Q1: is class1 have older onset patients ?
big = merge(kkk, x, by='MRN')
#all the patients, onset age, min(dxage) = 4.9 max(dxage) = 80, mean = 30.38
#class1: min = 11.4; max = 80; mean = 33.34 ; median = 30 
#class2: min = 9; max = 46; mean = 23.7; median = 22 
#class 3: min = 4.9 ; max = 68.255; mean = 31; median = 30 
big$age = as.Date('2018-10-29') - as.Date(big$BIRTHDAY.x)
big$age = as.numeric(big$age)/365
summary(big[big$class==1,])
#class1: min = 22, max = 107, mean = 53, median = 52.93
#class2: min = 22, max = 82, mean = 48.59, median = 48.5
#class3: minn = 21.59, max = 99.63, mean = 54.70, median = 54.79
cox <- coxph(Surv(ttoevent, as.factor(death_flag)) ~ DXAGE + class, data = big)


#number of probable lupus 
kkk$ACRSCORE = as.numeric(kkk$ACRSCORE)
for (i in (1:dim(kkk)[1]))
{
  if (kkk[i,'ACRSCORE'] > 3){
    kkk[i,'definite_SLE'] = 1
  }
  else
  {
    kkk[i,'definite_SLE'] = 0
    
  }
  
}






