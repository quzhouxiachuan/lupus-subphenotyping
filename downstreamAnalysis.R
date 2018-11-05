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
#delete duplicated patients; 
x$V1 = as.character(x$V1)
x= x[x$V1 != '',] # a complete script see medinfo code secured folder 
x = x[x$V1 != '',]
x = x[x$V1 != 'MRN',]
colns = colnames (read.csv('LUPUS_DATA1 1Dec2017.CSV'))
colnames(x) = colns 
#x = x[,6:16] 
#get onset age 
x$DXDATE <- as.Date(strptime(x$DXDATE, "%m/%d/%Y"))
x$BIRTHDAY <- as.Date(strptime(x$BIRTHDAY, "%m/%d/%Y")) 
x$DXAGE = x$DXDATE - x$BIRTHDAY
x$DXAGE = x$DXAGE/365
x$DXAGE = as.numeric(x$DXAGE)
#discretize lupus onset age 
x$DXAGE1 = cut(x$DXAGE, breaks = c(0,16,50,100),labels = FALSE)

x$ACRMAL = as.numeric(as.character(x$ACRMAL))
x$ACRDISC = as.numeric(as.character(x$ACRDISC)) 
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
x1 = x [,c('MRN', 'BIRTHDAY', 'SEX','RACE','DXDATE','ACRSCORE','rash','ACRPHOTO','ACRJOINT',
           'ACRSERO', 'ACRRENAL', 'ACRNEURO','ACRHEME','ACRANA','ACRIMMUN','ACRGI','DXAGE1')]


###############################################
#get patient cohort 1). having dx before 2014-03-21  2)no missing in the death data  

#enc = read.table('./death_data_oldedw.txt',sep = '\t')
enc= read.table('./death_data_oldedw.txt', head=F,
                colClasses=c('character','character','character'), sep = '\t')
enc$V1 = as.character(enc$V1)
enc[1,1] = '' # a complete analysis see medinfo code folder
x1$MRN = as.character(x1$MRN)
kk =merge(enc, x1, by.x= 'V1',by.y = 'MRN')
colnames(kk)[1:4] = c('MRN','death_flag','alive','death_date')
#delete patients without death date ; a complete script see medinfo code folder 
kk = kk[kk$MRN != '',]
kk = kk[kk$MRN != '',]
kk = kk[kk$MRN != '',]
kk = kk[kk$MRN != '',]
kk = kk[kk$MRN != ' ',]
kk = kk[kk$MRN != ' ',]
kk = kk[kk$MRN != ' ',]




## only keep individual who has DXDate before 2014-03-21
kk$DXDATE <- as.Date(kk$DXDATE)
kk$DX_timediff = as.Date('2014-03-21') - kk$DXDATE
#994 patients 
kk = kk[kk$DX_timediff>0,]
#982 patients left after filter out thos whose dx date is after 2014-03-21
kk = kk[!is.na(kk$DXDATE),]
kk = kk[kk$RACE!='', ]
colnames(kk)[19] = 'ACRULCER'

kk = kk[as.numeric(as.character(kk$ACRSCORE))>3,]
#####################################################perform LCA on the remaining cohort #########################

mydata = kk[,c(6,7,10:19)]
mydata$RACE = as.numeric(as.character(mydata$RACE)) 
mydata[mydata$RACE> 2,'RACE'] = 3 
mydata$SEX = as.numeric(as.character(mydata$SEX)) + 1 
mydata$ACRPHOTO = as.numeric(as.character(mydata$ACRPHOTO)) + 1 
mydata$ACRJOINT = as.numeric(as.character(mydata$ACRJOINT)) + 1 
mydata$ACRSERO = as.numeric(as.character(mydata$ACRSERO)) + 1 
mydata$ACRRENAL = as.numeric(as.character(mydata$ACRRENAL)) + 1 
mydata$ACRNEURO = as.numeric(as.character(mydata$ACRNEURO)) + 1 
mydata$ACRHEME = as.numeric(as.character(mydata$ACRHEME)) + 1 
mydata$ACRANA = as.numeric(as.character(mydata$ACRANA)) + 1 
mydata$ACRIMMUN = as.numeric(as.character(mydata$ACRIMMUN)) + 1 
mydata$ACRULCER = as.numeric(as.character(mydata$ACRULCER)) + 1 
mydata$rash = as.numeric(as.character(mydata$rash)) + 1 
mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)



f<-with(mydata, cbind(rash, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER)~1) 

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

#heatmap 
test = cbind(LCA_best_model$predclass,mydata)
test = test[order(test$`LCA_best_model$predclass`),]
test = sapply(test, function(x) as.numeric(as.character(x)))
matrix = data.matrix(test[,4:13])
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
big$RACE.x = as.numeric(as.character(big$RACE.x))
big[big$RACE.x>2,] = 3  

kkk$DXAGE = kkk$DXDATE - kkk$BIRTHDAY
kkk$DXAGE = kkk$DXAGE/365
kkk$SEX = factor(kkk$SEX)
kkk[as.numeric(as.character(kkk$RACE))>2, 'RACE'] = 3
kkk$RACE = factor(kkk$RACE)
cox <- coxph(Surv(ttoevent, as.numeric(death_flag)) ~ DXAGE + class + RACE+SEX, data = kkk)


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






