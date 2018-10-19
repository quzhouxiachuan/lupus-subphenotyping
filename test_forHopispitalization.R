library(dplyr)
library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr)

setwd('R:/IPHAM/CHIP/Data_Team/Projects/Lupus/YU/')
x = read.table('cld_full_acr.csv', sep = ',' )
y = read.table('inpatient_enc_noreasonforvisit.csv', sep = ',') #, fileEncoding = 'UTF-8-BOM')
mm = merge(x, y, by= 'V1') #number of patients? # this will delete patients without an inpatient enc 


colnames(mm) = c('MRN', 'BIRTHDAY', 'sex','race', 'DXdate', 'acrscore','acrmal','acrdisc'
                 , 'acrphoto','acrjoint','acrsero','acrrenal', 'acrneuro','acrheme', 'acrana'
                 ,'acrimmun','acrulcer', 'curpt_irid','enc_ir_id','admission_date','discharge_date'
                 ,'enc_startdate','enc_enddate','primary_dx')
mm$DXdate <- as.Date(strptime(mm$DXdate, "%m/%d/%Y"))
mm$admission_date <- as.Date(mm$admission_date)
mm$ttoevent = (mm$admission_date - mm$DXdate) #how many months difference 
mm = mm[mm$ttoevent>4 ,] #what about those with the NULL ??? 354 patients left 
#patient 07163561 has inpatient enc on 2007-12-18 for the reason of lupus 
df = mm %>%
  group_by(MRN) %>%
  mutate(my_ranks = order(order(ttoevent, decreasing=F)))
#kk = mm 

df = df[df$my_ranks==1,]
df = as.data.frame(df)
df = na.omit(df)
x = df 
x$acrmal = as.numeric(as.character(x$acrmal))
x$acrdisc = as.numeric(as.character(x$acrdisc))
x$rash = x$acrmal + x$acrdisc
x$rash[x$rash>1] =1 
x = x[,c('acrphoto', 'acrjoint','acrsero', 'acrrenal', 'acrneuro', 'acrheme', 'acrana', 'acrimmun', 'acrulcer','rash')]

x$acrphoto =as.numeric(as.character( x$acrphoto) )+ 1
x$acrjoint =as.numeric(as.character( x$acrjoint) )+ 1
x$acrsero =as.numeric(as.character( x$acrsero) )+ 1
x$acrrenal =as.numeric(as.character( x$acrrenal) )+ 1
x$acrheme =as.numeric(as.character( x$acrheme) )+ 1
x$acrana =as.numeric(as.character( x$acrana) )+ 1
x$acrimmun =as.numeric(as.character( x$acrimmun) )+ 1
x$acrulcer =as.numeric(as.character( x$acrulcer) )+ 1
x$rash =as.numeric(as.character( x$rash) )+ 1
mydata = x 

mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)

f<-with(mydata, cbind(rash, acrphoto, acrjoint,    
                      acrsero, acrrenal, acrneuro, acrheme,   
                      acrana, acrimmun, acrulcer)~1) 


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


set.seed(01125)
lc3<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000)
LCA_best_model = lc3 
kk = cbind(LCA_best_model$predclass, x)
kk$`LCA_best_model$predclass` = as.factor(kk$`LCA_best_model$predclass`)

fit <- survfit(Surv(time, enc) ~ class, data = kk)
ggsurvplot(fit, data = kk, risk.table = TRUE, pval = TRUE, conf.int = TRUE )







