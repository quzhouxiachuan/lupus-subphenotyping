library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr)
#setwd('R:/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/YU/')
x = read.table('cld_full_acr.csv', sep = ',')
y = read.table('./cld_slicc_yu.csv', sep = ',',)
mm = merge(x, y, by = 'V1')
#delete duplicated patients; 08506490 and 07923366
mm$V1 = as.character(mm$V1)
mm = mm[mm$V1 != '08506490',]
mm = mm[mm$V1 != '07923366',]
mm = mm[mm$V1 != 'MRN',]

colnames(mm) = c('MRN', 'BIRTHDAY', 'SEX','RACE','DXDATE','ACRSCORE','ACRMAL','ACRDISC','ACRPHOTO','ACRJOINT',
                 'ACRSERO', 'ACRRENAL', 'ACRNEURO','ACRHEME','ACRANA','ACRIMMUN','ACRULCER',
                 'slicc_acute', 'slicc_discoid', 'slicc_ulcer', 'slicc_alopecia', 
                 'slicc_arthritis', 'slicc_serositis', 'slicc_renal', 'slicc_neuro',
                 'slicc_hemanemia', 'slicc_leuk', 'slicc_throm', 'slicc_ana', 'slicc_dsdna'
                 ,'slicc_sm', 'slicc_apa', 'slicc_compliment', 'slicc_coombs',
                 'BIRTHDAY2', 'SEX2', 'RACE2', 'MARITAL', 'EDLEVEL', 'EMPLOY', 'DXDATE2')

#x = x[,6:16] 
#get onset age 
x = mm
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
x = x [,c(1:17,42, 43,44)]
x = x[,-c(2,6,7,8)]

inpt = read.table('inpatient_enc_noreasonforvisit.csv', sep = ',')
inpt$inpt_ind = 1
outpt = read.table('lupus_outpatient_enc.csv', sep = ',') 
outpt$outpt_ind = 0
#inout = merge(inpt, outpt, by.x = 'V1',by.y = 'V1', all = T)
#colnames(inout) = c('nmff_mrn', 'ir_id', 'enc_ir_id', 'admission_date_key' ,'discharge_date_key' 
                    ,'encounter_start_date_key', 'encounter_end_date_key', 'primary_diagnosis_key'
                    , 'inpt_ind',  'out_ir_id','out_enc_ir_id', 'outpatient_visit_date', 'outpt_id')

inpt = inpt[,c('V1','V4','inpt_ind')]
inpt = inpt[-1,]
outpt = outpt[,c('V1','V4','outpt_ind')]
outpt = outpt[-1,]
colnames(outpt) = c('MRN','date','event_ind')
colnames(inpt) = c('MRN','date','event_ind')
inout = rbind(inpt, outpt)
#find visits happened only after the diangosis of lupus 
pt = x[,c('MRN','DXDATE')]
df = merge(pt, inout, by = 'MRN', all.y = T)
df$DXDATE <- as.Date(df$DXDATE)
df$date <- as.Date(df$date) 
df$visitDX_timediff = df$date - df$DXDATE
df = df[df$visitDX_timediff>0,]
df = na.omit(df)
###cluster on patients who only have inpatient or outpatient visits after the dx of lupus to make sure that 
#patients still visit northwestern hospital 
MRN = unique(df$MRN)
x = x[x$MRN %in% MRN, ]
#mydata = x[,c(1,3,4,7:17,28)]
mydata = x[,c(1,2,3,5:13,15, 16)]
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


####################################################
# using the whole dataset for clustering 
mydata1 = mydata 
mydata = mydata[,-1]
mydata = sapply(mydata, function(x) as.factor(x))
mydata = as.data.frame(mydata)

f<-with(mydata, cbind(rash, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER, DXAGE1, SEX, RACE)~1) 

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

lc2<-poLCA(f, data=mydata, nclass=2, na.rm = FALSE, nrep=30, maxiter=3000)
LCA_best_model = lc2 

kk = cbind(LCA_best_model$predclass, x)
kk$`LCA_best_model$predclass` = as.factor(kk$`LCA_best_model$predclass`)

# get time-to-event data and event indicator for survial analysis 
df2 = df %>%
  group_by(MRN) %>%
  mutate(my_ranks = order(order(-as.numeric(event_ind), visitDX_timediff,  decreasing=FALSE)))

df2 = as.data.frame(df2)
df2 = df2[df2$my_ranks == 1, ]

for (i in (1:dim(df2)[1]))
{
  if (df2[i,'event_ind'] == 0){
    df2[i,'visitDX_timediff'] = as.Date('2018-10-29') - as.Date(df2$DXDATE[i])
      
  }
}

# merge kk and df; kaplan meier curve 
big = merge(kk, df2, by = 'MRN')
library(survival)
library(survminer)
colnames(big)[2] = 'class'
fit <- survfit(Surv(visitDX_timediff, event_ind) ~ class, data = big)
ggsurvplot(fit, data = big, risk.table = TRUE, pval = TRUE, conf.int = TRUE )

