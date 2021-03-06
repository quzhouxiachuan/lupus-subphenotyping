setwd('R:/IPHAM/CHIP/Projects/SLE/SLE data and codebooks/')
x = read.csv('./LUPUS_DATA2 16April2018.csv')
#x$DATE <- as.Date(strptime(x$DATE, "%m/%d/%Y")) 
#x$DATE2 <- as.Date(strptime('06/01/2000', "%m/%d/%Y")) 
#x$datediff = x$DATE2 - x$DATE
#x = x[x$datediff<=700,]
#x = x[x$datediff>=-700,]
#x = x[!is.na(x$datediff),]
#length(unique(x$MRN))

x = x[!duplicated(x$MRN),]
missing = sapply(x, function(x) sum(is.na(x)))
#find columns that has missing data less than 30%
missing = missing[missing<309]
x = x[,colnames(x) %in% c('MRN','TOBACUSE','FATIGUE','THRALGIA',
                          'JOINTSWE','ALOPECIA','MTHULCER','PHOTOS','MAL_RASH','DIS_RASH','HEADACHE',
                          'STROKE','PLEURISY','RAYNAUDS','THROMBUS','ASA','PREDNISO','PULSE_ST',
                          'PLAQUENL','DAPSONE','CYTOXAN','PULSE_CY','IMURAN','METHOTRX','REMIT',
                          'NSAID','VASODIL','ANTI_HYP','CARDIAK','GI_AGENT','ANTIBIOT','THYRHORM',
                          'DIABETES','CYCLSPRN')]

x = subset(x, select=-c(PREDNISO,ASA,PULSE_ST,PLAQUENL,DAPSONE,CYTOXAN,PULSE_CY,IMURAN,METHOTRX,CYCLSPRN)) 


library(mice)
dat = x[,-1] 
init = mice(dat, maxit=0) 
meth = init$method
predM = init$predictorMatrix
#predM[c('MRN')]=0


meth[c("TOBACUSE",'JOINTSWE','MTHULCER','PHOTOS','MAL_RASH','DIS_RASH','PLEURISY','RAYNAUDS','ANTI_HYP','CARDIAK','ANTIBIOT','THYRHORM')]="logreg" 
meth[c("FATIGUE",'THRALGIA','HEADACHE','THROMBUS','REMIT','NSAID','VASODIL','GI_AGENT','STROKE','ALOPECIA','DIABETES')]="polyreg"
set.seed(103)
imputed = mice(dat, method=meth, predictorMatrix=predM, m=5)



