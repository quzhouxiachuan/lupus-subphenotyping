library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr)
setwd('R:/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
#x = read.csv('./cld_full_slicc.csv')
x = read.csv('./cld_slicc_final.csv')
x = x[,2:18]

res = data.frame() 
mydata= x + 1 
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



lc_test = LCA_best_model
res = cbind(lc_test$probs$slicc_acute[5:8],lc_test$probs$slicc_discoid[5:8],lc_test$probs$slicc_ulcer[5:8],lc_test$probs$slicc_alopecia[5:8]
            , lc_test$probs$slicc_arthritis[5:8],lc_test$probs$slicc_serositis[5:8],lc_test$probs$slicc_renal[5:8],lc_test$probs$slicc_neuro[5:8]
            ,lc_test$probs$slicc_hemanemia[5:8],lc_test$probs$slicc_leuk[5:8],lc_test$probs$slicc_throm[5:8],lc_test$probs$slicc_ana[5:8],lc_test$probs$slicc_dsdna[5:8]
            ,lc_test$probs$slicc_sm[5:8], lc_test$probs$slicc_apa[5:8],lc_test$probs$slicc_compliment[5:8],lc_test$probs$slicc_coombs[5:8],
            lc_test$P)

colnames(res) = c("slicc_acute","slicc_discoid", "slicc_ulcer", "slicc_alopecia", "slicc_arthritis" , "slicc_serositis","slicc_renal", "slicc_neuro", "slicc_hemanemia" 
                  ,"slicc_leuk", "slicc_throm", "slicc_ana" , "slicc_dsdna" , "slicc_sm", "slicc_apa" , "slicc_compliment", "slicc_coombs"  
                  , 'population_share')


res = cbind(lc_test$probs$slicc_acute[4:6],lc_test$probs$slicc_discoid[4:6],lc_test$probs$slicc_ulcer[4:6],lc_test$probs$slicc_alopecia[4:6]
            , lc_test$probs$slicc_arthritis[4:6],lc_test$probs$slicc_serositis[4:6],lc_test$probs$slicc_renal[4:6],lc_test$probs$slicc_neuro[4:6]
            ,lc_test$probs$slicc_hemanemia[4:6],lc_test$probs$slicc_leuk[4:6],lc_test$probs$slicc_throm[4:6],lc_test$probs$slicc_ana[4:6],lc_test$probs$slicc_dsdna[4:6]
            ,lc_test$probs$slicc_sm[4:6], lc_test$probs$slicc_apa[4:6],lc_test$probs$slicc_compliment[4:6],lc_test$probs$slicc_coombs[4:6],
            lc_test$P)








  
