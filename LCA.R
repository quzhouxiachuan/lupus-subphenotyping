#install.packages("poLCA")
library("poLCA")

# By the way, for all examples in this article, you´ll need some more packages:
library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
#library("ggparallel")
#library("igraph")
#library("tidyr")

# select variables
setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
x = read.csv('edw_cld_acr_071218_DENG.csv')
mydata <- x[,-1]
mydata = mydata+1 
# define function
#f<-with(mydata, cbind(cld_acr_malar, cld_acr_discoid, cld_acr_photo, cld_acr_ulcer,    
 #                     cld_acr_arthritis, cld_acr_serositis, cld_acr_renal, cld_acr_neuro,   
 #                    cld_acr_heme, cld_acr_ana,  cld_acr_immun)~1) 


f<-with(mydata, cbind(ACRMAL, ACRDISC, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER)~1) 

#------ run a sequence of models with 1-10 classes and print out the model with the lowest BIC
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




#select data
#mydata <- x[,-1]

# define function
#f<-with(mydata, cbind(cld_acr_malar, cld_acr_discoid, cld_acr_photo, cld_acr_ulcer,    
 #                     cld_acr_arthritis, cld_acr_serositis, cld_acr_renal, cld_acr_neuro,   
  #                    cld_acr_heme, cld_acr_ana,  cld_acr_immun)~1)

f<-with(mydata, cbind(ACRMAL, ACRDISC, ACRPHOTO, ACRJOINT,    
                      ACRSERO, ACRRENAL, ACRNEURO, ACRHEME,   
                      ACRANA, ACRIMMUN, ACRULCER)~1) 



## models with different number of groups without covariates:
set.seed(01012)
lc1<-poLCA(f, data=mydata, nclass=1, na.rm = FALSE, nrep=30, maxiter=3000) #Loglinear independence model.
lc2<-poLCA(f, data=mydata, nclass=2, na.rm = FALSE, nrep=30, maxiter=3000)
lc3<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000)
lc4<-poLCA(f, data=mydata, nclass=4, na.rm = FALSE, nrep=30, maxiter=3000) 
lc5<-poLCA(f, data=mydata, nclass=5, na.rm = FALSE, nrep=30, maxiter=3000)
lc6<-poLCA(f, data=mydata, nclass=6, na.rm = FALSE, nrep=30, maxiter=3000)

# generate dataframe with fit-values

results <- data.frame(Modell=c("Modell 1"),
                      log_likelihood=lc1$llik,
                      df = lc1$resid.df,
                      BIC=lc1$bic,
                      ABIC=  (-2*lc1$llik) + ((log((lc1$N + 2)/24)) * lc1$npar),
                      CAIC = (-2*lc1$llik) + lc1$npar * (1 + log(lc1$N)), 
                      likelihood_ratio=lc1$Gsq)
results$Modell<-as.integer(results$Modell)
results[1,1]<-c("Modell 1")
results[2,1]<-c("Modell 2")
results[3,1]<-c("Modell 3")
results[4,1]<-c("Modell 4")
results[5,1]<-c("Modell 5")
results[6,1]<-c("Modell 6")

results[2,2]<-lc2$llik
results[3,2]<-lc3$llik
results[4,2]<-lc4$llik
results[5,2]<-lc5$llik
results[6,2]<-lc6$llik

results[2,3]<-lc2$resid.df
results[3,3]<-lc3$resid.df
results[4,3]<-lc4$resid.df
results[5,3]<-lc5$resid.df
results[6,3]<-lc6$resid.df

results[2,4]<-lc2$bic
results[3,4]<-lc3$bic
results[4,4]<-lc4$bic
results[5,4]<-lc5$bic
results[6,4]<-lc6$bic

results[2,5]<-(-2*lc2$llik) + ((log((lc2$N + 2)/24)) * lc2$npar) #abic
results[3,5]<-(-2*lc3$llik) + ((log((lc3$N + 2)/24)) * lc3$npar)
results[4,5]<-(-2*lc4$llik) + ((log((lc4$N + 2)/24)) * lc4$npar)
results[5,5]<-(-2*lc5$llik) + ((log((lc5$N + 2)/24)) * lc5$npar)
results[6,5]<-(-2*lc6$llik) + ((log((lc6$N + 2)/24)) * lc6$npar)

results[2,6]<- (-2*lc2$llik) + lc2$npar * (1 + log(lc2$N)) #caic
results[3,6]<- (-2*lc3$llik) + lc3$npar * (1 + log(lc3$N))
results[4,6]<- (-2*lc4$llik) + lc4$npar * (1 + log(lc4$N))
results[5,6]<- (-2*lc5$llik) + lc5$npar * (1 + log(lc5$N))
results[6,6]<- (-2*lc6$llik) + lc6$npar * (1 + log(lc6$N))

results[2,7]<-lc2$Gsq
results[3,7]<-lc3$Gsq
results[4,7]<-lc4$Gsq
results[5,7]<-lc5$Gsq
results[6,7]<-lc6$Gsq

#####Elbow plot ########################################
# plot 1

# Order categories of results$model in order of appearance
install.packages("forcats")
library("forcats")
results$model < - as_factor(results$model) 

#convert to long format
results2<-tidyr::gather(results,Kriterium,Guete,4:7)
results2

#plot
fit.plot<-ggplot(results2) + 
  geom_point(aes(x=Modell,y=Guete),size=3) +
  geom_line(aes(Modell, Guete, group = 1)) +
  theme_bw()+
  labs(x = "", y="", title = "") + 
  facet_grid(Kriterium ~. ,scales = "free") +
  theme_bw(base_size = 16, base_family = "") +   
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour="grey", size=0.5),
        legend.title = element_text(size = 16, face = 'bold'),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text=  element_text(size=16),
        axis.line = element_line(colour = "black")) # Achsen etwas dicker

# save 650 x 800
fit.plot


lc2<-poLCA(f, data=mydata, nclass=2, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)









