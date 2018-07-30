library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")
setwd('R:/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')

x = read.csv('cld_full_acr.csv')
x = x[,6:16] 
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
x = x[,3:12]
x= sapply(x, function(x) as.factor(x))
data = x 
data = as.data.frame(data)

set.seed(32961)    # this makes the example exactly reproducible
# this returns the distance matrix with Gower's distance:  
g.dist = daisy(data, metric="gower", type=list(symm=1:10))
pc = pamk(g.dist, krange=1:30, criterion="asw")
pc[2:3]


tsne_obj <- Rtsne(g.dist, is_distance = TRUE, dims= 3)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y",'Z')) %>%
  mutate(cluster = factor(lc3$predclass))

p <- plot_ly(tsne_data, x = ~X, y = ~Y, z = ~Z, color = ~cluster, colors = c('#BF382A', '#0C4B8E','black'),size=I(2)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X'),
                      yaxis = list(title = 'Y'),
                      zaxis = list(title = 'Z')))

p





##latent class analysis 


library("poLCA")

# By the way, for all examples in this article, you´ll need some more packages:
library("reshape2")
library("plyr")
library("dplyr")
library("poLCA")
library("ggplot2")

x = read.csv('cld_full_acr.csv')
x = x[,6:16] 
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
mydata = x + 1 

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
 set.seed(01012)
 lc3<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000) #Loglinear independence model.



