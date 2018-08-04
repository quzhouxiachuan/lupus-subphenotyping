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
x = read.csv('cld_full_acr.csv')
x = x[,6:16] 
x$rash = x$ACRMAL + x$ACRDISC
x$rash[x$rash>1] =1 
x = x[,3:12]
#x= sapply(x, function(x) as.factor(x))
# split x into training and testing dataset 
set.seed(234511)
smp_size <- floor(0.5 * nrow(x))
train_ind <- sample(seq_len(nrow(x)), size = smp_size)
train <- x[train_ind, ]
test <- x[-train_ind, ]

###################################################
# using train dataset for clustering 
mydata = train + 1  #514
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

set.seed(01125)
lc1<-poLCA(f, data=mydata, nclass=1, na.rm = FALSE, nrep=30, maxiter=3000) #Loglinear independence model.
lc2<-poLCA(f, data=mydata, nclass=2, na.rm = FALSE, nrep=30, maxiter=3000)
lc3<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000)
lc4<-poLCA(f, data=mydata, nclass=4, na.rm = FALSE, nrep=30, maxiter=3000) 
lc5<-poLCA(f, data=mydata, nclass=5, na.rm = FALSE, nrep=30, maxiter=3000)
lc6<-poLCA(f, data=mydata, nclass=6, na.rm = FALSE, nrep=30, maxiter=3000)
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



library("forcats")
#results$Modell < - as_factor(results$Modell) 

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
set.seed(01125)
lc_train<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)


##################################################
###########use testing dataset ###########
########################################################
mydata = test + 1  #515
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

set.seed(12345)
lc1<-poLCA(f, data=mydata, nclass=1, na.rm = FALSE, nrep=30, maxiter=3000) #Loglinear independence model.
lc2<-poLCA(f, data=mydata, nclass=2, na.rm = FALSE, nrep=30, maxiter=3000)
lc3<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000)
lc4<-poLCA(f, data=mydata, nclass=4, na.rm = FALSE, nrep=30, maxiter=3000) 
lc5<-poLCA(f, data=mydata, nclass=5, na.rm = FALSE, nrep=30, maxiter=3000)
lc6<-poLCA(f, data=mydata, nclass=6, na.rm = FALSE, nrep=30, maxiter=3000)
results <- data.frame(Modell=c("Modell 1"),
                      log_likelihood=lc1$llik,
                      df = lc1$resid.df,
                      BIC=lc1$bic,
                      AIC=  lc1$aic,
                     # CAIC = (-2*lc1$llik) + lc1$npar * (1 + log(lc1$N)), 
                      likelihood_ratio=lc1$Gsq)
results$Modell<-as.integer(results$Modell)
results[1,1]<-c("1")
results[2,1]<-c("2")
results[3,1]<-c("3")
results[4,1]<-c("4")
results[5,1]<-c("5")
results[6,1]<-c("6")
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

results[2,5]<-lc2$aic
results[3,5]<-lc3$aic
results[4,5]<-lc4$aic
results[5,5]<-lc5$aic
results[6,5]<-lc6$aic

#results[2,6]<- (-2*lc2$llik) + lc2$npar * (1 + log(lc2$N)) #caic
#results[3,6]<- (-2*lc3$llik) + lc3$npar * (1 + log(lc3$N))
#results[4,6]<- (-2*lc4$llik) + lc4$npar * (1 + log(lc4$N))
#results[5,6]<- (-2*lc5$llik) + lc5$npar * (1 + log(lc5$N))
#results[6,6]<- (-2*lc6$llik) + lc6$npar * (1 + log(lc6$N))

results[2,6]<-lc2$Gsq
results[3,6]<-lc3$Gsq
results[4,6]<-lc4$Gsq
results[5,6]<-lc5$Gsq
results[6,6]<-lc6$Gsq



library("forcats")
#results$Modell < - as_factor(results$Modell) 

#convert to long format
results2<-tidyr::gather(results,Kriterium,Guete,4:6)
results2

#plot
fit.plot<-ggplot(results2) + 
  geom_point(aes(x=Modell,y=Guete),size=3) +
  geom_line(aes(Modell, Guete, group = 1)) +
  theme_bw()+
  labs(x = "number of clusters", y="", title = "LCA Model Assessment") + 
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
set.seed(12345)
lc_test<-poLCA(f, data=mydata, nclass=3, na.rm = FALSE, nrep=30, maxiter=3000,graph=TRUE)







########################################################################
#############################################################################





#extract starting values from our previous best model (with 3 classes)
probs.start<-lc_train$probs.start

#re-run the model, this time with "graphs=TRUE"
#lc<-poLCA(f, mydata, nclass=3, probs.start=probs.start,graphs=TRUE, na.rm=TRUE, maxiter=3000)

# If you don´t like the order, reorder them (here: Class 1 stays 1, Class 3 becomes 2, Class 2 becomes 1)
new.probs.start<-poLCA.reorder(probs.start, c(2,1,3))

#run polca with adjusted ordering
lc<-poLCA(f, mydata, nclass=3, probs.start=new.probs.start,graphs=TRUE, na.rm=TRUE)
lc
#######################################################################
#######################use 3D plot for visualization #############

data = train 
set.seed(32961)    # this makes the example exactly reproducible
# this returns the distance matrix with Gower's distance:  
g.dist = daisy(data, metric="gower", type=list(symm=1:10))
pc = pamk(g.dist, krange=1:10, criterion="asw")
pc[2:3]
#$nc
#[1] 2
#$crit
#[1] 0.0000000 0.2846817 0.2441667 0.2359307 0.2197452
pc = pc$pamobject;  
pam=pam(g.dist,k=3,diss=TRUE)
#get average silinfo score
asw=pam$silinfo$avg.width
plot(pam)

############data summary #############
pam_results <- x %>%
  #dplyr::select(-MRN) %>%
  mutate(cluster = pam$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
#####################visualization ###########
tsne_obj <- Rtsne(g.dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam$clustering))


tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(lc_test$predclass))


ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

tsne_obj <- Rtsne(g.dist, is_distance = TRUE, dims= 3)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y",'Z')) %>%
  mutate(cluster = factor(lc_test$predclass))

p <- plot_ly(tsne_data, x = ~X, y = ~Y, z = ~Z, color = ~cluster, colors = c('#BF382A', '#0C4B8E', 'black'),size=I(1)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X'),
                      yaxis = list(title = 'Y'),
                      zaxis = list(title = 'Z')))

p





