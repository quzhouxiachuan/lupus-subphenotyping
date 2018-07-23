##  https://stats.stackexchange.com/questions/130974/how-to-use-both-binary-and-continuous-variables-together-in-clustering

library(cluster)  # we'll use these packages
library(fpc)
setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
x = read.csv('edw_cld_acr_071218_DENG.csv')
x = x[,1:12]
x$cld_acr_malar = as.factor(x$cld_acr_malar)
x$cld_acr_discoid = as.factor(x$cld_acr_discoid)
x$cld_acr_photo = as.factor(x$cld_acr_photo)
x$cld_acr_ulcer = as.factor(x$cld_acr_ulcer)
x$cld_acr_arthritis = as.factor(x$cld_acr_arthritis)
x$cld_acr_serositis = as.factor(x$cld_acr_serositis)
x$cld_acr_renal = as.factor(x$cld_acr_renal)
x$cld_acr_neuro = as.factor(x$cld_acr_neuro)
x$cld_acr_heme = as.factor(x$cld_acr_heme)
x$cld_acr_ana = as.factor(x$cld_acr_ana)
x$cld_acr_immun = as.factor(x$cld_acr_immun)
data = x 

#################################
###model development, PAM clustering 
#####################################
set.seed(3296)    # this makes the example exactly reproducible
# this returns the distance matrix with Gower's distance:  
g.dist = daisy(data, metric="gower", type=list(symm=2))
pc = pamk(g.dist, krange=1:5, criterion="asw")
pc[2:3]
#$nc
#[1] 2
#$crit
#[1] 0.0000000 0.2846817 0.2441667 0.2359307 0.2197452
pc = pc$pamobject;  
pam=pam(g.dist,k=2,diss=TRUE)
#get average silinfo score
asw=pam$silinfo$avg.width
pdf("lupus_2clusters.pdf")
plot(pam)
dev.off()

################################################
#######hierarchical clustering 
################################################
hc.m = hclust(g.dist, method="median")
plot(hc.m)


###############################################
########polCA latent class analysis 
#######################################










