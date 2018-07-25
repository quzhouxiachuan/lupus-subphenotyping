library(cluster)  # we'll use these packages
library(fpc)
library(Rtsne)
library(dplyr) # for data cleaning

setwd('/Volumes/fsmresfiles/IPHAM/CHIP/Data_Team/Projects/Lupus/Yu')
x = read.csv('edw_cld_acr_071218_DENG.csv')
x = x[,1:12]

#x= sapply(x, function(x) as.factor(x)))
       
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
data = data[,-1]
#################################
###model development 
set.seed(32961)    # this makes the example exactly reproducible
# this returns the distance matrix with Gower's distance:  
g.dist = daisy(data, metric="gower", type=list(symm=1:11))
pc = pamk(g.dist, krange=1:10, criterion="asw")
pc[2:3]
#$nc
#[1] 2
#$crit
#[1] 0.0000000 0.2846817 0.2441667 0.2359307 0.2197452
pc = pc$pamobject;  
pam=pam(g.dist,k=2,diss=TRUE)
#get average silinfo score
asw=pam$silinfo$avg.width
plot(pam)

############data summary #############
pam_results <- x %>%
  dplyr::select(-MRN) %>%
  mutate(cluster = pam$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
#####################visualization ###########
tsne_obj <- Rtsne(g.dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam$clustering),
         name = x$MRN)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

tsne_obj <- Rtsne(g.dist, is_distance = TRUE, dims= 3)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y",'Z')) %>%
  mutate(cluster = factor(pam$clustering),
         name = x$MRN)

p <- plot_ly(tsne_data, x = ~X, y = ~Y, z = ~Z, color = ~cluster, colors = c('#BF382A', '#0C4B8E'),size=I(3)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X'),
                      yaxis = list(title = 'Y'),
                      zaxis = list(title = 'Z')))

p




hc.m = hclust(g.dist, method="median")
plot(hc.m)
s=NULL
for(i in 2:10){
  k=cutree(hclust(g.dist, method="median"),i)
  s[i]=summary(silhouette(k,dist(data)))$si.summary[3]
}
######################################
#######gower clustering 
#########################################






















#################################################
##try new dataset 
#############################################
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













