#this code is adapted from https://carsonfarmer.com/2010/04/parallel-bootstrapping-with-r/ 

library(fpc)
library(poLCA)
library(doMC)
library(foreach)
# Jaccard coeficcient function (taken from package fpc)
clujaccard = function (c1, c2, zerobyzero = NA) {
  if (sum(c1) + sum(c2) - sum(c1 & c2) == 0)
    out = zerobyzero
  else
    out = sum(c1 & c2)/(sum(c1) + sum(c2) - sum(c1 & c2))
  return(out)
}


data(values)# load a dataset
g= values 
f <- cbind(A,B,C,D)~1 # compute original clustering
c1 <- poLCA(f,values,nclass=2) # compute original clustering
#noc = length(unique(moc)) # count the number original clusters
noc = 2 
bg = g # make a copy of g for bootstrapping
#clusters = foreach(i=seq(B), .combine=cbind) %dopar% {
#for (iloop in 1:200)
kk = data.frame()
for (iloop in 1:4)
{
  bsamp <- sample(1:nrow(g), replace = TRUE)  # resample the edge weights
  bg <- g[bsamp,]
  bc1 <- poLCA(f,bg,nclass=2, verbose = F) # compute bootstrapping clusters
  mbc = bc1$predclass # get membership for bootstrapping clusters
  moc = c1$predclass[bsamp] # get membership for orginal clustering 
  noc = length(unique(moc))
  nbc = length(unique(mbc)) # count the number new clusters
  bootresult = c()
  
  for (j in seq(1, noc)) { # for each of the original clusters...
    
    maxgamma = 0
    # print(j)}
    if (nbc > 0) {
      for (k in seq(1, nbc)) { # for each of the new clusters...
        bv = as.vector(mbc == k) #new cluster 
        ov = as.vector(moc == j) #orginal cluster
        jc = clujaccard(bv, ov, zerobyzero = 0)
        if (jc > maxgamma) # if these two clusters are most similar...
          {maxgamma = jc}
      }
      
    }
    bootresult = c(bootresult, maxgamma) # combine results
    print (bootresult)
    }
   kk = rbind(kk, bootresult)  # return the results of this iteration (and cbind with the rest)
} 
