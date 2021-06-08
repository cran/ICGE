INCAncinner <- function(d, K, method,pert){
############### Number of atypical individuals  #############################
# Input:
# d: distance matrix (n x n)
# K: partitions in 2:K clusters are considered
# method: method to perform the clustering. If method="partition" a collection
#         of M partitions must by given in pert.
# 
# Output:
# Totals : INCA index for k=2, ..., K
#############################################################################

d <- as.matrix(d)
n <- dim(d)[1]


pama <- function(d, K){
  ppp <- sapply(2:K, function(x){cluster::pam(d, k=x, diss=TRUE)$clustering})
  T <- c(NA, apply(ppp, 2, function(x){INCAindex(d, x)$Total}))
  return(T)
}  

dian <- function(d, K){
  aux <- cluster::diana(d, diss=TRUE)
  ppp <- sapply(2:K, function(x){stats::cutree(aux, k=x)})
  T <- c(NA, apply(ppp, 2, function(x){INCAindex(d, x)$Total}))
  return(T)
}  


fannya <- function(d, K){
  ppp <- sapply(2:K, function(x){cluster::fanny(d, k=x, diss=TRUE)$clustering})
  T <- c(NA, apply(ppp, 2, function(x){INCAindex(d, x)$Total}))
  return(T)
} 

bestela <- function(d, K, method){
   aux <- fastcluster::hclust(as.dist(d), method=method)
   ppp <- sapply(2:K, function(x){stats::cutree(aux, k=x)})
   T <- c(NA, apply(ppp, 2, function(x){INCAindex(d, x)$Total}))
   return(T)
}   

if (method=="partition"){
    pert <- as.matrix(pert)
    M <- dim(pert)[2]
    K <- max(pert)
    #      0 can not be a partitions name
    if (any(pert == 0))
      stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")
    T <- rep(NA, K)
#    tops <- rep(NA, M)
    for (m in 1:M){
       tops <- max(pert[,m])
#      populations must be named with numbers from 1 to k
       if (length(tabulate(as.factor(pert[,m]))) != tops)
          stop("Partitions must be named by factors or with numbers from 1 to k.")
       T[tops] <- INCAindex(d, pert[,m])$Total
    }  
} else{
   T <- switch(method,  pam=pama(d, K), diana=dian(d, K), fanny=fannya(d, K), bestela(d,K, method))
}
                                          

                                                                
                                                               
out <- list(INCAindex=T, method=method)
#class(out) <- "incanc"
return(out)

} #end of function
