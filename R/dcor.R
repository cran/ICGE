dcor<- function(x){
########## Correlation distance between pais of objects ##################
# Input:
# x: data matrix
#
# Output:
# d: distance matrix
###########################################################################
    d <- sqrt(1-cor(x))
    d <- as.dist(d)
    return(d)
}
