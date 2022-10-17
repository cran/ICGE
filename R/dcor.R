dcor<- function(x){
########## Correlation distance between pairs of objects ##################
# Input:
# x: data matrix
#
# Output:
# d: distance matrix
###########################################################################
    d <- sqrt(1-cor(t(x)))
    d <- as.dist(d)
    return(d)
}
