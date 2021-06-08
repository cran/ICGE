dbhatta <- function(x){
########## Bhattacharyya distance #######################
# Input:
# x: data matrix
# Output:
# d: bhattacharyya distance matrix
#########################################################
	x <- as.matrix(x)
	aux <- apply(x, 1, sum)
	if (any(aux != 1))
	{
	  stop("Rows must add up to 1.")
	}
	
	x1 <- sqrt(x)
	d <- acos(x1%*%t(x1))
	d <- as.dist(d)
	return(d)
}
