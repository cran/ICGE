estW_simple<-function(phi, var, delta, Uout=TRUE){
######################## Estimation of INCA W (faster version) ###############
# Input:
    # phi: proximity vector
    # var: vector of geometric varibilities
    # delta: matrix of distances between clusters
    # Uout: logical to indicate whether U values should be returned
# Output:
   # Wvalue: estimation of W
   # Uvalue: projections U1, ..., Uk
##############################################################################

k <- length(var)
U <-matrix(0, k, 1)
M <- matrix(0, k-1, k-1)
N <- matrix(0, k-1, 1)

for (i in 1:(k-1)){
    for (j in i:(k-1)){

        M[i,j] <- delta[i,k]+delta[j,k]-delta[i,j]
        M[j,i] <- M[i,j]
    }
    N[i]=delta[i,k]+phi[k]-phi[i]
}


alpha <- MASS::ginv(M)%*%N

aux <- sum(alpha)
alpha <- c(alpha, 1-aux)


waux <- 0
for (i in 1:(k-1)){
    for(j in (i+1):k){
        waux <- waux + alpha[i]*alpha[j]*delta[i,j]
    }
}

W <- sum(alpha*phi) - waux


if ( W < 0 ){
    W <- 0
}

U <- phi - W

if(Uout)
{
    out=list(Wvalue=W, Uvalue=U) 
}else
{
    out <- W
}

return(out)
}
