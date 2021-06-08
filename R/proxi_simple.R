proxi_simple<- function(dx0, var, pert, frec){
################### Proximity function (faster version)########################
# Input:  
# dx0: vector of dist. from one individual to the ind. of different clusters.
# var: geometric variabilities.
# pert: integer vector indicating the group each individual belongs to.
# frec: frequencies of pert
# Output: 
# phi: vector of proximities from one individual to each cluster.

n <- length(pert)
k<- max(pert)


dx0 <- dx0*dx0
phi <- matrix(0,k,1)

#frec <- tabulate(pert) # vector of frecuencies of individuals in each population


for (pob in 1:k){
    aux <- sum(dx0[pert==pob])
    phi[pob] <- aux/frec[pob] - var[pob]
}


return(phi)

}
