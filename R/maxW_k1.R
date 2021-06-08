maxW_k1 <- function(x_d,v_d){
######################### Used in INCAnumclu #################################
#                  When number of clusters k=1
# Input:
# x_d <- distance matrix of the partition considered as populations
# v_d <- distance from each individual of the testing cluster to the rest.
# pert <- partition of individuals considered as populations.
# Output:
# contar: number of individuals of the testing cluster considered as atipical or #         well classified.
##########################################################################
if (is.null(dim(v_d))){dim(v_d) <- c(1, length(v_d))}
nt <- dim(v_d)[1]
k <- 1
if (is.null(dim(x_d))){dim(x_d) <- c(1, length(x_d))}
nn <- dim(x_d)[1]
atipico <- matrix(0,nt,1)


# As there is an only cluster in population
pert <- matrix(1,nn,1)
vg <- vgeo(x_d, pert);


# Calculate W for ind. in xx/d_x
TW <- apply(x_d, 1, proxi_simple, var=vg, pert=pert, frec=nn)



# Calculate the maximum of W for ind. in xx
M <- max(TW)

############################################
# Calculate W for ind. in v_d and evaluter whether is atypical or not
W0 <- apply(v_d, 1, proxi_simple, var=vg, pert=pert, frec=nn)
atipico <- ifelse(W0 <= M, 0, 1)

out <- sum(atipico)

return(out)

} # end function
