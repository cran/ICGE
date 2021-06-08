maxW_k <- function(x_d,v_d, pert){
######################### Used in INCAnumclu #################################
#                  When number of clusters k>=2
# Input: 
# x_d <- distance matrix of the partition considered as populations
# v_d <- distance from each individual of the testing cluster to the rest.
# pert <- partition of individual considered as populations.
# Output:
# contar: number of individuals of the testing cluster considered as atipical or #         well classified.
##########################################################################

if (is.null(dim(v_d))){dim(v_d) <- c(1, length(v_d))}
nt <- dim(v_d)[1]
k <- max(pert)
if (is.null(dim(x_d))){dim(x_d) <- c(1, length(x_d))}
nn <- dim(x_d)[1]

vg <- vgeo(x_d, pert);


delta <- deltas_simple(x_d,vg,pert)

frec <- tabulate(pert)

# Calculate W for individuals in xx/d_x
phi_d <- apply(x_d, 1, proxi_simple, var=vg, pert=pert, frec=frec)
TW <- apply(phi_d, 2, estW_simple, var=vg, delta= delta, Uout=FALSE)

# Calculate  maxumum of  W for ind. in  xx
M <- max(TW)

############################################
# Calculate W for ind. in  v_d and evalutate whether it is atypical or not
phi_x <- apply(v_d, 1, proxi_simple, var=vg, pert=pert, frec=frec)
W0 <- apply(phi_x, 2, estW_simple, var=vg, delta= delta, Uout=FALSE)

atipico <- ifelse(W0 <= M, 0, 1)

out <- sum(atipico)

return(out)

} # end function







