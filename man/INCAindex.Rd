\name{INCAindex}
\alias{INCAindex}
\alias{print.incaix}
\alias{summary.incaix}
\alias{plot.incaix}
\title{ INCA index}
\description{
  \code{INCAindex} helps to estimate the number of clusters in a dataset.
}
\usage{
INCAindex(d, pert_clus)
}
\arguments{
  \item{d}{ a distance matrix or a  \code{dist} object with distance information between units.}
  \item{pert_clus}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are numbers  greater than or equal to 1 (for instance 1,2,3,4..., k). The default value indicates the presence of only one group in data.}
}
\value{Returns an object of class \code{incaix} which is a list containing the following components:
  \item{well_class}{a vector indicating the number of well classified units.}
  \item{Ni_cluster}{a vector indicating each cluster size.}
  \item{Total}{percentage of objects well classified in the partition defined by \code{pert_clus}.}
}
\references{ Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.

 Irigoien, I. and  Arenas, C. (2008). INCA: New statistic for estimating the number of clusters and identifying atypical units.
\emph{Statistics in Medicine}, \bold{27}(15), 2948--2973.}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.eus}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV/EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\note{For a correct geometrical interpretation it is convenient to verify whether the distance matrix d is Euclidean.
It admits the associated methods summary and plot. The first simply returns the percentage of well-classified units and the second offers a barchart with the percentages of well classified units for each group in the given partition.
}
\seealso{ \code{\link{estW}}, \code{\link{INCAtest}}}
\examples{
#generate 3 clusters, each of them with 20 objects in dimension 5.
mu1 <- sample(1:10, 5, replace=TRUE)
x1 <- matrix(rnorm(20*5, mean = mu1, sd = 1),ncol=5, byrow=TRUE)
mu2 <- sample(1:10, 5, replace=TRUE)
x2 <- matrix(rnorm(20*5, mean = mu2, sd = 1),ncol=5, byrow=TRUE)
mu3 <- sample(1:10, 5, replace=TRUE)
x3 <- matrix(rnorm(20*5, mean = mu3, sd = 1),ncol=5, byrow=TRUE)
x <- rbind(x1,x2,x3)

# Euclidean distance between units.
d <- dist(x)

# given the right partition, calculate the percentage of well classified objects.
partition <- c(rep(1,20), rep(2,20), rep(3,20))
INCAindex(d, partition)


# In order to estimate the number of cluster in data, try several 
#  partitions and compare the results
library(cluster)
T <- rep(NA, 5)
for (l in 2:5){
	part <- pam(d,l)$clustering
	T[l] <- INCAindex(d,part)$Total
}

plot(T, type="b",xlab="Number of clusters", ylab="INCA", xlim=c(1.5, 5.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
\keyword{ cluster}
