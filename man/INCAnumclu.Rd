\name{INCAnumclu}
\alias{INCAnumclu}
\alias{plot.incanc}
\alias{print.incanc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Estimation of Number of Clusters in Data }
\description{
 \code{INCAnumclu} helps to estimate the number of clusters in a
 dataset. The INCA index associated to different partitions with
 different number of clusters is calculated.
}
\usage{
INCAnumclu(d, K, method = "pam", pert, L= NULL, noise=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ a distance matrix or a \code{\link{dist}} object with distance information between units.}
  \item{K}{ the maximum number of cluster to be considered. For each k value ( k=2,..,K)  a partition with k clusters is calculated. }
  \item{method}{ character string defining the clustering method in
    order to obtain the partitions. The clustering method is performed
    via the functions \code{\link{pam}}, \code{\link{agnes}},
    \code{\link{diana}} and  \code{\link{fanny}} in 
package \pkg{\link{cluster}}. The available clustering methods are \code{pam} (default method), \code{average} (UPGMA), \code{single} (single linkage),
\code{complete} (complete linkage), \code{ward} (Ward's method),
\code{weighted} (weighted average linkage), \code{diana} (hierarchical
divisive) and \code{fanny} (fuzzy clustering). Nevertheless, the user can introduce
particular or custom  partitions indicating \code{method=partition} and specifying
the partitions in argument \code{pert}.}
  \item{pert}{only useful when parameter \code{method}="partition"; it
    is a matrix and each column contains a partition of the units. That means that each column is an n-vector that indicates which group each unit belongs to. Note that the expected values of each column of \code{pert} are numbers greater than or equal to 1 (for instance 1,2,3,4..., k). }
  \item{L}{ default value NULL, but when some units are considered by
    the user as noise units, \code{L} must be specified as follows: (a)
    \code{L} is greater than or equal to 1 and all units in clusters with a
    cardinal <= L  are considered noise units; (b) \code{L="custom"} when the user wants to specify which units are considered noise units. These units must be specified in argument \code{noise}.}
  \item{noise}{when \code{L="custom"}, it is a logical vector indicating
    the units considered by the user as noise units.}
}

\value{
 Returns an object of class \code{incanc} which is a numeric vector containing the INCA index associated to each of the k (k=2,...,K) partitions. When \code{noise} is no null, the function returns a list with the INCA index for each partition, which is calculated without noise units as well as with noise units. The associated \code{plot} returns INCA index plot, both, with and without noise.
}
\references{ Irigoien, I. and  Arenas, C. (2008). INCA: New statistic for estimating the number of clusters and identifying atypical units.
\emph{Statistics in Medicine}, \bold{27}(15), 2948--2973.


Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.
 }
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}

\seealso{  \code{\link{INCAindex}}, \code{\link{estW}} }
\examples{
#------- Example 1 --------------------------------------
#generate 3 clusters, each of them with 20 objects in dimension 5.
mu1 <- sample(1:10, 5, replace=TRUE)
x1 <- matrix(rnorm(20*5, mean = mu1, sd = 1),ncol=5, byrow=TRUE)
mu2 <- sample(1:10, 5, replace=TRUE)
x2 <- matrix(rnorm(20*5, mean = mu2, sd = 1),ncol=5, byrow=TRUE)
mu3 <- sample(1:10, 5, replace=TRUE)
x3 <- matrix(rnorm(20*5, mean = mu3, sd = 1),ncol=5, byrow=TRUE)
x <- rbind(x1,x2,x3)

# calculte euclidean distance between them
d <- dist(x)

# calculate the INCA index associated to partitions with k=2, ..., k=5 clusters.
INCAnumclu(d, K=5)
out <- INCAnumclu(d, K=5)
plot(out)

#------- Example 1 cont. --------------------------------
# With hypothetical noise elements
noiseunits <- rep(FALSE, 60)
noiseunits[sample(1:60, 20)] <- TRUE
out <- INCAnumclu(d, K=5, L="custom", noise=noiseunits)
plot(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ cluster}
