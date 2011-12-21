\name{dproc2}
\alias{dproc2}
\title{ Modified Procrustes distance}
\description{
 \code{dproc2} computes and returns all the pairwise procrustes distances between genes in a time course experiment, using their expression profile.
}
\usage{
dproc2(x, timepoints = NULL)
}
\arguments{
  \item{x}{ a  matrix containing, in its rows, the gene expression values at the T considered time points.}
  \item{timepoints}{a T-vector with the T observed time points. If  \code{timepoints=NULL} (default), then timepoints=1:T.}
}
\details{
 Each row i of matrix x is arranged in a two column matrix Xi. In Xi, the first column contains the time points and the second column the observed gene expression values (xi1...). 
}
\value{
 A \code{\link{dist}} object with distance information.
}
\references{Irigoien, I. , Vives, S. and Arenas, C. (2011). Microarray
  Time Course Experiments: Finding Profiles. \emph{IEEE/ACM Transactions
  on Computational Biology and Bioinformatics}, \bold{8}(2), 464--475.
  
  Gower, J. C. and Dijksterhuis, G. B. (2004)
  \emph{Procrustes Problems}. Oxford University Press.
  
  Sibson, R. (1978). Studies in the Robustness of Multidimensional
    Scaling: Procrustes statistic. \emph{Journal of the Royal Statistical
  Society, Series B}, \bold{40}, 234--238.
}

\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\seealso{\code{\link{dist}}, \code{\link{dmahal}}, \code{\link{dgower}}, \code{\link{dcor}} \code{\link{dbhatta}} }
\examples{
# Given  10  hypothetical time course profiles
# over 6 time points at 1, 2, ..., 6 hours.
x <- matrix(c(0.38, 0.39, 0.38, 0.37, 0.385, 0.375,
              0.99, 1.19, 1.50, 1.83, 2.140, 2.770,
              0.38, 0.50, 0.71, 0.72, 0.980, 1.010,
              0.20, 0.40, 0.70, 1.06, 2.000, 2.500,
              0.90, 0.95, 0.97, 1.50, 2.500, 2.990,
              0.64, 2.61, 1.51, 1.34, 1.330 ,1.140,
              0.71, 1.82, 2.28, 1.72, 1.490, 1.060,
              0.71, 1.82, 2.28, 1.99, 1.975, 1.965,
              0.49, 0.78, 1.00, 1.27, 0.590, 0.340,
              0.71,1.00, 1.50, 1.75, 2.090, 1.380), nrow=10, byrow=TRUE)

# Graphical representation
matplot(t(x), type="b")

# Distance matrix between them
d <- dproc2(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
