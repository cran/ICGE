\name{dmahal}
\alias{dmahal}
\title{ Mahalanobis Distance}
\description{
  \code{dmahal} computes and returns the Mahalanobis distance matrix between the rows of a data matrix.
}
\usage{
dmahal(datos, S)
}
\arguments{
  \item{datos}{ data matrix.}
  \item{S}{ covariance matrix.}
}

\value{
   A \code{\link{dist}} object with distance information.
}
\references{ Pe\~nA D. (2002) \emph{Analisis de Datos Multivariantes}. Mc Graw Hill. 

Everitt B. S. and Dunn G. (2001) \emph{ Applied Multivariate Data Analysis}. 2 edition, Edward Arnold, London.
}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.eus}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV/EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\note{There is a function \code{mahalanobis()} in stats package which can perform the Mahalanobis distance. While \code{mahalanobis()} calculates the Mahalanobis distance with respect  to given a center, function  \code{dmahal()} is designed to calculate the distance between each pair of units given a data matrix.
}
\seealso{ \code{\link{dist}}, \code{\link{dbhatta}},
  \code{\link{dgower}}, \code{\link{dcor}}, \code{\link{dproc2}} }
\examples{
#Generate 10 objects in dimension 2
mu <- rep(0, 2)
Sigma <- matrix(c(10,3,3,2),2,2)

x <- mvrnorm(n=10, rep(0, 2), Sigma)

d <- dmahal(x, Sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
