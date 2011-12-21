\name{dcor}
\alias{dcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Correlation Distance }
\description{
 \code{dcor} computes and returns the Correlation distance matrix between the rows of a data matrix. This distance is defined by \eqn{d=\sqrt{1-r}}.
}
\usage{
dcor(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a  numeric matrix. }
}

\value{
  A \code{\link{dist}} object with distance information.
}

\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\references{ Gower, J.C. (1985). Measures of similarity, dissimilarity and distance. 
In: \emph{Encyclopedia of Statistical Sciences}, volume \bold{5}, 397--405. J. Wiley and Sons.}
\seealso{  \code{\link{dist}}, \code{\link{dmahal}},
  \code{\link{dgower}}, \code{\link{dbhatta}}, \code{\link{dproc2}} }
\examples{
#Generate 10 objects in dimension 8
n <- 10
mu <- sample(1:10, 8, replace=TRUE)
x <- matrix(rnorm(n*8, mean=mu, sd=1), nrow=n, byrow=TRUE)


# Correlation distances between pairs 
d <- dcor(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

