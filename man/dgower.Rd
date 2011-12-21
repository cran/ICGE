\name{dgower}
\alias{dgower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gower Distance for Mixed Variables}
\description{
  \code{dgower} computes and returns the Gower distance matrix for mixed variables. 
}
\usage{
dgower(x, type = list())
}
\arguments{
  \item{x}{ data matrix.}
  \item{type}{ it is a list with components \code{cuant}, \code{bin}, \code{nom}. Each component
indicates the column position of the quantitative, binary or nominal variables, respectively.}
}
\details{
  The distance  between two
pairs of objects i and j is obtained as
\eqn{\sqrt{2(1-s_{ij})}} where \eqn{s_{ij}} is the Gower's similarity coefficient for mixed data. This function allows
to include missing values (as \code{NA}) and therefore calculates distances based on  Gower's weighted similarity coefficient.
}
\value{
  A \code{\link{dist}} object with distance information.
}
\references{ Gower, J.C. (1971). A general coefficient of similarity and some of its properties. \emph{Biometrics}, \bold{27}, 857--871.}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\note{There is the function \code{daisy()} in \code{cluster} package which can perform the Gower distance for mixed variables. The difference is that in \code{daisy()} the distance is calculated as \eqn{d(i,j)=1-s_{ij}} and in  \code{dgower()} it is calculated as \eqn{dij=sqrt(1-s_{ij})}. 
}
\seealso{ \code{\link{dist}}, \code{\link{dmahal}},
  \code{\link{dbhatta}}, \code{\link{dcor}}, \code{\link{dproc2}} }
\examples{
#Generate 10 objects in dimension 6
# Quantitative variables
mu <- sample(1:10, 2, replace=TRUE)
xc <- matrix(rnorm(10*2, mean = mu, sd = 1), ncol=2, byrow=TRUE)

# Binary variables
xb <- cbind(rbinom(10, 1, 0.1), rbinom(10, 1, 0.5), rbinom(10, 1, 0.9))

# Nominal variables
xn <- matrix(sample(1:3, 10, replace=TRUE), ncol=1)

x <- cbind(xc, xb, xn)

# Distances
d <- dgower(x, type=list(cuant=1:2, bin=3:5, nom=6))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

