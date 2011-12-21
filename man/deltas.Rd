\name{deltas}
\alias{deltas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Distance Between Groups}
\description{
 Assume that n units are divided into k groups C1,...,Ck . Function \code{deltas}  computes and returns the distance between each pair of groups. It uses the distances between pairs of units. 
}
\usage{
deltas(d, pert = "onegroup")
}
\arguments{
  \item{d}{ a distance matrix or a \code{dist} object with distance information between units.}
  \item{pert}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are numbers greater than or equal to 1 (for instance 1,2,3,4..., k). The default value indicates there is only one group in data.}
}
\value{
A matrix containing the distances between each pair of groups.
}
\references{ Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.

Cuadras, C.M., Fortiana, J. and Oliva, F. (1997). The proximity of an individual to a population with applications in
  discriminant analysis.
\emph{ Journal of Classification}, \bold{14}, 117--136.}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\seealso{
\code{\link{vgeo}}, \code{\link{proxi}}
}
\examples{
data(iris)
d <- dist(iris[,1:4])
deltas(d,iris[,5])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
