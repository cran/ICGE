\name{vgeo}
\alias{vgeo}
\title{ Geometric Variability}
\description{
 Assume that n units are divided into k groups C1,...,Ck. The function calculates the geometrical variability for each group in data.
}
\usage{
vgeo(d, pert = "onegroup")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ a distance matrix or a \code{dist} object with distance information between units.}
  \item{pert}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are numbers greater than or equal to 1  (for instance 1,2,3,4..., k). The default value indicates there is only one group in data.}
}
\value{
  It is a matrix containing the geometric variability for each group.
}
\references{ Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.

Cuadras, C.M. (1992). Some examples of distance based discrimination. \emph{Biometrical Letters}, \bold{29}(1), 3--20.
}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\seealso{
\code{\link{deltas}}, \code{\link{proxi}}
}

\examples{
data(iris)
d <- dist(iris[,1:4])
vgeo(d,iris[,5])


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
