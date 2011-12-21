\name{proxi}
\alias{proxi}
\title{ Proximity Function}
\description{
   Assume that n units are divided into k groups C1,...,Ck. The function calculates the proximity function from a specific unit x0 to the groups Cj. 
}
\usage{
proxi(d, dx0, pert = "onegroup")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ a distance matrix or a \code{dist} object with distance information between units.}
  \item{dx0}{ an n-vector containing the distances from x0 to the other units.}
  \item{pert}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are numbers greater than or equal to 1 (for instance 1,2,3,4..., k). The default value indicates there is only one group in data.}
}
\value{
  k-vector containing the proximity function value from x0 to each group.
}
\references{ Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.

Cuadras, C.M., Fortiana, J. and Oliva, F. (1997). The proximity of an individual to a population with applications in
  discriminant analysis.
\emph{ Journal of Classification}, \bold{14}, 117--136. 
}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\seealso{
\code{\link{vgeo}}, \code{\link{deltas}}
}



\examples{
data(iris)
d <- dist(iris[,1:4])

# xo contains a unit from one group, as for example group 1.
x0 <- c(5.3, 3.6, 1.1, 0.1) 
# distances between x0 and the other units.
dx0 <- rep(0,150)
for (i in 1:150){
	dif <-x0-iris[i,1:4]
	dx0[i] <- sqrt(sum(dif*dif))
}

proxi(d, dx0, iris[,5])


# xo contains a unit from one group, as for example group 2.
x0 <- c(6.4, 3.0, 4.8, 1.3) 
# distances between x0 and the other units.
dx0 <- rep(0,150)
for (i in 1:150){
	dif <-x0-iris[i,1:4]
	dx0[i] <- sqrt(sum(dif*dif))
}

proxi(d, dx0, iris[,5])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
