\name{estW}
\alias{estW}
\alias{print.incaest}
%\alias{summary.incaest}
\title{ INCA Statistic}
\description{
Assume that n units are divided into k clusters C1,...,Ck, and consider a fixed unit x0. Function \code{estW} calculates the INCA statistic \eqn{W(x0)} and the related \eqn{U_i} statistics.  
}
\usage{
estW(d, dx0, pert = "onegroup")
}
\arguments{
  \item{d}{ a distance matrix or a \code{dist} object with distance information between units.}
  \item{dx0}{ an n-vector containing the distances d0j between x0 and unit j.}
  \item{pert}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are consecutive integers bigger or equal than 1 (for instance 1,2,3,4..., k). The default value indicates the presence of only one group in data. }
}
\value{ The function returns an object of class \code{incaest} which is a list containing the following components:
  \item{Wvalue}{is the INCA statistic \eqn{W(x_0)}.}
  \item{Uvalue }{is a vector containing the statistics \eqn{U_i}. }
}
\references{Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.

 Irigoien, I. and  Arenas, C. (2008). INCA: New statistic for estimating the number of clusters and identifying atypical units.
\emph{Statistics in Medicine}, \bold{27}(15), 2948--2973.
}
\author{Itziar Irigoien \email{itziar.irigoien@ehu.eus}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV/EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\note{ For a correct geometrical interpretation it is convenient to verify whether the distance matrix d is Euclidean.
}
\seealso{ \code{\link{vgeo}}, \code{\link{proxi}} , \code{\link{deltas}} }
\examples{
data(iris)
d <- dist(iris[,1:4])

# characteristics of a specific flower (likely group 1)
x0 <- c(5.3, 3.6, 1.1, 0.1) 
# distances between  flower x0 and the rest of flowers in iris
dx0 <- rep(0,150)
for (i in 1:150){
	dif <-x0-iris[i,1:4]
	dx0[i] <- sqrt(sum(dif*dif))
}
estW(d, dx0, iris[,5])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
