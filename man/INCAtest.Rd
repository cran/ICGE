\name{INCAtest}
\alias{INCAtest}
\alias{print.incat}
\alias{summary.incat}
\title{ INCA Test}
\description{
 Assume that n units are divided into k groups C1,...,Ck. Function \code{INCAtest} performs the typicality INCA test. Therein, the null hypothesis that a new unit x0 is a typical unit with respect to a previously fixed partition is tested versus the alternative hypothesis that the unit is atypical.
}
\usage{
INCAtest(d, pert, d_test, np = 1000, alpha = 0.05, P = 1)
}
\arguments{
  \item{d}{a distance matrix or a \code{dist} object with distance information between units.}
  \item{pert}{ an n-vector that indicates which group each unit belongs to. Note that the expected values of \code{pert} are numbers greater than or equal to 1 (for instance 1,2,3,4..., k). The default value indicates there is only one group in data.}
  \item{d_test}{ an n-vector containing the distances from x0 to the other units.}
  \item{np}{ sample size for the bootstrap sample for the bootstrap procedure.}
  \item{alpha}{ fixed level for the test.}
  \item{P}{ Number of times the bootstrap procedure is repeated.}
}
\value{A list with class "incat" containing the following components:
  \item{StatisticW0 }{value of the INCA statistic.}
  \item{ProjectionsU }{values of statistics measuring the projection from the specific object to each considered group.}
  \item{pvalues}{p-values obtained in the \code{P} times repeated bootstrap procedure. Note: If \code{P}>1, it is printed the number of times the p-values were smaller than \code{alpha}.}
  \item{alpha }{specified value of the level of the test.}
}
\references{ Irigoien, I. and  Arenas, C. (2008). INCA: New statistic for estimating the number of clusters and identifying atypical units.
\emph{Statistics in Medicine}, \bold{27}(15), 2948--2973.


Arenas, C. and Cuadras, C.M. (2002). Some recent statistical methods based on distances.\emph{ Contributions to Science}, \bold{2}, 183--191.
 }
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.es}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\note{ 
To obtain the INCA statistic distribution, under the null hypothesis, the program can consume long time.
For a correct geometrical interpretation it is convenient to verify whether the distance matrix d is Euclidean.
}
\seealso{\code{\link{estW}}, \code{\link{INCAindex}} }
\examples{#generate 3 clusters, each of them with 20 objects in dimension 5.
mu1 <- sample(1:10, 5, replace=TRUE)
x1 <- matrix(rnorm(20*5, mean = mu1, sd = 1),ncol=5, byrow=TRUE)
mu2 <- sample(1:10, 5, replace=TRUE)
x2 <- matrix(rnorm(20*5, mean = mu2, sd = 1),ncol=5, byrow=TRUE)
mu3 <- sample(1:10, 5, replace=TRUE)
x3 <- matrix(rnorm(20*5, mean = mu3, sd = 1),ncol=5, byrow=TRUE)
x <- rbind(x1,x2,x3)

# Euclidean distance between units in matrix x.
d <- dist(x)
# given the right partition
partition <- c(rep(1,20), rep(2,20), rep(3,20))

# x0 contains a unit from one group, as for example group 1.
x0 <-  matrix(rnorm(1*5, mean = mu1, sd = 1),ncol=5, byrow=TRUE)

# distances between x0 and the other units.
dx0 <- rep(0,60)
for (i in 1:60){
	dif <-x0-x[i,]
	dx0[i] <- sqrt(sum(dif*dif))
}

INCAtest(d, partition, dx0, np=10)


# x0 contains a unit from a new group.
x0 <-  matrix(rnorm(1*5, mean = sample(1:10, 5, replace=TRUE),
        sd = 1), ncol=5, byrow=TRUE)

# distances between x0 and the other units in matrix x.
dx0 <- rep(0,60)
for (i in 1:60){
	dif <-x0-x[i,]
	dx0[i] <- sqrt(sum(dif*dif))
}

INCAtest(d, partition, dx0, np=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ cluster}
