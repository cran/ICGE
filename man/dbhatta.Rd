\name{dbhatta}
\alias{dbhatta}
\title{ Bhattacharyya Distance }
\description{
  \code{dbhatta}  computes and returns the Bhattacharyya distance matrix between the rows of a data matrix. This distance is defined between two units \eqn{i=(p_{i1},...,p_{im})} and \eqn{j=(p_{j1},...,p_{jm})} being \eqn{p_{kl}} frequencies with \eqn{p_{kl}>=0} and \eqn{p_{k1}+...+p_{km}=1}.
}
\usage{
dbhatta(x)
}
\arguments{
  \item{x}{ a matrix containing, in its rows, the frequencies for each unit.}
}

\value{
  A \code{\link{dist}} object with distance information.
}
\references{ Bhattacharyya, A. (1946). On a measure of divergence of two multinomial populations.
  \emph{Sankhya: The Indian Journal of Statistics, Series A}. \bold{14}, 177-136.}
\author{  Itziar Irigoien \email{itziar.irigoien@ehu.eus}; Konputazio Zientziak eta Adimen Artifiziala, Euskal Herriko Unibertsitatea (UPV-EHU), Donostia, Spain.

Conchita Arenas \email{carenas@ub.edu}; Departament d'Estadistica, Universitat de Barcelona, Barcelona, Spain.
}
\seealso{  \code{\link{dist}}, \code{\link{dmahal}},
  \code{\link{dgower}}, \code{\link{dcor}}, \code{\link{dproc2}} }
\examples{
#Generate 10 objects in dimension 4. Number of trials M=30. 
#Frequencies obtained using tabulate function.
n <- 10
f <- matrix(0, n, 4)
M <- 30
for (i in 1:n){
    f[i,] <- tabulate(sample(1:4, M, replace=TRUE))/M
}


# Bhattacharyya distances between pairs 
d <- dbhatta(f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

