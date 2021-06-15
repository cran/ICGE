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
  \item{x}{ a matrix containing, in its rows, the frequencies for each unit. Note: check that each row adds up to 1}
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
#5 individuals represented by their relative frequencies of 4 characteristics (M1-M4):
f <- matrix(c(0.36, 0.21, 0.23, 0.20,
              0.66, 0.18, 0.11, 0.05,
              0.01, 0.24, 0.62, 0.13,
              0.43, 0.38, 0.08, 0.11,
              0.16, 0.07, 0.09, 0.68), 
              byrow=TRUE, nrow=5, dimnames=list(1:5, paste("M", 1:4, sep="")))

# Bhattacharyya distances between pairs 
d <- dbhatta(f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

