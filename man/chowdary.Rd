\name{chowdary}
\alias{chowdary}
\docType{data}
\title{
Chowdary Database
}
\description{
The original authors compared pairs of snap-frozen and RNAlater
preservative-suspended tissue from lymph node-negative breast tumors (B)
and Dukes' B colon tumors (C). The actual data set, by de Souto et. al
(2008), is build with purpose of separating B from C.
}
\usage{data(chowdary)}
\format{
  Data frame with 183 rows and  104 columns.
}
\source{
Original source from `National Center for Biotechnology Information' from the United States of America, query GSE3726.
}
\references{
 de Souto MCP,  Costa IG, de  Araujo DSA, Ludermir TB, and  Schliep A
 (2008). Clustering Cancer Gene Expression Data: a Comparative
 Study. \emph{BMC Bioinformatics}, \bold{8}, 497--511.

  
Chowdary D, Lathrop J, Skelton J, Curtin K, Briggs T, Zhang Y, Yu J,
Wang X, and  Mazumder A (2006). Prognostic gene expression signatures
can be measured in tissues collected in RNAlater
preservative. \emph{Journal Molecular Diagnosis}, \bold{8}, 31--39.
}
\examples{
data(chowdary)

tumor <- as.factor(as.matrix(chowdary[1,]))
x <- as.matrix(chowdary[-1,])
mode(x) <- "numeric"

s <- sample(row.names(x),1)
boxplot( x[s,] ~ tumor , ylab=s)
}
\keyword{datasets}
