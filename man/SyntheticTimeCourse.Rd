\name{SyntheticTimeCourse}
\alias{SyntheticTimeCourse}
\docType{data}
\title{
Synthetic Time Course data
}
\description{
Sythetic time course data  where 210 genes profiles along 6 time points
are reported and where the genes are drawn from 8 different populations.
}
\usage{data(SyntheticTimeCourse)}
\format{
  Data frame with 120 rows and  7 columns.
}
\details{
  Attribute information:
  Column cl: the class that the gen belongs to.
  Columns t1 - t6: gene's expression along the t1, ..., t6 time
  points considered.
}

\examples{
data(SyntheticTimeCourse)
x <- SyntheticTimeCourse[, 2:7]
cl <- SyntheticTimeCourse[, 1]
par(mfrow=c(3,3))
for (g in 1:8){ 
   xx <- t(x[cl==g,] )
   yy <- matrix(c(1:6 ), nrow=6, ncol=15, byrow=FALSE)
   matplot(yy,xx,  pch=21, type="b", axes=FALSE,  
        ylim=c(0,3.5), xlim=c(0.5,6.5), xlab="", ylab="", col="black", main=paste("G",g)) 
   abline(h=0)  
   abline(v=0.5) 
   mtext("Time", side=1) 
   mtext("Expression", side=2)  
}
}
\keyword{datasets}
