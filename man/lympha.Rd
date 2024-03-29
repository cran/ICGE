\name{lympha}
\alias{lympha}
\docType{data}
\title{
Limphatic Database
}
\description{
This lymphography domain was obtained from the University Medical Centre,
   Institute of Oncology, Ljubljana, Yugoslavia.  Thanks go to M. Zwitter and 
   M. Soklic for providing the data. The data are available at the  UCI
   KDD data repository (Hettich S and Bay SD, 1999).
   
   The data set consists of 148 instances presenting 18 different
   mixed attributes (1 cuantitative, 9 binaries and 9 nominals), and a class variable indicating the diagnostic.
There are not missing values.  
}
\usage{data(lympha)}
\format{
 Data frame with 148 instances and 19 features. 
}
\details{Attribute information:
  
    --- NOTE: All attribute values in the database have been entered as
              numeric values corresponding to their index in the list
              of attribute values for that attribute domain as given
        below.
	
    1. class: normal find, metastases, malign lymph, fibrosis
	
    2. lymphatics: normal, arched, deformed, displaced
	
    3. block of affere: no, yes
	
    4. bl. of lymph. c: no, yes
	
    5. bl. of lymph. s: no, yes
	
    6. by pass: no, yes
	
    7. extravasates: no, yes
	
    8. regeneration of: no, yes
	
    9. early uptake in: no, yes
	
   10. lym.nodes dimin: 0-3
	
   11. lym.nodes enlar: 1-4
	
   12. changes in lym.: bean, oval, round
	
   13. defect in node: no, lacunar, lac. marginal, lac. central
	
   14. changes in node: no, lacunar, lac. margin, lac. central
	
   15. changes in stru: no, grainy, drop-like, coarse, diluted, reticular, 
   stripped, faint
   
   16. special forms: no, chalices, vesicles
   
   17. dislocation of: no, yes
   
   18. exclusion of no: no, yes
   
   19. no. of nodes in: 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, >=70
}
\source{
The UCI KDD Archive.
}
\references{
Hettich S and Bay SD (1999). The UCI KDD Archive. Department of
Information and Computer Science. University of California at Irvine,
Irvine, USA.
}
\examples{
data(lympha)
aux <- table(lympha[,1])
barplot(aux, names.arg=c("normal", "metastases", "malign lymph", "fibrosis"))
}
\keyword{datasets}
