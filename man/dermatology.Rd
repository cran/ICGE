\name{dermatology}
\alias{dermatology}
\docType{data}
\title{  Dermatology Database }
\description{
   Data from a dermatology study provided by H.A. Guvenir (Dpt. Computer Engineering and Information Science, Bilkent University, Turkey).The data set contains 366 instances presenting 34 different clinical attributes (12 clinical features as age or family history and 22 histopathological features obtained from a biopsy), and a class variable indicating the disease.
There are 8 missing values. 
This data set has been used extensively for classification tasks. 
}
\usage{data(dermatology)}
\details{Attribute information obtained from the  UCI KDD data repository:
  
       Clinical Attributes: (they take values 0, 1, 2, 3, unless otherwise
       indicated)
       
      1: erythema;
      2: scaling;
      3: definite borders;
      4: itching;
      5: koebner phenomenon;
      6: polygonal papules;
      7: follicular papules;
      8: oral mucosal involvement;
      9: knee and elbow involvement;
     10: scalp involvement;
     11: family history, (0 or 1);
     34: Age.

     Histopathological Attributes: (they take values 0, 1, 2, 3)
     
     12: melanin incontinence;
     13: eosinophils in the infiltrate;
     14: PNL infiltrate;
     15: fibrosis of the papillary dermis;
     16: exocytosis;
     17: acanthosis;
     18: hyperkeratosis;
     19: parakeratosis;
     20: clubbing of the rete ridges;
     21: elongation of the rete ridges;
     22: thinning of the suprapapillary epidermis;
     23: spongiform pustule;
     24: munro microabcess;
     25: focal hypergranulosis;
     26: disappearance of the granular layer;
     27: vacuolisation and damage of basal layer;
     28: spongiosis;
     29: saw-tooth appearance of retes;
     30: follicular horn plug;
     31: perifollicular parakeratosis;
     32: inflammatory monoluclear inflitrate;
     33: band-like infiltrate.

     The considered diseases are: 1 - psoriasis, 2 - seboreic dermatitis, 3- lichen planus, 4 - pityriasis rosea, 5 - chronic dermatitis, 6 - pityriasis rubra pilaris. 
}
\format{
  Matrix with 366 rows.
}
\source{
  The UCI KDD Archive.
}
\references{
  Guvenir H, Demiroz G, Ilter N (1998). Learning differential
diagnosis of erythemato-squamous diseases using voting feature
intervals. \emph{Artificial Intelligence in Medicine}, \bold{13},
147--165.

Irigoien I, Arenas C (2008). INCA: New statistic for estimating the
number of clusters and identifying atypical units. \emph{Statistics
in Medicine}, \bold{27}, 2948--2973.
}
\examples{
data(dermatology)
x <- dermatology[, 1:34]
group <- as.factor(dermatology[,35])

plot(group)
}
\keyword{datasets}
