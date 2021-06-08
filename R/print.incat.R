print.incat <- function(x, ...){
# x must be a list with three components
   cat ("   --INCA test--    \n")
   k <- length(x$ProjectionsU)
   P <- length(x$pvalues)
   if(P > 1)
   {
      cat(" INCA statistic value =", x$StatisticW0,  "\n")
      cat("\n Number of significative tests for alpha= " , x$alpha,": ", 
          format(sum(x$pvalues < x$alpha), digits=3), "\n")
   }else
   {
      cat("\n INCA statistic value =", x$StatisticW0, ", p-value= " , x$pvalues, "\n")
   }

   cat("\n U projections values: \n")
   for (i in 1:k){
        cat("   U_", i," = ", x$ProjectionsU[i], " \n", sep="")
   }
}