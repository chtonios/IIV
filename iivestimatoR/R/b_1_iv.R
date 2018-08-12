#' Beta value calculation for simple IIV
#'
#' This function calculates the beta as seen in the estimation equations
#' 
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @keywords iiv, beta value, bound
#' @return the coefficient value for single iiv estimation
#' @export
#' @examples
#' B.1.IV()


B.1.IV <- function(dep.var, endo.regressor, instrument)
{
  sigma.x <- sd(endo.regressor)
  sigma.z <- sd(instrument)
  sigma.y <- sd(dep.var)
  sigma.zy <- cov(instrument, dep.var)
  sigma.xy <- cov(endo.regressor, dep.var)
  sigma.xz <- cov(endo.regressor, instrument)
  
  beta.value <- ((sigma.x*sigma.zy - sigma.z*sigma.xy)/
                   (sigma.x*sigma.xz - sigma.z*sigma.x*sigma.x))
  return(beta.value)
}
