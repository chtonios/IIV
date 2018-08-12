#' Beta value calculation for multiple iv regressions
#'
#' This function calculates the criterion as seen in the paper.
#' 
#' @param X.wiggle the regressor matrix with the variance of the exogenous regressors filtered out
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @keywords iiv, beta value, bound
#' @return the criteria needed for multiple iiv regression
#' @export
#' @examples
#' Multiple.Direction()


Multiple.Direction <- function(X.wiggle, endo.regressor, instrument)
{
  sigma.xxwiggle <- cov(X.wiggle, endo.regressor)
  sigma.zxwiggle <- cov(X.wiggle, instrument)
  sigma.z <- sd(instrument)
  sigma.x <- sd(endo.regressor)
  criteria <- (sigma.xxwiggle*sigma.z - sigma.x*sigma.zxwiggle)*sigma.zxwiggle
  return(criteria)
}
