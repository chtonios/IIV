#' Imperfect instrumental variable estimator with one imperfect instrument. Additionally the assumption 4
#' from the paper is imposed (the instrument is less endogenous than the endogenous regressor)
#'
#' This estimator calculates the bounds of an IIV with one endogenous regressor, with A4 as well
#'
#' @param data a dataframe that contains all data
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @param cor.direction (string) direction of the correlation
#' @keywords iiv, simple, one endogenous regressor, instrument less endogenous
#' @return the estimator bounds as well as an indication about the nature of the bounds under a4
#' @export
#' @examples
#' Simple.Linear.IIV.A4()


Simple.Linear.IIV.A4 <- function(data,
                                 dep.var,
                                 endo.regressor,
                                 instrument,
                                 cor.direction = "xxxx")
{
  Beta.One <- B.1.IV(dep.var, endo.regressor, instrument)
  IVZ <- ivreg(dep.var ~ endo.regressor | instrument, data = data)
  IV.Coef <- unname(IVZ$coefficients[2])

  if(cor(endo.regressor, instrument) < 0)
  {
    if(cor.direction == "weak.positive")
    {
      return(paste("Identification Region:","[", IV.Coef,
                   "; ", Beta.One, "]", sep = ""))
    }
    if(cor.direction == "weak.negative")
    {
      return(paste("Identification Region:","[", Beta.One,
                   "; ", IV.Coef, "]", sep = ""))
    }
  }
  if( cor(endo.regressor, instrument) >= 0)
  {
    if(cor.direction == "weak.positive")
    {
      Minimum <- min(Beta.One, IV.Coef)
      return(paste("One sided bound, coeff. is =<:",
                   Minimum , sep = " "))
    }
    if(cor.direction == "weak.negative")
    {
      Maximum <- max(Beta.One, IV.Coef)
      return(paste("One sided bound, coeff. is >= :",
                   Maximum , sep = " "))
    }
  }
}
