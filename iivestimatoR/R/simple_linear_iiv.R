#' Imperfect instrumental variable estimator with one imperfect instrument
#'
#' This estimator calculates the bounds of an IIV with one endogenous regressor
#' 
#' @param data a dataframe that contains all data
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @param cor.direction (string) direction of the correlation 
#' @keywords iiv, simple, one endogenous regressor
#' @return the estimator bounds as well as an indication about the nature of the bounds
#' @export
#' @examples
#' Simple.Linear.IIV()


Simple.Linear.IIV <- function(data,
                              dep.var,
                              endo.regressor,
                              instrument,
                              cor.direction = "xxxx")
{
  OLS <- lm(dep.var ~ endo.regressor, data = data)
  IVZ <- ivreg(dep.var ~ endo.regressor | instrument, data = data)
  OLS.Coef <- unname(OLS$coefficients[2])
  IV.Coef <- unname(IVZ$coefficients[2])
  outcome <- list()
  
  if(cor(endo.regressor, instrument) < 0)
  {
    if(cor.direction == "weak.negative")
    {
      return(paste("Identification Region:","[", IV.Coef, "; ",
                   OLS.Coef, "]", sep = ""))
    }
    if(cor.direction == "weak.positive")
    {
      return(paste("Identification Region:","[", OLS.Coef, "; ", 
                   IV.Coef, "]", sep = ""))
    }
  }
  if( cor(endo.regressor, instrument) >= 0)
  {
    if(cor.direction == "weak.positive")
    {
      Minimum <- min(OLS.Coef, IV.Coef)
      return(paste("One sided bound, coeff. is =< :", 
                   Minimum , sep = " "))
    }
    if(cor.direction == "weak.negative")
    {
      Maximum <- max(OLS.Coef, IV.Coef)
      return(paste("One sided bound, coeff. is >= :", 
                   Maximum , sep = " "))
    }
  }
}