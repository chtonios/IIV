#' Imperfect instrumental variable with one or more imperfect instruments. 
#'
#' This estimator calculates the bounds of the ivv coefficients when one ore more regressors are exogenous
#' and/or can be estimated with valid instruments
#' 
#' @param data a dataframe that contains all data
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @param exog.regressors the regressors that are not correlated with dep.var and the error term
#' @param valid.instruments instruments that can be considered independent of the error term
#' @param covar.direction (string) direction of the covariance 
#' @keywords iiv, simple, one endogenous regressor
#' @return beta estimator as well as indication if it is a one or two-sided bound
#' @export
#' @examples
#' Multiple.Linear.IIV.beta()


Multiple.Linear.IIV.beta <- function(data,
                                     dep.var,
                                     endo.regressor,
                                     instrument,
                                     exog.regressors = c(),
                                     valid.instruments = c(),
                                     covar.direction = "xxxx")
{
  #### Netting out effects of exogenous regressors ####
  IV.Reg1 <- ivreg(as.formula(paste("endo.regressor ~ ",
                                    paste(exog.regressors,collapse = "+"),
                                    "|", paste(valid.instruments, collapse = "+"))),
                   data = data)
  IV.Reg2 <- ivreg(as.formula(paste("dep.var ~ ",
                                    paste(exog.regressors,collapse = "+"),
                                    "|", paste(valid.instruments, collapse = "+"))),
                   data = data)
  
  #### Calculating X and Y wiggle ####
  X.coef <- matrix(unname(IV.Reg1$coefficients),
                   nrow = length(unname(IV.Reg2$coefficients)), ncol = 1)
  Y.coef <- matrix(unname(IV.Reg2$coefficients),
                   nrow = length(unname(IV.Reg2$coefficients)), ncol = 1)
  
  Exog.wo.const <- as.matrix(unname(data[,exog.regressors]))
  constant <- rep(1, dim(Exog.wo.const)[1])
  W.matrix <- cbind(constant, Exog.wo.const)
  
  X.wiggle <- endo.regressor - W.matrix %*% X.coef
  Y.wiggle <- dep.var - W.matrix %*% Y.coef
  
  #### Getting the IV bounds ####
  IV.Reg3 <- ivreg(Y.wiggle ~  0 + X.wiggle | 0 + instrument, data = data)
  IV.Coef <- unname(IV.Reg3$coefficients[1])
  NuO1 <- ((sd(endo.regressor)*instrument)-(sd(instrument)*endo.regressor))
  Addreg <- ivreg(Y.wiggle ~ 0 + X.wiggle | 0 + NuO1, data = data)
  Beta.One <- unname(Addreg$coefficients[1])
  
  #### Applying the results ####
  outcome <- matrix(ncol = 2, nrow=1)
  if(Multiple.Direction(X.wiggle, endo.regressor, instrument) < 0)
  {
    if(covar.direction == "negative")
    {
      colnames(outcome) <- c("Two sided identification region",":")
      outcome[1,1] <- IV.Coef
      outcome[1,2] <- Beta.One
      return(outcome)
    }
    if(covar.direction == "weak.positive")
    {
      colnames(outcome) <- c("Two sided identification region",":")
      outcome[1,1] <- Beta.One
      outcome[1,2] <- IV.Coef
      return(outcome)
    }
  }
  if(Multiple.Direction(X.wiggle, endo.regressor, instrument) > 0)
  {
    if(covar.direction == "negative")
    {
      Maximum <- max(Beta.One, IV.Coef)
      colnames(outcome) <- c("One sided bound, Beta Coefficient is",">=")
      outcome[1,1] <- Maximum
      return(outcome)
    }
    if(covar.direction == "weak.positive")
    {
      Minimum <- min(Beta.One, IV.Coef)
      outcome[[1]] <- c("One sided bound, Beta Coefficient is" ,"<=")
      outcome[1,1] <- Minimum
      return(outcome)
    }
  }
}