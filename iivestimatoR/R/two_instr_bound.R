#' Imperfect instrumental variable estimation when more than one imperfect instrument is available
#' 
#' This estimator calculates the bounds of the ivv coefficients when more than one imperfect instrument is 
#' available. 
#' 
#' @param data a dataframe that contains all data
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param instrument1 the (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @param instrument2 the second (imperfect) instrument used to estimate the effect of the endogenous regressor
#' @param exog.regressors the regressors that are not correlated with dep.var and the error term
#' @param valid.instruments instruments that can be considered independent of the error term
#' @param gamme.value the weight given to the two instruments (0,1). The closed interval would mean one instrument is not used, therefore it does not make sense to used
#' @keywords iiv, simple, one endogenous regressor
#' @return the beta bounds calculated using two iiv
#' @export
#' @examples
#' Two.Instr.Bounds()

Two.Instr.Bounds <- function(data,
                             endo.regressor,
                             dep.var,
                             exog.regressors,
                             valid.instruments,
                             instrument1,
                             instrument2,
                             gamma.value)
{
  ## creating X and Y wiggle as usual
  IV.X <- ivreg(as.formula(paste("endo.regressor ~ ",
                                 paste(exog.regressors, collapse = "+"),
                                 "|", paste(valid.instruments, collapse = "+"))),
                data = data)
  IV.Y <- ivreg(as.formula(paste("dep.var ~ ",
                                 paste(exog.regressors, collapse = "+"),
                                 "|", paste(valid.instruments, collapse = "+"))),
                data = data)
  X.wiggle <- as.matrix(IV.X$residuals)
  Y.wiggle <- as.matrix(IV.Y$residuals)
  
  ## creating variables given inputs
  ## see the term paper or Nevo/Rosens paper for details
  omega <- ((gamma.value*instrument2) - ((1-gamma.value)*instrument1))
  Nu1 <- ((sd(endo.regressor)*instrument1) - (sd(instrument1)*endo.regressor))
  Nu2 <- ((sd(endo.regressor)*instrument2) - (sd(instrument2)*endo.regressor))
  Nu.star <- ((sd(endo.regressor)*omega) - (sd(omega)*endo.regressor))
  
  ## Creating matrix needed for inference later on
  inferencematrix <- matrix(nrow = 5, ncol = 2)
  
  ## IV regressions
  IVNu1 <- ivreg(Y.wiggle ~ X.wiggle | Nu1)
  IVNu2 <- ivreg(Y.wiggle ~ X.wiggle | Nu2)
  IVZ1 <- ivreg(Y.wiggle ~ X.wiggle |  instrument1, data=data)
  IVZ2 <- ivreg(Y.wiggle ~ X.wiggle |  instrument2, data=data)
  IVOmeg <- ivreg(Y.wiggle ~ X.wiggle | omega)
  IVNu.star <- ivreg(Y.wiggle ~ X.wiggle |  Nu.star)
  
  ## Isolating coefficients and std.errs
  ## needed to construct bounds and inference
  B.z1 <- unname(IVZ1$coefficients[2])
  B.z2 <- unname(IVZ2$coefficients[2])
  B.Nu.star <- unname(IVNu.star$coefficients[2])
  Nu1 <- unname(IVNu1$coefficients[2])
  Nu2 <- unname(IVNu2$coefficients[2])
  B.Omeg <- unname(IVOmeg$coefficients[2])
  
  ## pack everything neatly into the matrix
  ## coefs
  inferencematrix[1,1] <- B.z1
  inferencematrix[2,1] <- B.z2
  inferencematrix[3,1] <- B.Nu.star
  inferencematrix[4,1] <- Nu1
  inferencematrix[5,1] <- Nu2
  ## std.errs
  inferencematrix[1,2] <- sqrt(vcovHC(IVZ1)[2,2])
  inferencematrix[2,2] <- sqrt(vcovHC(IVZ2)[2,2])
  inferencematrix[3,2] <- sqrt(vcovHC(IVNu.star)[2,2])
  inferencematrix[4,2] <- sqrt(vcovHC(IVNu1)[2,2])
  inferencematrix[5,2] <- sqrt(vcovHC(IVNu2)[2,2])
  rownames(inferencematrix) <- c("Instr.1", "Instr.2", "Nu.star", "Nu1", "Nu2")
  colnames(inferencematrix) <- c("Estimate", "Std. Error")
  
  ## returning the values
  Bounds <- matrix(nrow=1, ncol=2)
  Bounds[1,1] <- B.Omeg
  Bounds[1,2] <- min(Nu1, Nu2, B.z1, B.z2, B.Nu.star)
  colnames(Bounds) <- c("Lower Bound", "Upper bound")
  
  list(coefs = Bounds , matrix = inferencematrix, B.Omeg = B.Omeg, IVOmeg = IVOmeg)
  
}