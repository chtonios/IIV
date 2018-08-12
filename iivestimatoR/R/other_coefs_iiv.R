#' Other coefficients for a multiple imperfect instrumental variable regression 
#'
#' This estimator calculates the bounds and standard errors for the other, non endogenous variables in the 
#' iiv model. 
#' 
#' @param data a dataframe that contains all data
#' @param dep.var the dependend variable (y)
#' @param endo.regressor the regressor that is correlated with the error term(x)
#' @param exog.regressors the regressors that are not correlated with dep.var and the error term
#' @param valid.instruments instruments that can be considered independent of the error term
#' @param beta.low the lower bound of the estimation for the endogenous variable
#' @param beta.upp the upper bound of the estimation for the endogenous variable
#' @keywords iiv, simple, one endogenous regressor
#' @return estimators and std. errors for the exogenous variables in the iiv equation
#' @export
#' @examples
#' Multiple.Linear.IIV.other()


Multiple.Linear.IIV.other <- function(data,
                                      dep.var,
                                      endo.regressor,
                                      exog.regressors = c(),
                                      valid.instruments = c(),
                                      beta.low,
                                      beta.upp)
{
  ## Netting out the effects of all other instruments to get 
  ## the W.wiggle matrix
  W.wiggle <- matrix(nrow = dim(data)[1], ncol = length(exog.regressors))
  for(j in 1:length(exog.regressors))
  {
    IV <- ivreg(as.formula(paste(exog.regressors[j],"~",
                                 paste(exog.regressors[-j],collapse = "+"),
                                 "|", paste(valid.instruments[-j], collapse = "+"))),
                data = data)
    W.wiggle[,j] <- as.matrix(IV$residuals)
    
  }
  ##Add the constant term
  constant <- matrix(rep(1, times= dim(data)[1]), ncol = 1)
  IV.Const <- ivreg(as.formula(paste("constant ~ ","0","+", paste(exog.regressors,
                                                                  collapse = "+"),
                                     "|","0","+", paste(valid.instruments, collapse = "+"))),
                    data=data)
  Const.Wiggle <- as.matrix(IV.Const$residuals)
  W.matrix <- cbind(Const.Wiggle,W.wiggle)
  
  ## Getting X and Y wiggle the same way
  ## Note that we need to calculate a Y.wiggle for 
  ## each exogenous regressor and cannot use a single
  ## Y wiggle as it is suggested in one equation in the paper.
  ## The proof in the appendix correctly states the need
  ## to use a seperate estimate for Y wiggle for each 
  ## delta coefficient
  
  Y.wiggle <- matrix(nrow = dim(data)[1], ncol = length(exog.regressors))
  for(t in 1:length(exog.regressors))
  {
    IV <- ivreg(as.formula(paste("dep.var","~",
                                 paste(exog.regressors[-t],collapse = "+"),
                                 "|", paste(valid.instruments[-t], collapse = "+"))),
                data = data)
    
    
    Y.wiggle[,t] <- IV$residuals 
    
  }
  ## Add constant
  IV.ConstY <- ivreg(as.formula(paste("dep.var","~","0","+",
                                      paste(exog.regressors,collapse = "+"),
                                      "|", "0","+", paste(valid.instruments, collapse = "+"))),
                     data = data)
  
  Const.YWiggle <- as.matrix(IV.ConstY$residuals)
  Y.matrix <- cbind(Const.YWiggle,Y.wiggle)
  
  
  X.wiggle <- matrix(nrow = dim(data)[1], ncol = length(exog.regressors))
  for(w in 1:length(exog.regressors))
  {
    IV <- ivreg(as.formula(paste("endo.regressor","~",
                                 paste(exog.regressors[-w],collapse = "+"),
                                 "|", paste(valid.instruments[-w], collapse = "+"))),
                data = data)
    
    
    X.wiggle[,w] <- as.matrix(IV$residuals) 
    
  }
  ## Add constant
  IV.ConstX <- ivreg(as.formula(paste("endo.regressor","~","0","+",
                                      paste(exog.regressors,collapse = "+"),
                                      "|", "0","+", paste(valid.instruments, collapse = "+"))),
                     data = data)
  
  Const.XWiggle <- as.matrix(IV.ConstX$residuals)
  X.matrix <- cbind(Const.XWiggle,X.wiggle)
  
  ## Calculate the delta coefficients (the coefs for all
  ## exogenous or consistenly estimatable coefs)
  Coefmatrix <- matrix(nrow=2, ncol = (length(valid.instruments)))
  Sdmatrix <- matrix(nrow=2, ncol = (length(valid.instruments)))
  for(k in 2:dim(Y.matrix)[2])
  {
    Delta.IV.Y <- ivreg(Y.matrix[,k] ~ 0 + W.matrix[,k] | 0 + as.matrix(unname(data[valid.instruments[(k-1)]])))
    Delta.IV.X <- ivreg(X.matrix[,k] ~ 0 + W.matrix[,k] | 0 + as.matrix(unname(data[valid.instruments[(k-1)]])))
    
    Lower <- unname(Delta.IV.Y$coefficients - Delta.IV.X$coefficients*beta.low)
    Upper <- unname(Delta.IV.Y$coefficients - Delta.IV.X$coefficients*beta.upp)
    
    Ysd <- sqrt(vcov(Delta.IV.Y))
    Xsd <- sqrt(vcov(Delta.IV.X))
    
    Coefmatrix[1,(k-1)] <- Lower
    Coefmatrix[2,(k-1)] <- Upper
    
    Sdmatrix[1,(k-1)] <- Ysd
    Sdmatrix[2,(k-1)] <- Xsd
  }
  ## seperate regression for the constant term
  Constmat <- matrix(nrow = 2, ncol = 1)
  Constmatsd <- matrix(nrow = 2, ncol = 1)
  
  Delta.IV.Y.Const <- ivreg(Y.matrix[,1] ~ 0 + W.matrix[,1] | 0 + constant)
  Delta.IV.X.Const <- ivreg(X.matrix[,1] ~ 0 + W.matrix[,1] | 0 + constant)
  
  Lower <- unname(Delta.IV.Y.Const$coefficients - Delta.IV.X.Const$coefficients*beta.low)
  Upper <- unname(Delta.IV.Y.Const$coefficients - Delta.IV.X.Const$coefficients*beta.upp)
  
  Ysd <- sqrt(vcov(Delta.IV.Y.Const))
  Xsd <- sqrt(vcov(Delta.IV.X.Const))
  
  Constmat[1,1] <- Lower
  Constmat[2,1] <- Upper
  
  Constmatsd[1,1] <- Ysd
  Constmatsd[2,1] <- Xsd
  
  ## combine and add names to the output matrix 
  Coefs <- cbind(Constmat, Coefmatrix)
  Stds <- cbind(Constmatsd, Sdmatrix)
  colnames(Coefs) <- c("Constant",paste(exog.regressors))
  rownames(Coefs) <- c("Bound A", "Bound B")
  rownames(Stds) <- c("YStds", "XStds")
  
  list(Coefs = Coefs, Stderrs = Stds)
}