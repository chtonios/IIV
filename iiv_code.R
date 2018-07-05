##---------------------------------------------------##
## Author: Chtonios
## Title: IV Estimation with imperfect instruments
## 2017/2018
## Preliminary version, markdown available
##---------------------------------------------------##

#### Import Packages ####

rm(list = ls())
library(tidyr)
library(AER)
library(readr)
library(MASS)
library(haven)
library(sandwich)

#### Simple Linear IIV ####

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

#### Multiple variables ####

Multiple.Direction <- function(X.wiggle, endo.regressor, instrument)
{
  sigma.xxwiggle <- cov(X.wiggle, endo.regressor)
  sigma.zxwiggle <- cov(X.wiggle, instrument)
  sigma.z <- sd(instrument)
  sigma.x <- sd(endo.regressor)
  criteria <- (sigma.xxwiggle*sigma.z - sigma.x*sigma.zxwiggle)*sigma.zxwiggle
  return(criteria)
}

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

#### Two side bounded estimation ####

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

