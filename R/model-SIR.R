#' The SIR epidemic class
#' @export
SIRclass <- setClass(
  "SIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
#' The iSIR epidemic class
#' @export
newISIRclass <- setClass(
  "iSIR",
  contains = "SIR"
)
#' The rSIR epidemic class
#' @export
newRSIRclass <- setClass(
  "rSIR",
  contains = "SIR"
)
#' Creates an object of class SIR, iSIR or rSIR
#'
#' Generates and compiles SIR epidemic model.
#'
#' @param newI Time series of new infections
#' @param newR Time series of new recoveries
#' @param N The size of the population
#' @param Frequency True/False value indicating if the model has frequency based tranmission
#' @param t.step The amount of time in each step of the model, actual value depends on the units of Beta and Gamma
#' @return Object of class SIR, iSIR or rSIR which will contain the compiled nimble model
#' @export
SIR <- function(newI = NULL,
                newR = NULL,
                N = NULL,
                t.step = 1,
                Frequency = TRUE){
  #calculating initial values from given dataS
  if(is.null(N)){
    print("Error: N must be specified")
    return(NA)
  }
  #Setting up models based on nulls
  tempCode <- nimbleCode({
    # Set priors
    Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
    Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
    # likelihood
    S[1] <- Pop - 1
    I[1] <- 1
    for(i in 1:TimePeriod){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Beta*t.step/(Pop^Frequency)))
      newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newR[i]
    }
  })
  if(is.null(newI)){
    return(newISIRclass(
      Model = compileNimble(
        nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newR)),
          data = list(newR = newR,
                      t.step = t.step,
                      Pop = N,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1,
                      Frequency = Frequency),
          inits = list(Beta = 1,
                       Gamma = 1,
                       newI = c(sum(newR) - 1, rep(0, length(newR) - 1))
          ),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA
    )
    )
  }
  else if(is.null(newR)){
    return(newRSIRclass(
      Model = compileNimble(
        nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newI)),
          data = list(newI = newI,
                      Pop = N,
                      t.step = t.step,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1,
                      Frequency = Frequency),
          inits = list(Beta = 1,
                       Gamma = 1,
                       newR = c(rep(0, length(newI) - 1), sum(newI) + 1)
          ),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA
    )
    )
  }
  else{
    tempCode <- nimbleCode({
      # Set priors
      Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
      Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
      # likelihood
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step/(Pop^Frequency)))
        newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
      }
    })
    return(SIRclass(
      Model = compileNimble(
        nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newI)),
          data = list(I = I,
                      S = S,
                      newI = newI,
                      newR = newR,
                      t.step = t.step,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1,
                      Frequency = Frequency,
                      Pop = N),
          inits = list(Beta = 1,
                       Gamma = 1
          ),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA
    )
    )
  }
}
#' Builds the MCMC for the SIR model. Applies a block-RWM to Beta and Gamma.
#' @param epiModel An object of the SIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return a complied MCMC
#' @export
buildMCMCInternal.SIR <- function(epiModel, hyperParameters, showCompilerOutput){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  return(compileNimble(
    buildMCMC(
      output
    ),
    project = epiModel@Model,
    resetFunctions = TRUE,
    showCompilerOutput = showCompilerOutput
  ))
}
#' Method to initialize the SIR model, sets Beta/Gamma by their given initial
#' values.
#' @param epiModel An object of the SIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param showCompilerOutput Whether compileNimble should prince the compiler output
#' @return SIR class object with the initial values
#' @export
initialValues.SIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Shape
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Rate
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  return(
    epiModel
  )
}
