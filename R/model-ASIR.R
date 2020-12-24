#' The expanded SIR epidemic class
#' @export
ASIRclass <- setClass(
  "ASIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
#' Function to create an epidemic model of the expanded SIR class.
#' @param newR Time series that acts as the observed recoveries
#' @param N The population
#' @param t.step The time step of the model
#' @param Frequency TRUE/FALSE value specifying if the model has frequency based transmission
#' @param ChangePoint The time at which we switch Beta values
#' @param TotalInfections The full size of the epidemic, including unobserved infections
#' @return An object of ASIR class with the compiled model code
#' @export
ASIR <- function(newR,
                N,
                t.step = 1,
                Frequency = TRUE,
                ChangePoint){
  tempCode <- nimbleCode({
    # Set priors
    UGamma ~ dgamma(shape = UGammaShape, rate = UGammaRate)
    DGamma ~ dgamma(shape = DGammaShape, rate = DGammaRate)
    Betas[1] ~ T(dnorm(mean = R0Means[1]*(UGamma + DGamma), sd = R0SDs[1]), 0, Inf)
    Betas[2] ~ T(dnorm(mean = R0Means[2]*(UGamma + DGamma), sd = R0SDs[2]), 0, Inf)
    # likelihood
    #S <- integer(length = Timepe)
    S[1] <- Pop - 1
    I[1] <- 1
    for(i in 1:(ChangePoint-1)){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Betas[1]*t.step/(Pop^Frequency)))
      newDR[i] ~ dbinom(size = I[i], prob =  probGen(DGamma*t.step))
      newUR[i] ~ T(dbinom(size = I[i], prob =  probGen(UGamma*t.step)), 0, I[i] - newDR[i])
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newDR[i] - newUR[i]
    }
    for(i in ChangePoint:TimePeriod){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Betas[2]*t.step/(Pop^Frequency)))
      newDR[i] ~ dbinom(size = I[i], prob =  probGen(DGamma*t.step))
      newUR[i] ~ T(dbinom(size = I[i], prob =  probGen(UGamma*t.step)), 0, I[i] - newDR[i])
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newDR[i] - newUR[i]
    }
  })
  return(ASIRclass(
    Model = compileNimble(
      nimbleModel(
        code = tempCode,
        constants = list(TimePeriod = length(newR),
                         ChangePoint = ChangePoint),
        data = list(newDR = newR,
                    t.step = t.step,
                    Pop = N,
                    Frequency = Frequency,
                    #priors
                    UGammaShape = 1,
                    UGammaRate = 1,
                    DGammaShape = 1,
                    DGammaRate = 1,
                    R0Means = rep(1, 2),
                    R0SDs = rep(1, 2)
                    ),
        inits = list(Betas = rep(1, 2),
                     UGamma = 1,
                     DGamma = 1,
                     newI = rep(0, length(newR)),
                     newUR = rep(0, length(newR))
        ),
        calculate = FALSE
      )
    ),
    MCMC = NA,
    Samples = NA
  )
  )
}
#' Initialization method for the expanded SIR model.
#' Sets up provided initial values and runs the NpmDelta algorithm to estimate
#' newI and newUR for those values.
#' @param epiModel An object of the ASIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return ASIR class with the initial values
#' @export
initialValues.ASIR <- function(epiModel, hyperParameters){
  epiModel@Model$Betas <- hyperParameters$`Initial Values`$Betas
  epiModel@Model$UGamma <- hyperParameters$`Initial Values`$UGamma
  epiModel@Model$DGamma <- hyperParameters$`Initial Values`$DGamma
  epiModel@Model$UGammaShape <- hyperParameters$Priors$RecoveryRate$Shape
  epiModel@Model$UGammaRate <- hyperParameters$Priors$RecoveryRate$Rate
  epiModel@Model$DGammaShape <- hyperParameters$Priors$DetectionRate$Shape
  epiModel@Model$DGammaRate <- hyperParameters$Priors$DetectionRate$Rate
  epiModel@Model$R0Means <- hyperParameters$Priors$R0$Means
  epiModel@Model$R0SDs <- hyperParameters$Priors$R0$SDs

  epiModel@Model$newUR <- round(epiModel@Model$newDR*hyperParameters$ProportionUndetected)
  firstRow <- epiModel@Model$newUR[1] + epiModel@Model$newDR[1]
  epiModel@Model$newI <- epiModel@Model$newUR[-1] + epiModel@Model$newDR[-1]
  epiModel@Model$newI[1] <- epiModel@Model$newI[1] + firstRow
  return(
    epiModel
  )
}
#' Method to build an MCMC for the expanded SIR class.
#' Applies a block RWM to the beta and gamma parameters and two NpmDelta algorithms
#' on newI and newUR.
#' @param epiModel An object of the ASIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param showCompilerOutput Whether compileNimble should prince the compiler output
#' @return a complied MCMC
#' @export
buildMCMCInternal.ASIR <- function(epiModel, hyperParameters, showCompilerOutput){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Betas[1]', 'Betas[2]', 'UGamma', 'DGamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addSampler(target = "newUR",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors('Betas')
  output$addMonitors(c('UGamma', 'DGamma'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE, showCompilerOutput=showCompilerOutput)
  return(output)
}
