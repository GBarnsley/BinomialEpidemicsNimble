#' The SEIR epidemic class
#' @export
SEIRclass <- setClass(
  "SEIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
#' Function to create an epidemic model of the SEIR class.
#' @param newR Time series that acts as the recoveries
#' @param N The population
#' @param t.step The time step of the model
#' @param Frequency TRUE/FALSE value specifying if the model has frequency based transmission
#' @return An object of SEIR class with the compiled model code
#' @export
SEIR <- function(newR = NULL,
                 N,
                 t.step = 1,
                 Frequency = TRUE){
  if(!Frequency){
    Freq <- 1
  }else{
    Freq <- 1/N
  }
  inits <- list(Beta = 1,
                Gamma = 1,
                k = 1,
                newI = rep(0, length(newR)),
                newE =  rep(0, length(newR))
  )
  data <- list(newR = newR,
              #priors
              GammaShape = 1,
              GammaRate = 1,
              kShape = 1,
              kRate = 1,
              R0Mean = 1,
              R0SD = 1
  )
  tempCode <- nimbleCode({
    # Set priors
    Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
    Beta ~ T(dnorm(mean = R0Mean[1]*Gamma*Pop*Freq, sd = R0SD[1]), 0, Inf)
    k ~ dgamma(shape = kShape, rate = kRate)
    # likelihood
    S[1] <- Pop - 1
    E[1] <- 0
    I[1] <- 1
    for(i in 1:TimePeriod){
      newE[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Beta*t.step*Freq))
      newI[i] ~ dbinom(size = E[i], prob =  probGen(k*t.step))
      newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
      S[i+1] <- S[i] - newE[i]
      E[i+1] <- E[i] + newE[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newR[i]
    }
  })
  return(SEIRclass(
    Model = compileNimble(
      nimbleModel(
        code = tempCode,
        constants = list(TimePeriod = length(newR),
                         Pop = N,
                         Freq = Freq,
                         t.step = t.step),
        data = data,
        inits = inits,
        calculate = FALSE
      )
    ),
    MCMC = NA,
    Samples = NA
  )
  )
}
#' Method to initialize the SEIR model, sets Beta/Gamma/k by their given initial
#' values then runs the NpmDelta algorithm on newI and newE.
#' @param epiModel An object of the SEIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return SEIR class object with the initial values
#' @export
initialValues.SEIR <- function(epiModel, hyperParameters){
  #inital values
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$k <- hyperParameters$`Initial Values`$k
  #priors
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  epiModel@Model$kShape <- hyperParameters$Priors$k$Shape
  epiModel@Model$kRate <- hyperParameters$Priors$k$Rate
  epiModel@Model$R0Mean <- hyperParameters$Priors$R0$Mean
  epiModel@Model$R0SD <- hyperParameters$Priors$R0$SD
  #generating starting newI, newE values
  first <- sum(epiModel@Model$newR[1:2])
  epiModel@Model$newE <- c(0,epiModel@Model$newR[c(-1,-2)],0)
  epiModel@Model$newE[2] <- epiModel@Model$newE[2] + first
  first <- epiModel@Model$newE[1]
  epiModel@Model$newI <- c(epiModel@Model$newE,0)
  epiModel@Model$newI[1] + epiModel@Model$newI[1] + first - 1
  return(
    epiModel
  )
}
#' Builds the MCMC for the SEIR model. Applies a block-RWM to Beta, Gamma and k, and
#' the NpmDelta sampler to newE and newI.
#' @param epiModel An object of the SEIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param showCompilerOutput Whether compileNimble should prince the compiler output
#' @return a complied MCMC
#' @export
buildMCMCInternal.SEIR <- function(epiModel, hyperParameters, showCompilerOutput){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma', 'k'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newE",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors(c('Beta', 'k', 'Gamma'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE, showCompilerOutput = showCompilerOutput)
  return(output)
}
