#' Method to initialize the rSIR model, sets Beta/Gamma by their given initial
#' values then runs the NpmDelta algorithm on newR.
#' @param epiModel An object of the rSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return rSIR class object with the initial values
#' @export
initialValues.rSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Shape
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Rate
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  epiModel@Model$newR <- rep(0, length(epiModel@Model$newI))
  epiModel@Model$newR[length(epiModel@Model$newI)] <- sum(epiModel@Model$newI) + 1
  #inital censored values
  last <- epiModel@Model$newI[length(epiModel@Model$newI)]
  epiModel@Model$newR <- c(0,epiModel@Model$newI[-length(epiModel@Model$newI)])
  epiModel@Model$newR[length(epiModel@Model$newR)] <- epiModel@Model$newR[length(epiModel@Model$newR)] + last
  return(
    epiModel
  )
}
#' Builds the MCMC for the rSIR model. Applies a block-RWM to Beta and Gamma and
#' the NpmDelta sampler to newR.
#' @param epiModel An object of the rSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param showCompilerOutput Whether compileNimble should prince the compiler output
#' @return a complied MCMC
#' @export
buildMCMCInternal.rSIR <- function(epiModel, hyperParameters, showCompilerOutput){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newR",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors(c('Beta', 'Gamma', 'newR'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE, showCompilerOutput = showCompilerOutput)
  return(output)
}
