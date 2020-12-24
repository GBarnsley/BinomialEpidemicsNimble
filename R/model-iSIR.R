#' Method to initialize the iSIR model, sets Beta/Gamma by their given initial
#' values then runs the NpmDelta algorithm on newI.
#' @param epiModel An object of the iSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return iSIR class object with the initial values
#' @export
initialValues.iSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Shape
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Rate
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  #initial Censored values
  first <- epiModel@Model$newR[1]
  epiModel@Model$newI <- c(epiModel@Model$newR[-1],0)
  epiModel@Model$newI[1] <- epiModel@Model$newI[1] + first
  return(
    epiModel
  )
}
#' Builds the MCMC for the iSIR model. Applies a block-RWM to Beta and Gamma and
#' the NpmDelta sampler to newI
#' @param epiModel An object of the iSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param showCompilerOutput Whether compileNimble should prince the compiler output
#' @return a complied MCMC
#' @export
buildMCMCInternal.iSIR <- function(epiModel, hyperParameters, showCompilerOutput){
  if(!is.null(hyperParameters$newI$Bounded)){
    if(hyperParameters$newI$Bounded){
      sampler <- nimbleFunction(
        contains = sampler_BASE,
        setup = stepSampler_setup,
        run = stepSampler_run_bounded,
        methods = list(
          reset = function() {}
        )
      )
    }}
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors(c('Beta', 'Gamma', 'newI'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE, showCompilerOutput = showCompilerOutput)
  return(output)
}
