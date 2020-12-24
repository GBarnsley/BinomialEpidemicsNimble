#' Generic function that calls the model specific MCMC set-up methods.
#' The naive version of the NpmDelta Algorithm is generated here.
#' @param epiModel An object of the class one of the specific epidemic model
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return An object of the given epidemic class with a compiled MCMC
#' @export
buildMCMCInternal <- function(epiModel, hyperParameters, showCompilerOutput){
  sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = stepSampler_setup,
    run = stepSampler_run,
    methods = list(
      reset = function() {}
    )
  )
  UseMethod("buildMCMCInternal")
}
