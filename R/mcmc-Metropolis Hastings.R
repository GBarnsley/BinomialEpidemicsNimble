#' A function to run the Metropolis-Hastings-within-Gibbs Algorithm, on a
#' given epidemic model. All functions called are generic so this function works
#' for any specified epidemic model.
#' @param epiModel An object of the class one of the specific epidemic model
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @param samples The target number of samples to generate
#' @param burnin The amount of samples to be dropped as burnin
#' @param keep The post-burnin thinning interval
#' @return An object of the same epidemic class, with the produced samples and MCMC in its final state
#' @export
metropolisHastings <- function(epiModel,
                               hyperParameters,
                               samples = 1000,
                               burnin = 500,
                               thin = 10,
                               showCompilerOutput = FALSE){
  epiModel <- initialValues(epiModel, hyperParameters)
  epiModel@MCMC <- buildMCMCInternal(epiModel, hyperParameters, showCompilerOutput)
  if(samples>0){
    epiModel@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    epiModel@Samples <- as.matrix(epiModel@MCMC$mvSamples)
   }
  return(epiModel)
}
