# BinomialEpidemicsNimble
A framework for fitting Binomial Epidemic Models with MCMCs in Nimble

An extension of my masters disseration, developing a frame work for applying Bayesian Data Augmention to Binomial Epidemic Models.

Most of this package is wrappers for Nimble code that sets up a model, assigns samplers and then run an MCMC to simulate parameter inference.
The general gist is that you set up a model object (all based on the SIR), normally providing a time-series of Recoveries, and then give this to the MetropolisHastings() function which compiles an MCMC and runs it.
Current Models included are:
  SIR
  iSIR(SIR with missing infection times)
  rSIR(SIR with missing recovery times)
  SEIR(missing infection times)
  ASIR(a more involved SIR, that includes unobserved recoveries and changepoints for infectious weight. the time series provided represents detections)
  COVIDUK(a far more detailed model of the COVID-19 outbreak in the UK)
  
Priors and initial values must be provided in the HyperParameters, alongside the epidemic model, in a list of lists of a certain format, an example for an iSIR is:
hyperParameters <- list(
  `Initial Values` = list(
    Beta = 1,
    Gamma = 0.01,
    Runs = 2000 #some set up generate initial values by running the data augmentation part with fix parameters, this is the number of runs that is used. This may change
  ),
  Priors = list(
    Beta = list( #ignoring the values, these are prior values where Beta and Gamma are from the Gamma Distribution, the notation is the same as for the eqivalent R functions
      Shape = 0.0001,
      Rate = 0.0001
    ),
    Gamma = list(
      Shape = 0.0001*0.01,
      Rate = 0.0001
    )
  ),
  RWM = list( #In some set up a Random Walk-Metropolis is used, this feeds through the standard Nimble options
    propCov = matrix(c(1,0,0,0.01), nrow = 2)
  ),
  `N+-Delta` <- list( #The control options for the time-series sampler, these greatly effect the efficiency of the MCMC
  R = 1,              #The number of times this algorithm runs per update of the population parameters. These occur in sequence
  TMax = 4,           #The maximum number of time-points it will propose moving in a single run
  DeltaMax = 19       #The maximum number of events it will propose moving in a single run
)
)
