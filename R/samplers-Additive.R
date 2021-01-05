#' A function used internally in the Additive Sampler.
#' Serves to setup values used in the running of the sampler and to provide the
#' dependent nodes of the nodes being sampled.
#' @export
additiveSampler_setup <- function(model, mvSaved, target, control) {
  addSubMax <- as.integer(control$addSubMax)
  runs <- control$R
  #getting dependant nodes
  calcNodes <- model$getDependencies(target)
  targetLength <- as.integer(length(model[[target]]))
}
#' A function used internally in the Additive Sampler.
#' Adds or subtracts a value at random and runs the Metropolis hastings algorithm to accept or reject the move
#' @export
additiveSampler_run <- function() {
  model_lp_initial <- getLogProb(model, calcNodes)
  subPositions <- which(model[[target]]!=0)
  for(i in 1:runs){
  
    #Choose to add or subtract
    pos <- rcat(n = 1, prob = c(0.5,0.5))
    addSub <- c(1,-1)[pos]
    
    #Choose the position
    if(addSub == 1){
      #if adding choose any position
      pos <- rcat(n = 1, prob = rep(1/targetLength, targetLength))
      position <- 1:targetLength[pos]
      #determining maximum size
      maxAmount <- addSubMax
    }else{
      #choose a position with values to remove
      pos <- rcat(n = 1, prob = rep(1/length(subPositions), subPositions))
      position <- subPositions[pos]
      #determining maximum size
      maxAmount <- min(model[[target]][position], addSubMax)
    }
    
    #Choosing the size of addition/subtraction
    amount <- rcat(n = 1, prob = rep(1/maxAmount, maxAmount))
    
    sampler_lp_proposed <-
      (-
         #Choosing that time point
         log(if(addSub==1){targetLength}else{length(subPositions)}) -
         #choosing that that amount point of points to move
         log(maxAmount)
      )
    
    model[[target]][position] <<- model[[target]][position] + addSub*amount ###CARRYON FROM HERE
    model_lp_proposed <- calculate(model, calcNodes)
    
    sampler_lp_initial <-
      (-
         #Choosing that time point
         log(sum(model[[target]]!=0)) -
         #choosing that that amount point of points to move
         log(
           min(maxChange, model[[target]][newPosition])
         ) -
         #choosing the size of the move
         log(min(
           maxStep,
           (newPosition - 1)*(-direction == -1) +
             (length(model[[target]]) - newPosition)*(-direction == 1)
         )) -
         #choosing direction
         log(1 + 1*(newPosition != targetLength & newPosition != 1))
      )
    log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)
    
    u <- runif(1, 0, 1)
    if(u < exp(log_MH_ratio)){
      model_lp_initial <- model_lp_proposed
      positions <- which(model[[target]]!=0)
      jump <- TRUE
    }else{
      jump <- FALSE
    }
    ## if we accepted the proposal, then store the updated
    ## values and logProbs from 'model' into 'mvSaved'.
    ## if the proposal was not accepted, restore the values
    ## and logProbs from 'mvSaved' back into 'model'.
    if(jump) copy(from = model, to = mvSaved, row = 1,
                  nodes = calcNodes, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,
              nodes = calcNodes, logProb = TRUE)
  }
}
