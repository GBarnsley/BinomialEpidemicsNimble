#' A function used internally in the NpmDelta Sampler.
#' Serves to setup values used in the running of the sampler and to provide the
#' dependent nodes of the nodes being sampled.
#' @export
stepSampler_setup <- function(model, mvSaved, target, control) {
  maxStep <- as.integer(control$TMax)
  maxChange <- as.integer(control$DeltaMax)
  #Lists calc nodes for each time point, since change after time t don't effect 
  #changes before time t 
  calcNodesRaw <- list()
  for(i in 1:length(model[[target]])){
    calcNodesRaw[[i]] <- nimbleType(name=as.character(i),type="character", dim = 1)
  }
  calcNodes <- nimbleList(calcNodesRaw)
  calcNodes <- calcNodes$new()
  nodeNames <- c()
  for(i in 1:length(model[[target]])){
    nodeNames[i] <- as.character(i)
    calcNodes[[as.character(i)]] <- model$getDependencies(paste0(target,"[",i,":",length(model[[target]]),"]"))
  }
  #This method requires two posterior calculations for each repeat but involves 
  #less calculations than calculating every node but only once per repeat (+1),
  #up to the point where the number of calculations (nodes) > 10^3 and the number
  #of repeats >10^3
  runs <- as.integer(control$R)
}
#' A function used internally in the NpmDelta Sampler.
#' Runs the NpmDelta algorithm, steps are annotated in the code.
#' @export
stepSampler_run <- function() {
  positions <- which(model[[target]]!=0)
  for(i in 1:runs){
    #Choose the original position
    pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
    position <- positions[pos]
    #Choosing the direction of the move
    if(position == length(model[[target]])){
      #if the position is the last time point we can only go backward
      directions <- c(-1)
      direction <- -1
    }else if(position == 1){
      #if it is the first time point we have to go forward
      directions <- c(1)
      direction <- 1
    }else{
      #in all other cases we choose with uniform probability either direction
      directions <- c(-1,1)
      pos <- rcat(n = 1, prob = c(0.5,0.5))
      direction <- directions[pos]
    }
    #Choosing the number of places to move
    sizes <- min(maxStep, (position - 1)*(direction == -1) +
                   (length(model[[target]]) - position)*(direction == 1))
    #prevents moving further than the time frame of the epidemic
    size <- rcat(n = 1, prob = rep(1/sizes, sizes))
    
    #the index of the position we are moving to
    newPosition <- position + direction*size
    
    #choosing the number of points to move
    amounts <- min(maxChange, model[[target]][position])
    #stops us from moving more than the number of points at that time
    amount <- rcat(n = 1, prob = rep(1, amounts))
    
    sampler_lp_proposed <-
      (-
         #Choosing that time point
         log(length(positions)) -
         #choosing that that amount point of points to move
         log(amounts) -
         #Choosing the size of the move
         log(sizes) -
         #choosing direction
         log(length(directions))
      )
    
    calcIndex <- nodeNames[position]
    if(direction == -1){
      calcIndex <- nodeNames[newPosition]
    }
    
    model_lp_initial <- getLogProb(model, calcNodes[[calcIndex]])
    model[[target]][position] <<- model[[target]][position] - amount
    model[[target]][newPosition] <<- model[[target]][newPosition] + amount
    model_lp_proposed <- calculate(model, calcNodes[[calcIndex]])
    
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
         log(1 + 1*(newPosition != length(model[[target]]) & newPosition != 1))
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
#' A function used internally in the NpmDelta Sampler.
#' Runs the NpmDelta algorithm (bounded version 2, for newI), steps are
#' annotated in the code.
#' @export
stepSampler_run_bounded <- function() {
  pointsToMove <- model[[target]]!=0
  canMoveBackward <- rep(TRUE, length(model[[target]]))
  canMoveBackward[1] <- FALSE
  canMoveForward <- rep(FALSE, length(model[[target]]))
  for(i in 1:(length(model[[target]])-1)){
    canMoveForward[i] <- model[["I"]][i+1] > 1
  }
  positions <- which(pointsToMove * (canMoveBackward + canMoveForward - canMoveBackward*canMoveForward))
  jump <- FALSE
  
  for(i in 1:runs){
    #Choose the original position
    pos <- rcat(n = 1, prob = rep(1, length(positions)))
    position <- positions[pos]
    
    #Choosing the direction of the move
    directions <- canMoveBackward[position] + canMoveForward[position]
    direction <- c(-1,1)[rcat(n=1, prob = c(canMoveBackward[position], canMoveForward[position]))]
    
    #Choosing the number of places to move
    if(direction == 1){
      sizes <- min(maxStep, length(model[[target]]) - position)
      #now we check to make sure we can't propose a step where the amount <= 0
      minAmounts <- rep(0, sizes)
      for(i in 1:sizes){
        minAmounts[i] <- min(model[["I"]][(position + 1):(position + i)]) - 1
      }
      #if there is a value where amount <= 0
      if(min(minAmounts) <= 0){
        sizes <- min(which(minAmounts == 0)) - 1
        #then the step size before this occurs is our new maximum
      }
    }
    else{
      sizes <- min(maxStep, position - 1)
    }
    size <- rcat(n = 1, prob = rep(1, sizes))
    
    #the index of the position we are moving to
    newPosition <- position + direction*size
    #choosing the number of points to move
    amounts <- min(maxChange, model[[target]][position])
    #stops us from moving more than the number of points at that time
    if(direction == 1){
      amounts <- min(amounts, minAmounts[size])
      #also bounded by the minimum amounts for that size, which has already been calculated
    }
    amount <- rcat(n = 1, prob = rep(1, amounts))
    #the probability of generating this proposal
    sampler_lp_proposed <-
      (-
         #Choosing that time point
         log(length(positions)) -
         #choosing that that amount point of points to move
         log(amounts) -
         #Choosing the size of the move
         log(sizes) -
         #choosing direction
         log(directions)
      )
    
    #generating some values for the reverse probability
    #some of these are easier to find when we haven't updated the values
    positionsRev <- length(positions) + (model[[target]][newPosition] == 0) - (amount == model[[target]][position])
    
    
    calcIndex <- nodeNames[position]
    if(direction == -1){
      calcIndex <- nodeNames[newPosition]
    }
    
    
    model_lp_initial <- getLogProb(model, calcNodes[[calcIndex]])
    model[[target]][position] <<- model[[target]][position] - amount
    model[[target]][newPosition] <<- model[[target]][newPosition] + amount
    model_lp_proposed <- calculate(model, calcNodes[[calcIndex]])
    
    #Calculating the rest of the values for the reverse probability
    if(newPosition == length(model[[target]])){
      directionsRev <- 1
      #if we are at the last timepoint we can only go backward
    }
    else{
      directionsRev <-
        (newPosition != 1) +
        #backwards move always possible unless at the start
        (model[["I"]][newPosition + 1] > 1)
      #is a forward move possible
    }
    
    amountsRev <- min(maxChange, model[[target]][newPosition])
    if(direction == -1){
      amountsRev <- min(amountsRev, min(model[["I"]][(newPosition+1):(newPosition + size)]) - 1)
    }
    if(direction == -1){
      sizesRev <- min(maxStep, length(model[[target]]) - newPosition)
      if(sizesRev > size){
        #we know that it must be possible for the size of the proposed jump
        minAmounts <- rep(Inf, sizesRev)
        for(i in (size+1):sizesRev){
          minAmounts[i] <- min(model[["I"]][(newPosition + size + 1):(newPosition + i)]) - 1
        }
        #if there is a value where amount <= 0
        if(min(minAmounts) == 0){
          sizesRev <- min(which(minAmounts == 0)) - 1
          #then the step size before this occurs is our new maximum
        }
      }
    }
    else{
      sizesRev <- min(maxStep, newPosition - 1)
    }
    
    #probability of generating the initial values from the proposed state
    sampler_lp_initial <-
      (-
         #Choosing that time point
         log(positionsRev) -
         #choosing that that amount point of points to move
         log(amountsRev) -
         #choosing the size of the move
         log(sizesRev) -
         #choosing direction
         log(directionsRev)
      )
    log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)
    
    u <- runif(1, 0, 1)
    if(u < exp(log_MH_ratio)){
      model_lp_initial <- model_lp_proposed
      pointsToMove <- model[[target]]!=0
      for(i in 1:(length(model[[target]])-1)){
        canMoveForward[i] <- model[["I"]][i+1] > 1
      }
      positions <- which(pointsToMove * (canMoveBackward + canMoveForward - canMoveBackward*canMoveForward))
      jump <- TRUE
    }
    else{
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
