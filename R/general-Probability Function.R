#' Probability Generator
#'
#' Converts a given value into a value in the region (0,1).
#' Uses the exponential cdf. This function is defined
#' separately so that it can be easily changed if needed.
#'
#' @param x some non-negative value
#' @return a value on the range (0,1)
#' @export
probGen <- nimble::nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))
    return(1 - exp(-x))
  }
)
#' Multi Probability Generator
#'
#' Converts n given values into into n + 1 values in the region (0,1).
#' Additionally these values will sum to 1.
#' Uses the exponential cdf to calculate each probability for each weight
#' and uses e^-(sum of weights) to calculate the probability that its not 
#' any of them. The probabilities are then divided by their sum so that they
#' will add to 1.
#'
#' @param x some non-negative values
#' @return a set of values on the range (0,1) that sum to 1
#' @export
multiProbGen <- function(values){
  failValue <- sum(values)
  values <- c(probGen(values),exp(-failValue))
  return(values/sum(values))
}
