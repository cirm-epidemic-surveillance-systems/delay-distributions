#' @title Get the mean of a Normal distribution for a random variable Y
#' needed to ensure that the distribution of X = exp(Y) (which is Log-Normal)
#' has a specified mean and sd.
#' @description
#'  see arithmetic moments here
#' https://en.wikipedia.org/wiki/Log-normal_distribution
#'
#' @param mean target mean for the Log-Normal distribution of X
#' @param sd target sd for the Log-Normal distribution X
#'
#' @return corresponding mean for the underlying Normal
#' distribution of Y = log(X).

#' @export
convert_to_logmean <- function(mean, sd) {
  logmean <- log(mean^2 / sqrt(sd^2 + mean^2))
  return(logmean)
}



#' @title Get the sd of a Normal distribution for a random variable Y
#' needed to ensure that the distribution of X = exp(Y) (which is Log-Normal)
#' has a specified mean and sd.
#' @description see arithmetic moments here
#' https://en.wikipedia.org/wiki/Log-normal_distribution
#'
#' @param mean target mean for the Log-Normal distribution of X
#' @param sd target sd for the Log-Normal distribution of X
#'
#' @return corresponding sd for the underlying Normal distribution of Y = log(X)
#' @export
convert_to_logsd <- function(mean, sd) {
  logsd <- sqrt(log(1 + (sd^2 / mean^2)))
  return(logsd)
}

to_simplex <- function(x) {
  # Handle non-positive values (optional)
  x[x < 0] <- 0
  
  # If all values are 0, return equal proportions
  if(sum(x) == 0) {
    return(rep(1/length(x), length(x)))
  }
  
  # Normalize to sum to 1
  x / sum(x)
}
