

#' Title
#'
#' @param Fn.hat.value Estimated CDF from bjr()
#' @param supp.points Points on the support corresponding to Fn.hat.value
#' @param Fb A cumulative density function such as pnorm()
#' @param alpha The "contaminant's" proportion of the mixture distribution
#'
#' @return An object
#' @export
#'
#' @examples
#' 1
patra.sen.bjr <- function(Fn.hat.value, supp.points, Fb, alpha) {
  
  # Modified code based on code provided with
  # "Estimation of a Two-component Mixture Model with Applications to Multiple Testing"
  # https://users.stat.ufl.edu/~rohitpatra/research.html
  
  Fn <- stepfun(supp.points, c(0, Fn.hat.value))
  Fn.1 <- Fn(supp.points) ## Empirical DF of the data at the data points
  pava.freq <- diff(c(0,Fn.1))
  
  F.hat <- (Fn.1 - (1 - alpha) * Fb(supp.points)) / alpha ## Computes the naive estimator of F_s
  Fs.hat <- Iso::pava(F.hat, pava.freq, decreasing = FALSE) ## Computes the Isotonic Estimator of F_s
  Fs.hat[which(Fs.hat <=0 )] <- 0
  Fs.hat[which(Fs.hat >= 1)] <- 1
  
  data.frame(Fs.hat = Fs.hat, F.hat = F.hat)
  
}



