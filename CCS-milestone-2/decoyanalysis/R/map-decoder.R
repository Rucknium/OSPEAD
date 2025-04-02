

#' MAP Decoder Effectiveness
#'
#' @param f_S Probability mass function (PMF) of the real spend distribution.
#' This should be a vector of probabilities summing to one.
#' @param f_D PMF of the decoy distribution. Must be the same length as `f_S`.
#' @param n.decoys Number of decoys.
#' @param subsample A number between 0 and 1. Takes a subsample of the support
#' of `f_S`, as a proportion, to speed up computation. The subsample is generated
#' by weighted random sampling without replacement, using the implementation of
#' `wrswoR::sample_int_expj()`. Recommended not to go below 0.1. Alternatively,
#' a vector of support points. This alternative specification can be used
#' if the user desires to specify the support subsample directly.
#' @param rounding.digit Optional. If specified, round an intermediate
#' vector to the nearest digit. Helps speed up computation, with minimal
#' loss of accuracy of the final result. Setting to 5 seems to work well.
#' @param return.only.success.probability.vector If TRUE, returns only the
#' `success.probability.vector` object
#' 
#' @details
#' Computes equation 13 of Aeeneh et al. (2021).
#' 
#' Let \eqn{Z} be the total number of outputs eligible to be spent in a RingCT
#' ring. Let \eqn{f_S} be the real spend probability mass function (PMF). Let
#' \eqn{f_D} be a proposed decoy PMF. Let \eqn{M} be the number of decoys in a
#' ring. Let \eqn{\mathbf{1}\{\}} be the indicator function that takes value 1
#' when the statement within the braces is true and 0 otherwise. The average
#' success probability of the MAP Decoder attack is
#' \deqn{L_{MAP\:Decoder}=\sum_{i=1}^{Z}f_{S}\left(i\right)\left(\sum_{j=1}^{Z}f_{D}\left(j\right)\mathbf{1}\left\{ \tfrac{f_{S}\left(j\right)}{f_{D}\left(j\right)}<\tfrac{f_{S}\left(i\right)}{f_{D}\left(i\right)}\right\} \right)^{M}}
#' 
#'
#' @return A list with the following elements:
#' - `mean.success.probability`: A scalar, between 0 and 1. The value of
#' \eqn{L_{MAP\:Decoder}}.
#' - `success.probability.vector`: A vector of the same length as `sub.supp`.
#' This vector is the computation of the expression inside the parentheses of
#' the formula that is taken to the power \eqn{M}. Note that this vector is
#' "raw" and has not yet been taken to the power \eqn{M}. This vector is useful
#' for analyzing which output ages are most vulnerable to the MAP Decoder
#' attack.
#' - `sub.supp`: The subsample of the support. If the `subsample` argument is
#' left at its default value of 1, then this is just a vector of the integers
#' from 1 to the length of `f_S`.
#' If `return.only.success.probability.vector` is specified as TRUE, only
#' the `success.probability.vector` vector is returned, without a list.
#' @export
#'
#' @examples
#' n.outputs <- 1000
#' 
#' real.spend <- dnbinom(seq_len(n.outputs), size = 1, mu = 500)
#' decoy <- dnbinom(seq_len(n.outputs), size = 5, mu = 500)
#' # Negative binomial distributions
#' 
#' real.spend <- real.spend/sum(real.spend)
#' decoy <- decoy/sum(decoy)
#' 
#' plot(decoy, type = "l", ylim = c(0, max(c(real.spend, decoy))),
#'   ylab = "Probability mass function")
#' lines(real.spend, col = "red")
#' legend("topright", legend = c("Decoy", "Real spend"),
#'   col = c("black", "red"), lty = 1)
#'   
#' result <- map.decoder.success.prob(f_S = real.spend, f_D = decoy,
#'   n.decoys = 15)
#' 
#' print(result$mean.success.probability)
#' 
map.decoder.success.prob <- function(f_S, f_D, n.decoys, subsample = 1,
  rounding.digit = NULL, return.only.success.probability.vector = FALSE) {
  
  stopifnot(length(f_S) == length(f_D))
  
  stopifnot(isTRUE(all.equal(sum(f_S), 1)))
  stopifnot(isTRUE(all.equal(sum(f_D), 1)))
  
  if (length(subsample) == 1) {
    
    stopifnot(subsample %between% c(0, 1))
    
    if (subsample != 1) {
      sub.supp <- wrswoR::sample_int_expj(length(f_S),
        ceiling(length(f_S) * subsample), prob = f_S)
      # wrswoR::sample_int_expj is much faster than sample() when the prob vector is long
      sub.supp <- sort(sub.supp)
      f_S <- f_S[sub.supp]
      f_S <- f_S/sum(f_S)
      f_D <- f_D[sub.supp]
      f_D <- f_D/sum(f_D)
    } else {
      if (! return.only.success.probability.vector) {
        sub.supp <- seq_along(f_S)
      }
    }
    
  } else {
    sub.supp <- subsample
    f_S <- f_S[sub.supp]
    f_S <- f_S/sum(f_S)
    f_D <- f_D[sub.supp]
    f_D <- f_D/sum(f_D)
  }
  
  # START MAIN COMPUTATION
  
  # This implementation assumes that every output selected in a ring
  # is unique (i.e. does not handle multi-output aggregates like 
  # pseudo-blocks). The inequality comparison is strict.
  
  cut.vector <- f_S/f_D
  # rm(f_S)
  if (!is.null(rounding.digit)) {
    cut.vector <- Rfast::Round(cut.vector, digit = rounding.digit)
  }
  
  y <- data.table(f_D = f_D, cut.vector.var = cut.vector)
  
  setorder(y, cut.vector.var)
  
  cut.vector.copied <- cut.vector
  cut.vector <- data.table(cut.vector = cut.vector, cut.vector.ind = seq_along(cut.vector))
  cut.vector <- collapse::funique(cut.vector, cols = "cut.vector", sort = TRUE)
  
  y[, `:=`(cut.vector.cut, cut.vector$cut.vector.ind[.bincode(cut.vector.var,
    c(-1, cut.vector$cut.vector), right = FALSE)])]
  
  setorder(y, cut.vector.var)
  y[, cut.vector.var := NULL]
  
  y[, f_D := cumsum(f_D)]
  rm(f_D)
  y <- collapse::flast(y, collapse::GRP(y, by = "cut.vector.cut"))
  setnames(y, "f_D", "success.prob")
  
  y <- merge(y, data.table(cut.vector.cut = cut.vector$cut.vector.ind,
    cut.vector = cut.vector$cut.vector))
  y[, cut.vector.cut := NULL]
  
  y <- merge(y, data.table(cut.vector.cut.names = seq_along(cut.vector.copied),
    cut.vector = unname(cut.vector.copied)), all = TRUE, by = "cut.vector")
  y[, cut.vector := NULL]
  setorder(y, cut.vector.cut.names)
  y[, cut.vector.cut.names := NULL]
  
  setDF(y)
  
  y$success.prob <- collapse::replace_na(y$success.prob, value = 0)
  # At the point(s) where f_S/f_D is at a minimum, the attack would always
  # choose another block height, so the attack success probability is zero
  
  y <- y$success.prob
  
  # END MAIN COMPUTATION
  
  stopifnot(length(y) == length(f_S))
  
  if (return.only.success.probability.vector) {
    return(y)
  }
  
  weight <- f_S
  
  mean.success.probability <- sum(weight * y^n.decoys)
  
  list(
    mean.success.probability = mean.success.probability,
    success.probability.vector = y,
    sub.supp = sub.supp)
  
}
