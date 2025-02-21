


#' Title
#'
#' @param N Sample size
#' @param T Number of variables (repeated measurements) 
#' @param K Number of components
#' @param theta1 Vector of mean parameters
#' @param theta2 Vector of standard deviation parameters
#' @param omega  Vector of mixing proportions
#'
#' @return N x T matrix
#' @export
#'
#' @examples
#' N = 500 # Sample size
#' T = 3   # Number of variables (repeated measurements) 
#' K = 2   # Number of components
#' theta1 = c(0.00, 3.00) # Mean parameters
#' theta2 = c(1.00, 1.00) # Standard deviation parameters
#' omega  = c(0.30, 0.70) # Mixing proportions
#' 
#' set.seed(314)
#' y <- gen.bjr.normal.data(N = N, T = T, K = K, theta1 = theta1, theta2 = theta2, omega = omega)
gen.bjr.normal.data <- function(N, T, K, theta1, theta2, omega) {
  stopifnot(all(
    c(length(theta1), length(theta2), length(omega)) == K)
  )
  y = matrix(0, nrow = N, ncol = T)
  u = sample(1:K, N, replace = TRUE, prob = omega)
  for (t in 1:T) { 
    y[, t] =  rnorm(N, mean = theta1[u], sd = theta2[u])
  } 
  y
}




#' Title
#'
#' @param n.rings S
#' @param ring.size S
#' @param mixing.proportions S
#' @param decoy.dists S
#' @param real.spend.dists S
#'
#' @return Matrix
#' @export
#'
#' @examples
#' 1
#' 
gen.ring.data <- function(n.rings, ring.size, mixing.proportions, decoy.dists, real.spend.dists) {
  stopifnot( length(unique(length(mixing.proportions),
    length(decoy.dists), length(real.spend.dists))) == 1)
  # Check that the length of these arguments are all the same (same number of distribution components)
  stopifnot(ring.size >= 3)
  # bjr() can only be estimated when ring size is 3 or greater
  dist.type.draws <- c(stats::rmultinom(1, size = n.rings, prob = mixing.proportions))
  # This will not give exactly the same proportion as the theoretical mixing.proportions
  # parameter, but it is consistent with sampling from a population
  
  rings <- vector("list", length(mixing.proportions))
  decoy.label <- vector("list", length(mixing.proportions))
  for (dist.type in seq_along(mixing.proportions)) {
    n.decoy.draws <- (ring.size - 1) * dist.type.draws[dist.type]
    decoy.draws <- decoy.dists[[dist.type]]( n.decoy.draws )
    rings.component <- matrix(decoy.draws, nrow = dist.type.draws[dist.type], ncol = (ring.size - 1))
    n.real.spend.draws <- dist.type.draws[dist.type]
    rings.component <- cbind(rings.component, real.spend.dists[[dist.type]]( n.real.spend.draws ))
    rings[[dist.type]] <- t(apply(rings.component, 1, sample))
    # Randomize observations within the ring just to be safe.
  }
  rings <- do.call(rbind, rings)
  random.row.order <- sample(nrow(rings))
  rings <- rings[random.row.order, ]
  # Randomize the order of rings in the matrix.
  component.label <- rep(seq_along(dist.type.draws), times = dist.type.draws)
  component.label <- component.label[random.row.order]
  return(list(rings = rings, component.label = component.label, decoy.label = decoy.label))
  # TODO: Complete the decoy label
}


#' Title
#'
#' @return Matrix
#' @export
#'
#' @examples
#' 1
#' 
gen.standard.bjr.test.dataset <- function() {
  
  mixing.proportions <- c(0.30, 0.70)
  
  decoy.dists <- list(
    function(x) { rnorm(x, 0, 1)},
    function(x) { rnorm(x, 3, 1)}
  )
  real.spend.dists <- list(
    function(x) { rnorm(x, 0, 1)},
    function(x) { rnorm(x, 3, 1)}
  )
  
  N = 500 # Sample size
  M = 3   # Number of variables (repeated measurements) 
  K = length(mixing.proportions)   # Number of components
  
  set.seed(314)
  ring.data <- gen.ring.data(n.rings = N, ring.size = M, mixing.proportions, decoy.dists, real.spend.dists)
  y <- ring.data$rings
  y
}










