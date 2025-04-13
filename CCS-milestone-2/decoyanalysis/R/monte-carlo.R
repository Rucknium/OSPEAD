


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





#' Create functions to draw random values from old and new DSAs and real spend distribution
#'
#' @details
#' Creates functions to draw random values from the old and new Decoy Selection
#' Algorithms (DSAs) and the real spend distribution.
#' 
#' @return A list with 3 elements: `draw.decoy.old.dsa`, `draw.decoy.new.dsa`,
#' and `draw.real.spend`. These elements are functions that take a single
#' argument: the number of random values to be drawn.
#' @export
#'
#' @examples
#' 1
#' 
gen.old.new.dsa.distributions <- function() {
  
  v <- 0.882831895764746
  z <- 122325532
  
  GAMMA_SHAPE = 19.28
  GAMMA_RATE = 1.61
  DIFFICULTY_TARGET_V2 = 120
  RECENT_SPEND_WINDOW = 15 * DIFFICULTY_TARGET_V2
  
  G_star <- function(x) {
    decoyanalysis::wallet2_gamma_cdf(as.numeric(x), v, z, GAMMA_SHAPE, GAMMA_RATE, RECENT_SPEND_WINDOW)
  }
  
  get.decoy.pmf <- function(cdf, v, z, sub.supp, log.trans = FALSE, ...) {
    
    G <- function(x) {
      cdf(x, ...)
    }
    
    if (log.trans) {
      G_star <- function(x) {
        G(log1p(x*v))/G(log1p(z*v))
      }
    } else {
      G_star <- function(x) {
        G(x*v)/G(z*v)
      }
    }
    
    G_star(sub.supp) - G_star(sub.supp - 1)
    
  }
  
  f_D.lgamma <- function(param, v, z, sub.supp, get.decoy.pmf) {
    list(decoy.pmf =  get.decoy.pmf(actuar::plgamma, v, z, sub.supp = sub.supp, shapelog = exp(param[1]), ratelog = exp(param[2])),
      tail.beyond.support = 1 - actuar::plgamma(z*v, shapelog = exp(param[1]), ratelog = exp(param[2])))
  }
  
  
  f_D.fun <- f_D.lgamma
  fitted.par <- log(c(5.3915659, 0.5428289))
  
  message("Generating new decoy PMF...")
  new_decoy <- stepfun(1:z, cumsum(c(0, f_D.fun(fitted.par, v, z, sub.supp = 1:z, get.decoy.pmf)$decoy.pmf)))
  new_decoy.pmf <- diff(c(0, new_decoy(1:z)))
  
  message("Generating old decoy PMF...")
  G_star.pmf <- diff(c(0, G_star(1:z)))
  
  new_decoy.pmf <- new_decoy.pmf/sum(new_decoy.pmf)
  G_star.pmf <- G_star.pmf/sum(G_star.pmf)
  
  
  real.spend.cdf.points <- list(
    supp.points.reg = c(0.000575456426584267, 0.0109378622177646, 
      0.0203469960769643, 0.0311601866522138, 0.0467999151458647, 0.0609107611266039, 
      0.0734149022024337, 0.0847949320361706, 0.0950801101525983, 0.104437640059852, 
      0.114196913907047, 0.122645929365097, 0.130745360385555, 0.138179897372298, 
      0.145299567138383, 0.151215168847502, 0.161301523281243, 0.173764241657995, 
      0.184505786879567, 0.193736849257658, 0.203277785466477, 0.212444406541027, 
      0.220600235751461, 0.231854557849399, 0.2429001913787, 0.252959566218011, 
      0.262425090373765, 0.270224400197801, 0.28108079969358, 0.290922234060782, 
      0.300530983952483, 0.310446268346829, 0.320966815447017, 0.329995942638808, 
      0.337746083109144, 0.347471960694376, 0.360273584613406, 0.371891512808412, 
      0.380192777892175, 0.390910493636414, 0.400857411840069, 0.409347284375441, 
      0.420328182706099, 0.430372573878364, 0.440387061627405, 0.450586209405976, 
      0.459244604578224, 0.470313868463874, 0.480061324374701, 0.490286511298928, 
      0.500549865607883, 0.509700581274239, 0.520393267841975, 0.530294241189403, 
      0.539798228578179, 0.549794019799906, 0.560485279197046, 0.570328132477984, 
      0.580241854348751, 0.590225419808097, 0.600534288716839, 0.609673187275634, 
      0.620265104850399, 0.630138413433524, 0.63871319966383, 0.650137420030894, 
      0.659387496450059, 0.669888740906512, 0.680975161460898, 0.690310725025491, 
      0.700348534742927, 0.710153413121885, 0.719669981407517, 0.730395994265596, 
      0.738745011501996, 0.750004317146264, 0.759935875879115, 0.770046712519537, 
      0.779982135325032, 0.789653447570464, 0.799933844981753, 0.809973240889287, 
      0.819447173761549, 0.830177595920535, 0.840627487373578, 0.850325129294852, 
      0.859750591041841, 0.869738362244709, 0.879746333670099, 0.890192934500753, 
      0.899727591486718, 0.909528938716729, 0.919305653325935, 0.930043059255659, 
      0.939757820503097, 0.949472202879246, 0.959924191462532, 0.96834416415755, 
      0.976260609542887, 0.987475513621487, 0.999031813759897),
    Fs.hat.transformed.reg = c(2, 
      16, 30, 47, 75, 103, 130, 156, 182, 208, 234, 259, 285, 311, 
      337, 362, 407, 458, 510, 562, 614, 666, 719, 791, 871, 952, 1033, 
      1115, 1224, 1334, 1444, 1578, 1733, 1890, 2041, 2252, 2490, 2734, 
      2982, 3269, 3616, 3999, 4460, 4904, 5341, 5835, 6389, 7099, 7979, 
      8923, 9927, 10897, 11977, 13283, 14757, 16673, 19038, 22072, 
      25210, 29200, 33811, 39785, 46357, 54352, 61647, 69708, 81037, 
      96580, 113202, 128582, 143390, 162596, 181177, 205854, 227296, 
      242268, 259145, 298222, 343149, 392645, 451052, 510761, 592332, 
      674226, 760117, 851197, 1077279, 1528320, 1957672, 2449500, 2971191, 
      3993817, 5368590, 7066115, 18753479, 24951735, 37608506, 42692673, 
      48938384, 56748829, 116047780))
  
  real_spend <- suppressWarnings(
    with(real.spend.cdf.points,
      binsmooth::splinebins(
        bEdges = c(Fs.hat.transformed.reg, max(Fs.hat.transformed.reg) + 1),
        bCounts = c(diff(c(0, supp.points.reg)), 0))$splineCDF
    )
  )
  
  message("Generating real spend PMF...")
  real_spend.pmf <- diff(c(0, real_spend(1:z)))
  
  real_spend.pmf[real_spend.pmf <= 0] <- .Machine$double.eps
  
  real_spend.pmf <- real_spend.pmf/sum(real_spend.pmf)
  
  draw.real.spend <- function(n) {
    sample(1:z, size = n, replace = TRUE, prob = real_spend.pmf)
  }
  
  draw.decoy.new.dsa <- function(n) {
    sample(1:z, size = n, replace = TRUE, prob = new_decoy.pmf)
  }
  
  draw.decoy.old.dsa <- function(n) {
    sample(1:z, size = n, replace = TRUE, prob = G_star.pmf)
  }
  
  list(
    draw.decoy.old.dsa = draw.decoy.old.dsa,
    draw.decoy.new.dsa = draw.decoy.new.dsa,
    draw.real.spend = draw.real.spend
  )
  
}





gen.nonfungible.rings <- function(draw.my.dsa, draw.their.dsa, draw.real.spend,
  n, beta, C, n.rings, check.data.sorted = FALSE) {
  
  non.change.real.spend <- sample(c(1L, 0L), size = n.rings, replace = TRUE, prob = c(beta, 1 - beta))
  # Draw n.rings elements with replacement from the set of {1, 0} with probability beta of drawing 1 and
  # probability 1-beta of drawing 0. The 1 is a real spend output that has the defect.
  
  real.spend <- ifelse(
    sample(c(TRUE, FALSE), size = n.rings, replace = TRUE, prob = c(C, 1 - C)),
    1L,
    non.change.real.spend)
  # With probability C, the ring spends change. The change has the defect by assumption. With probability
  # 1-C, the user spend a non-change output that has beta probability (above) of having the defect.
  
  rings <- matrix(c(
    sample(c(1L, 0L), size = n.rings * (n - 1), replace = TRUE, prob = c(beta, 1 - beta)),
    real.spend),
    nrow = n.rings, ncol = n)
  # Create a matrix n.rings x n in size. The first n-1 columns are filled with decoys. With probability
  # beta these decoys have the defect. The last column is the real spend, created above.
  
  
  my.dsa.rings.to.produce <- sum(c(rings) == 1)
  their.dsa.rings.to.produce <- sum(c(rings) == 0)
  
  draw.my.dsa.all <- draw.my.dsa(n.rings * (n - 1) + my.dsa.rings.to.produce * (n - 1))
  draw.real.spend.all <- draw.real.spend(n.rings + my.dsa.rings.to.produce + their.dsa.rings.to.produce)
  # Drawing all at once is much faster
  
  
  # Final column will always be the real spend
  ring.distribution <- cbind(
    matrix(draw.my.dsa.all[1:(n.rings * (n - 1))], nrow = n.rings),
    matrix(draw.real.spend.all[1:n.rings], nrow = n.rings)
  )
  
  antecedent.ring.distribution.mine <- rbind(
    matrix(draw.my.dsa.all[n.rings * (n - 1) + seq_len(my.dsa.rings.to.produce * (n - 1))], ncol = my.dsa.rings.to.produce),
    matrix(draw.real.spend.all[n.rings + seq_len(my.dsa.rings.to.produce)], ncol = my.dsa.rings.to.produce)
  )
  
  rm(draw.my.dsa.all)
  
  antecedent.ring.distribution.theirs <- rbind(
    matrix(draw.their.dsa(their.dsa.rings.to.produce * (n - 1)), ncol = their.dsa.rings.to.produce),
    matrix(draw.real.spend.all[n.rings + my.dsa.rings.to.produce + seq_len(their.dsa.rings.to.produce)], ncol = their.dsa.rings.to.produce)
  )
  
  
  rm(draw.real.spend.all)
  
  antecedent.ring.distribution.index <- replicate(n, rings)
  # replicate() will bind the duplicate matrices together along the third dimension
  # Don't need Reduce(abind, list(rings, rings,...))
  
  antecedent.ring.distribution.dims <- dim(antecedent.ring.distribution.index)
  
  antecedent.ring.distribution.index <- c(antecedent.ring.distribution.index)
  
  antecedent.ring.distribution <- vector("integer", length(antecedent.ring.distribution.index))
  
  stopifnot(length(antecedent.ring.distribution) ==
      length(antecedent.ring.distribution.mine) + length(antecedent.ring.distribution.theirs))
  
  antecedent.ring.distribution[antecedent.ring.distribution.index == 1L] <- antecedent.ring.distribution.mine
  antecedent.ring.distribution[antecedent.ring.distribution.index == 0L] <- antecedent.ring.distribution.theirs
  
  stopifnot( ! any(is.na(antecedent.ring.distribution)))
  
  
  antecedent.ring.distribution <- array(antecedent.ring.distribution, dim = antecedent.ring.distribution.dims)
  
  antecedent.ring.distribution <- apply(antecedent.ring.distribution, 2, "c")
  
  antecedent.ring.distribution <- Rfast::rowSort(antecedent.ring.distribution)
  
  antecedent.ring.distribution <- array(c(antecedent.ring.distribution), dim = antecedent.ring.distribution.dims)
  
  antecedent.ring.distribution <- aperm(antecedent.ring.distribution, perm = c(1, 3, 2))
  
  if (check.data.sorted) {
    stopifnot(all( ! c(apply(antecedent.ring.distribution, c(1,3), is.unsorted))))
  }
  
  
  ring.order <- Rfast::rowOrder(ring.distribution)
  
  ring.distribution <- Rfast::rowSort(ring.distribution)
  
  
  for (i in seq_len(nrow(ring.order))) {
    antecedent.ring.distribution[i, , ] <-
      antecedent.ring.distribution[i, , ring.order[i, ]]
  }
  
  
  antecedent.ring.distribution.for.classification <- t(apply(antecedent.ring.distribution, 1, "c"))
  
  grid.names <- expand.grid(0:(n-1), 0:(n-1))
  
  colnames(antecedent.ring.distribution.for.classification) <- 
    paste0("Inputs.0.Previous_Tx_Block_Num_Delta.", grid.names[, 2], ".", grid.names[, 1])
  # Need "slowest" varying label first, which is the 2nd column of grid.names
  # First is the ring member slot of the tx of interest.
  # Second is the ring member slot of the antecedent tx
  
  colnames(ring.distribution) <- 
    paste0("Inputs.0.Previous_Tx_Block_Num_Delta.self.", 0:(n-1))
  
  result <- cbind(ring.distribution, antecedent.ring.distribution.for.classification)
  
  if (check.data.sorted) {
    stopifnot(all( ! apply(result[, paste0("Inputs.0.Previous_Tx_Block_Num_Delta.self.", 0:(n-1))],
      1, is.unsorted)))
    stopifnot(all( ! apply(result[, paste0("Inputs.0.Previous_Tx_Block_Num_Delta.0.", 0:(n-1))],
      1, is.unsorted)))
  }
  
  list(X = result, y = apply(ring.order, 1, FUN = function(x) which(x == n)) - 1)
  # y is now indexed from zero
  
}






#' Generate rings from old and new DSAs and their antecedents
#'
#' @param n.rings Number of rings to generate
#' @param beta Share of rings that are created by the new DSA
#' @param old.dsa Function that draws a random variable from the old DSA. It
#' should have a single argument: the number of values to draw.
#' @param new.dsa Function that draws a random variable from the new DSA. 
#' @param real.spend Function that draws a random variable from the real spend
#' distribution.
#' @param ring.size Ring size. 16 by default.
#' @param change.probability Probability that a ring spends the change in the
#' wallet from a prior transaction.
#' @param check.data.sorted Double-check if the ring data is properly sorted?
#' 
#' @details
#' 
#' 
#' For more information on variable definitions and the behavioral model,
#' see "Formula for Accuracy of Guessing Monero Real Spends Using Fungibility
#' Defects":
#' \url{https://github.com/Rucknium/misc-research/blob/main/Monero-Fungibility-Defect-Classifier/pdf/classify-real-spend-with-fungibility-defects.pdf}
#'
#' @return A list with three elements:
#' 
#' `X`: a matrix with `n.rings` rows and `ring.size + ring.size^2` columns of
#' the output ages of each ring member and its antecedents. The column
#' "Inputs.0.Previous_Tx_Block_Num_Delta.self.i" is the age of the `i`th ring
#' member of the transaction of interest.
#' "Inputs.0.Previous_Tx_Block_Num_Delta.j.k" is the age of the `k`th
#' ring member of the antecedent transaction referenced in the `j`th ring
#' member slot of the transaction of interest.
#' 
#' `y`: a vector of length `n.rings` containing the ring member position of
#' the real spend, indexed from zero.
#' 
#' `z`: a vector of length `n.rings` with zero where the ring is an old-DSA
#' ring and one where the ring is a new-DSA ring.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' distributions <- gen.old.new.dsa.distributions()
#' 
#' ring.data <- gen.old.new.dsa.rings(n.rings = 10000, beta = 0.10,
#'   old.dsa = distributions$draw.decoy.old.dsa,
#'   new.dsa = distributions$draw.decoy.new.dsa,
#'   real.spend = distributions$draw.real.spend)
#' 
#' write.csv(ring.data$X, file = "X.csv", row.names = FALSE)
#' cat(ring.data$y, file = "y.csv", sep = "\n")
#' cat(ring.data$z, file = "z.csv", sep = "\n")
#' }
gen.old.new.dsa.rings <- function(n.rings, beta, old.dsa, new.dsa, real.spend,
  ring.size = 16, change.probability = 0.40, check.data.sorted = FALSE) {
  
  if (beta <= 0 || beta >= 1) {
    stop("beta must be strictly between 0 and 1.")
  }
  
  n.rings.old <- ceiling(n.rings * (1 - beta) )
  n.rings.new <- ceiling(n.rings * beta)
  
  message("Generating old-DSA rings...")
  rings.old <- gen.nonfungible.rings(draw.my.dsa = old.dsa, draw.their.dsa = new.dsa,
    draw.real.spend = real.spend, n = ring.size, beta = 1 - beta,
    C = change.probability, n.rings = n.rings.old,
    check.data.sorted = check.data.sorted)
  
  message("Generating new-DSA rings...")
  rings.new <- gen.nonfungible.rings(draw.my.dsa = new.dsa, draw.their.dsa = new.dsa,
    draw.real.spend = real.spend, n = ring.size, beta = beta,
    C = change.probability, n.rings = n.rings.old,
    check.data.sorted = check.data.sorted)
  
  list(
    X = rbind(rings.old$X, rings.new$X),
    y = c(rings.old$y, rings.new$y),
    z = c(rep(0, n.rings.old), rep(1, n.rings.new))
  )
  
}
































