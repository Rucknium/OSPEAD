# install.packages("data.table")
# install.packages("actuar")
# install.packages("distributionsrd")
# Un-comment lines to install packages

library(data.table)

xmr <- read.csv("output_age_data.csv", stringsAsFactors = FALSE)
# Obtain data from https://github.com/monero-project/monero/files/6968268/output_age_data.zip

xmr$Observed.pdf <- xmr$Observed/sum(xmr$Observed)
# Convert f(x) to a "probability density function", more or less
xmr$Current.decoy.selection.algo.pdf <- xmr$Current.decoy.selection.algo / 
  sum(xmr$Current.decoy.selection.algo)
# Do the same for f_M(x)

alpha <- 10/11

xmr$f_S <- (1/(1-alpha)) * 
  (xmr$Observed.pdf - alpha * xmr$Current.decoy.selection.algo.pdf)
# Construct f_S(x) according to equation (2) in the paper

xmr$f_S[xmr$f_S < 0] <- xmr$Current.decoy.selection.algo.pdf[xmr$f_S < 0]
# Sometimes the observed can be below what the idealized algorithm would have 
# selected since there is noise in the long tails. 

xmr$f_S <- xmr$f_S / sum(xmr$f_S)
# Make f_S(x) be a proper PDF



guess.prob.efficient <- function(f_S, f_D) {
  
  cut.vector <- f_S/f_D
  names(cut.vector) <- as.character(1:length(cut.vector))
  
  y <- data.table(f_D = f_D, cut.vector.var = cut.vector)
  
  setorder(y, cut.vector.var)
  
  cut.vector <- sort(cut.vector)
  
  cut.vector.name.unique <- rownames(unique(matrix(cut.vector, dimnames = list(names(cut.vector), NULL))))
  # Use matrix() since unique() doesnt carry names otherwise
  
  y[, cut.vector.cut := cut(cut.vector.var, c(-1, unique(cut.vector)), 
    labels = cut.vector.name.unique, right = FALSE)]
  
  setorder(y, cut.vector.var)
  
  y[, f_D.cumsum := cumsum(f_D)]
  y <- y[, .(success.prob = f_D.cumsum[.N]), by = cut.vector.cut]
  y <- merge(y, data.table(cut.vector.cut = cut.vector.name.unique, cut.vector = unique(cut.vector)))
  
  y <- merge(y, data.table( cut.vector.cut.names = names(cut.vector), cut.vector = cut.vector ), 
    all = TRUE, by = "cut.vector")
  y[, cut.vector.cut.names := as.numeric(cut.vector.cut.names)]
  setorder(y, cut.vector.cut.names)
  
  y$success.prob[is.na(y$success.prob)] <- 0
  # At the point(s) where f_S/f_D is at a minimum, the attack would always
  # choose another block height, so the attack success probability is zero
  
  setDF(y)
  
  y
  
}


guess.prob.inefficient <- function(f_S, f_D) {
  y <- vector(mode = "numeric", length = length(f_D))
  for (i in seq_along(f_D)) {
    y[i] <- sum(f_D  * (f_S/f_D < (f_S/f_D)[i]))
  }
  y
}



x.support <- seq_len(1e5)

f_D <- actuar::dlgamma(x.support + 1, shapelog = 6.48, ratelog = 0.894)
f_S <- distributionsrd::drightparetolognormal(x.support, 
  shape2 = 0.235, meanlog = 4.12, sdlog = 1.38)
# Simple simulation with some parameters from "dry run" FGT, alpha = 1

f_D <- f_D/sum(f_D)
f_S <- f_S/sum(f_S)
# Make into valid probability measures

system.time(values.guess.prob.efficient <- guess.prob.efficient(f_S, f_D))
system.time(values.guess.prob.inefficient <- guess.prob.inefficient(f_S, f_D))

n.decoys <- 10

sum(f_S * (values.guess.prob.efficient$success.prob)^n.decoys)
# [1] 0.1017509

sum(f_S * (values.guess.prob.inefficient)^n.decoys)
# [1] 0.1017509

sum(f_S * (values.guess.prob.efficient$success.prob)^n.decoys) -
  sum(f_S * (values.guess.prob.inefficient)^n.decoys)
# 0
# The two implementations give the same value



xmr$Current.decoy.selection.algo.pdf[xmr$Current.decoy.selection.algo.pdf == 0] <- .Machine$double.eps
# Give decoy at least a nonzero, negligible probability at every block height

guess.prob.Moser.data <- guess.prob.efficient(
  f_S = xmr$f_S, f_D = xmr$Current.decoy.selection.algo.pdf)

sum(xmr$f_S * (guess.prob.Moser.data$success.prob)^n.decoys)
# [1] 0.3543693
# Very close to the 0.354512 estimated by Monte Carlo simulation
# in my HackerOne submission

n.decoys <- 15

sum(xmr$f_S * (guess.prob.Moser.data$success.prob)^n.decoys)
# [1] 0.3039381
# For ring size = 16



