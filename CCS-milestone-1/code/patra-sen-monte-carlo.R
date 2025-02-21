# install.packages("Iso")
# install.packages("kde1d")
# Un-comment lines to install packages


mc.ring.size <- c(11, 16, 128)
names(mc.ring.size) <- mc.ring.size

n.rings <- 100000
n.monte.carlo <- 1000

decoy.param <- c(2, 150)
real.param <- c(3, 100)

Fb <- function(x) {
  stats::pnbinom(x, size = decoy.param[1], mu = decoy.param[2])
}

discard.above.support <- 500

mc.results <- as.list(mc.ring.size)

for (ring.size.i in names(mc.results)) {
  
  n.decoys <- as.numeric(ring.size.i) - 1
  alp.hat <- 1/(n.decoys + 1)
  
  intermediate.results <- list()
  
  set.seed(314)
  
  for (intermediate.i in seq_len(n.monte.carlo)) {
    
    decoys <- stats::rnbinom(n.rings * n.decoys, size = decoy.param[1], mu = decoy.param[2])
    real <- stats::rnbinom(n.rings, size = real.param[1], mu = real.param[2])
    
    data.mix <- c(decoys, real)
    
    # Modified code based on code provided with
    # "Estimation of a Two-component Mixture Model with Applications to Multiple Testing"
    # https://users.stat.ufl.edu/~rohitpatra/research.html
    
    data.mix <- sort(data.mix) ## Sorts the data set
    data.1 <- unique(data.mix) ## Finds the unique data points
    Fn <- ecdf(data.mix) ## Computes the empirical DF of the data
    Fn.1 <- Fn(data.1) ## Empirical DF of the data at the data points
    pava.freq <- diff(c(0,Fn.1))
    
    F.hat <- (Fn.1 - (1 - alp.hat) * Fb(data.1)) / alp.hat ## Computes the naive estimator of F_s
    Fs.hat <- Iso::pava(F.hat, pava.freq, decreasing = FALSE) ## Computes the Isotonic Estimator of F_s
    Fs.hat[which(Fs.hat <=0 )] <- 0
    Fs.hat[which(Fs.hat >= 1)] <- 1
    
    ## End code based on "Estimation of a Two-component..."
    
    kernel.points <- rep(min(data.1):(length(diff(Fs.hat)) - min(data.1)), 
      times = c(ifelse( round(diff(Fs.hat) * n.rings) > 0, round(diff(Fs.hat) * n.rings), 1), 0))
    
    est.kernel.dens <- kde1d::kde1d(ordered(kernel.points))
    est.kernel.dens.eval <- 
      kde1d::dkde1d(sort(unique(ordered(kernel.points)))[1:discard.above.support], est.kernel.dens)
    
    intermediate.results[[intermediate.i]] <- 
      data.frame(x = intersect(1:discard.above.support, kernel.points), y = est.kernel.dens.eval)
    
    cat(ring.size.i, intermediate.i, base::date(), "\n")
    
  }
  
  intermediate.results.merging <- intermediate.results[[1]]
  
  for (i in seq_along(intermediate.results)[-1]) {
    intermediate.results.merging <- merge(intermediate.results.merging, 
      intermediate.results[[i]], all = TRUE, by = "x")
  }
  
  intermediate.results.merging <- intermediate.results.merging[order(as.numeric(intermediate.results.merging$x)), ]
  
  mc.results[[ring.size.i]] <- intermediate.results.merging
  
}


support.x <- 20:200

estimated.CI <- lapply(mc.results, FUN = function(x) {
  x[is.na(x)] <- 0
  apply(x[support.x, -1], 1, FUN = quantile, na.rm = FALSE, probs = c(0.05, 0.95))
})

estimated.CI <- do.call(rbind, estimated.CI)


median.CI.width.11 <- median(estimated.CI[2, ] - estimated.CI[1, ])
median.CI.width.16 <- median(estimated.CI[4, ] - estimated.CI[3, ])
median.CI.width.128 <- median(estimated.CI[6, ] - estimated.CI[5, ])


median.CI.width.16 / median.CI.width.11
median.CI.width.128 / median.CI.width.16


png(filename = "patra-sen-monte-carlo.png", width = 1000, height = 1200)

par(mfrow = c(2, 2), mar = c(2, 4, 2, 2) + 0.1)

for (support.x in list(1:20, 20:100, 100:200, 200:500)) {
  
  true.density <- stats::dnbinom(support.x, size = real.param[1], mu = real.param[2])
  decoy.density <- stats::dnbinom(support.x, size = decoy.param[1], mu = decoy.param[2])
  
  estimated.CI <- lapply(mc.results, FUN = function(x) {
    x[is.na(x)] <- 0
    apply(x[support.x, -1], 1, FUN = quantile, na.rm = FALSE, probs = c(0.05, 0.95))
  })
  
  estimated.CI <- do.call(rbind, estimated.CI)
  
  stopifnot(length(mc.ring.size) == 3)
  # NOTE: Will have to adjust line.colors if length(mc.ring.size) changes
  line.colors <- c("red", "blue", "darkgreen")
  
  plot(x = support.x, true.density, type = "l", ylim = c(0, max(estimated.CI)),
    ylab = "Density", xlab = NULL,
    main = paste0("Support ", min(support.x), "-", max(support.x)))
  abline(h = 0)
  lines(x = support.x, decoy.density, lty = 1, col = "violet")
  matlines(x = support.x, t(estimated.CI), lty = rep(2, nrow(estimated.CI)),
    col = rep(line.colors, each = 2))
  
  if (1 %in% support.x) {
    legend("topleft", 
      legend = c("Real spend distribution",
        "Decoy distribution",
        paste0("90% C.I. with ring size ", mc.ring.size)),
      col = c("black", "violet", line.colors),
      lty = c(1, 1, rep(2, length(mc.ring.size))),
      bty = "n")
  }
  
}

dev.off()


