
# install.packages("distr")
# install.packages("fitdistrplus")
# install.packages("ks")
# library(ks)
# library(distr)
# library(fitdistrplus)

set.seed(314)
n.obs <- 100000
x.trunc <- 5

mixed.dist <- distr::UnivarMixingDistribution(
  distr::Exp(20), 
  distr::Chisq(df = 2),
  mixCoeff = c(0.2, 0.8))

r.mixed.dist <- distr::r(mixed.dist)
d.mixed.dist <- distr::d(mixed.dist)

x <- r.mixed.dist(n.obs)
x <- x[x <= x.trunc]
length(x)
summary(x)
hist(x, breaks = 200, freq = FALSE, col = "blue", border = NA)


mle.gamma <- fitdistrplus::fitdist(x, "gamma", start = list(shape = 2, scale = 2))

png("Histogram-and-fitted-gamma-PDF.png", width = 600, height = 600)

hist(x, breaks = 200, freq = FALSE, col = "blue", border = NA,
  main = "Histogram and gamma MLE fitted PDF")

lines(
  seq(0, x.trunc, by = 0.01),
  dgamma(seq(0, x.trunc, by = 0.01), shape = mle.gamma$estimate["shape"], scale = mle.gamma$estimate["scale"]), 
  col = "red", lwd = 2)

dev.off()



kde.est <- ks::kde(x, gridsize = 5000) 

x.seq <- seq(0, x.trunc, length.out = 500)


png("Comparison-of-PDF-estimates.png", width = 600, height = 600)

plot(0, 0, xlim = c(0, x.trunc), 
  ylim = c(0, 2),
  col = "transparent",
  xlab = "x",
  ylab = "Density",
  main = "Comparison of PDF estimates")

legend( "topright",
  legend = c("True PDF", "Gamma MLE fitted PDF", "Kernel density estimate PDF"),
  col = c("black", "red", "green"),
  lty = 1, lwd = 2)

lines(x.seq, d.mixed.dist(x.seq), col = "black", lwd = 2)
lines(x.seq,
  dgamma(x.seq, shape = mle.gamma$estimate["shape"], scale = mle.gamma$estimate["scale"]), 
  col = "red", lwd = 2)
lines(kde.est$eval.points[kde.est$eval.points > 0], 
  kde.est$estimate[kde.est$eval.points > 0], col = "green", lwd = 2)


dev.off()


ks.test.gamma<- ks.test(x, pgamma, shape = mle.gamma$estimate["shape"], scale = mle.gamma$estimate["scale"])
ks.test.gamma

png("Comparison-of-CDF-gamma.png", width = 600, height = 600)

plot(ecdf(x), do.points = FALSE, col = "black", lwd = 2, 
  main = "Comparison of empirical cdf of\nsimulated data and gamma MLE fitted CDF",
  sub = paste0("K-S test statistic: ", round(ks.test.gamma$statistic, digits = 5))) 

lines(seq(0, x.trunc, by = 0.01), col = "red", lwd = 2,
  pgamma(seq(0, x.trunc, by = 0.01), shape = mle.gamma$estimate["shape"], scale = mle.gamma$estimate["scale"]))

legend( "bottomright",
  legend = c("Empirical CDF", "Gamma MLE fitted CDF"),
  col = c("black", "red"),
  lty = 1, lwd = 2)

dev.off()


kde.cdf.est <- ks::kcde(x, gridsize = 5000)
kde.cdf.est.fun <- approxfun(kde.cdf.est$eval.points, kde.cdf.est$estimate, yleft = 0, yright = 0)


ks.test.kde <- ks.test(x, kde.cdf.est.fun)
ks.test.kde

png("Comparison-of-CDF-kernel-density.png", width = 600, height = 600)

plot(ecdf(x), do.points = FALSE, col = "black", lwd = 2, 
  main = "Comparison of empirical cdf of\nsimulated data and kernel density estimate CDF",
  sub = paste0("K-S test statistic: ", round(ks.test.kde$statistic, digits = 5))) 

lines(seq(0, x.trunc, by = 0.01),
  kde.cdf.est.fun(seq(0, x.trunc, by = 0.01)), col = "red", lwd = 2) 

legend( "bottomright",
  legend = c("Empirical CDF", "Kernel density estimate CDF"),
  col = c("black", "red"),
  lty = 1, lwd = 2)

dev.off()


