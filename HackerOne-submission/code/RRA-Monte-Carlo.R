
xmr <- read.csv("output_age_data.csv", stringsAsFactors = FALSE)
# From https://github.com/monero-project/monero/files/6968268/output_age_data.zip

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

xmr$Rucknium.Ratio <- xmr$f_S / xmr$Current.decoy.selection.algo.pdf

png("Rucknium-Ratio.png", width = 600, height = 600)

par(mar =  c(5, 6, 4, 2) + 0.1)
#  c(5, 4, 4, 2) + 0.1

plot(xmr$Rucknium.Ratio[1:10000], log = "x", cex = 0.5, 
  col = rgb(0, 0, 0, alpha = c(rep(1, 100), rep(0, 9900))),
  main = "Rucknium Ratio\nDots are partially transparent for better visualization",
  xlab = "Age of blocks corresponding to ring members (log scale)",
  ylab = expression("Rucknium Ratio   " * frac(f[S](x), f[M](x))) )
lines(x = 1:100, y = xmr$Rucknium.Ratio[1:100])
points(x = 101:10000, y = xmr$Rucknium.Ratio[101:10000], col = rgb(0, 0, 0, alpha = 0.2), pch = ".")
abline(h = 1, lty = 2, col = "red")
axis(2, 1, "1")

dev.off()

par(mar =  c(5, 4, 4, 2) + 0.1 )


sum(xmr[1:10, "f_S"])
# [1] 0.1380584
sum(xmr[1:10, "Current.decoy.selection.algo.pdf"])
# [1] 0.02727423

sum(xmr[1:3, "f_S"])
# [1] 0.05523403
sum(xmr[1:3, "Current.decoy.selection.algo.pdf"])
# [1] 0.005761155


set.seed(314)

n.rings <- 10000000
# 10 million

simulation.mixin.quantity <- 10

rings <- matrix(c(
  sample(nrow(xmr), size = simulation.mixin.quantity * n.rings, replace = TRUE, 
    prob = xmr$Current.decoy.selection.algo.pdf),
  sample(nrow(xmr), size = 1 * n.rings, replace = TRUE, prob = xmr$f_S)
), byrow = FALSE, ncol = simulation.mixin.quantity + 1)


attack.prob.11 <- apply(rings, 1, FUN = function(x) {
  which.max(xmr$Rucknium.Ratio[x])
})
# Note that which.max() returns the index of the "first" maximum if the maximum
# is not unique. Therefore, the real spend was chosen to be the "last"
# element of the set so that it would be clear that the real spend, if
# guessed, was the result of unique guess

t(t(100 * prop.table(table(attack.prob.11)) ))
# Main result


##############################################################
# Simulation for different ring sizes follows
# Note: This takes about an hour to run
##############################################################

set.seed(314)

n.rings <- 1000000
# 1 million

max.ring.size <- 256

guessibility.results <- list()

for (i in seq_len(max.ring.size - 1)) {
  
  simulation.mixin.quantity <- i
  
  rings <- matrix(c(
    sample(nrow(xmr), size = simulation.mixin.quantity * n.rings, replace = TRUE, prob = xmr$Current.decoy.selection.algo.pdf),
    sample(nrow(xmr), size = 1 * n.rings, replace = TRUE, prob = xmr$f_S)
  ), byrow = FALSE, ncol = simulation.mixin.quantity + 1)
  
  xmr$Rucknium.Ratio <- xmr$f_S / xmr$Current.decoy.selection.algo.pdf
  
  
  attack.prob <- apply(rings, 1, FUN = function(x) {
    which.max(xmr$Rucknium.Ratio[x])
  })
  
  guessibility.results[[i]] <- t(t(100 * prop.table(table(attack.prob)) ))
  
  print( guessibility.results[[i]])
  
}

guessibility <- sapply( guessibility.results, FUN = function(x) {
  x[nrow(x), ]
})

guessibility <- unname(guessibility)

png("Guessibility-by-ring-size.png", width = 600, height = 600)

plot(seq_along(guessibility) + 1, guessibility, ylim = c(0, max(guessibility)), 
  cex = 0.25,  
  main = "Simulated guessibility as a function of ring size",
  xlab = "Ring size",
  ylab = "Pr(guess real spend correctly | R)")
abline(v = 11, lty = 2, col = "red")
axis(1, 11, "11")

dev.off()

guessibility.data.frame <- data.frame(ring.size = seq_along(guessibility) + 1, guessibility = guessibility)

summary(lm(log(guessibility) ~ ring.size, data = guessibility.data.frame))

write.csv(guessibility.data.frame, 
  file = "Guessibility-by-ring-size-data.csv", row.names = FALSE)

