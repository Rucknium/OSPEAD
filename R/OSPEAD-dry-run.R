
# install.packages("readr")

xmr.Moser <- as.data.frame(readr::read_csv("data/log_spend_times.csv.xz"))
# Data from Moser et al. (2018) "An Empirical Analysis of Traceability in the Monero Blockchain"

xmr.Moser$spend_times_seconds <- exp(xmr.Moser$log_spend_times)
xmr.Moser$spend_times_blocks <- round(xmr.Moser$spend_times_seconds/120)
xmr.Moser <- xmr.Moser[xmr.Moser$spend_times_blocks > 0, ]

theta_i <- as.data.frame(table(xmr.Moser$spend_times_blocks))
theta_i$Freq <- theta_i$Freq/sum(theta_i$Freq)
theta_i$Var1 <- as.numeric(as.character(theta_i$Var1))
theta_i <- merge(data.frame(Var1 = 1:max(theta_i)), theta_i, all.x = TRUE)
theta_i <- theta_i$Freq
theta_i[is.na(theta_i)] <- 0
support <- 1:max(xmr.Moser$spend_times_blocks)


f_D.gamma <- function(param, support) {
  stats::dgamma(log(support), shape = exp(param[1]), rate = exp(param[2]))
}

f_D.lnorm <- function(param, support) {
  stats::dlnorm(log(support), meanlog = param[1], sdlog = exp(param[2]))
}

f_D.f <- function(param, support) {
  stats::df(log(support), df1 = exp(param[1]), df2 = exp(param[2]), ncp = exp(param[3]) )
}

L_FGT <- function(param, f_D, support, flavor = 1) {
  alpha <- flavor
  a_i <- f_D(param, support)
  a_i <- a_i/sum(a_i)
  sum( (a_i < theta_i & theta_i > 0) * theta_i * ((theta_i - a_i)/theta_i)^alpha , na.rm = TRUE )
}


u_CRRA <- function(x, eta) {
  if (eta != 1) {
    return( (x^(1 - eta)) / (1 - eta)  )
  } else {
    return(log(x))
  }
}

L_Welfare <- function(param, f_D, support, flavor = 1) {
  eta <- flavor
  a_i <- f_D(param, support + 1)
  a_i <- a_i/sum(a_i)
  (-1) * sum(  
    (theta_i > 0) * theta_i * u_CRRA(ifelse(a_i/theta_i < 1, a_i/theta_i, 1), eta = eta )
    , na.rm = TRUE )
}


run.iters <- expand.grid(
  f_D = c(f_D.gamma = f_D.gamma, f_D.lnorm = f_D.lnorm, f_D.f = f_D.f), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare)
)

run.iters$flavor[names(run.iters$L) == "L_Welfare"] <- 
  rep(c(1, 5), each = sum(names(run.iters$L) == "L_Welfare")/2)

start.params <- list(
  f_D.gamma = c(1, 1), 
  f_D.lnorm = c(0, 2), 
  f_D.f = c(2, 2, 2))

results <- list()

for (iter in 1:nrow(run.iters)) {
  
  f_D.fun <- run.iters$f_D[[iter]]
  L.fun <- run.iters$L[[iter]]
  flavor <- run.iters$flavor[[iter]]
  
  results[[iter]] <- optim(
    start.params[[names(run.iters$f_D[iter])]],
    L.fun, 
    f_D = f_D.fun,
    support = support, 
    flavor = flavor,
    control = list(trace = 6))
  
}

run.iters.results <- unique(run.iters[, c("L", "flavor")])
run.iters.results$L <- names(run.iters.results$L)

run.iters.results$f_D.gamma <- NA
run.iters.results$f_D.lnorm <- NA
run.iters.results$f_D.f <- NA

for (i in 1:nrow(run.iters.results)) {
  run.iters.results$f_D.gamma[i] <- results[[
    which(names(run.iters$f_D) == "f_D.gamma")[i]
  ]]$value
  
  run.iters.results$f_D.lnorm[i] <- results[[
    which(names(run.iters$f_D) == "f_D.lnorm")[i]
  ]]$value
  
  run.iters.results$f_D.f[i] <- results[[
    which(names(run.iters$f_D) == "f_D.f")[i]
  ]]$value
}


run.iters.obj.fn <- unique(run.iters[, c("L", "flavor")])
run.iters.obj.fn$L <- names(run.iters.obj.fn$L)
run.iters.obj.fn$title <- vector(mode = "list", length = nrow(run.iters.obj.fn))

run.iters.obj.fn$title[run.iters.obj.fn$L == "L_FGT" & run.iters.obj.fn$flavor == 1][[1]] <-
  bquote(L[FGT][alpha] ~ "," ~ alpha == 1)
run.iters.obj.fn$title[run.iters.obj.fn$L == "L_FGT" & run.iters.obj.fn$flavor == 2][[1]] <-
  bquote(L[FGT][alpha] ~ "," ~ alpha == 2)

run.iters.obj.fn$title[run.iters.obj.fn$L == "L_Welfare" & run.iters.obj.fn$flavor == 1][[1]] <-
  bquote(L[Welfare][eta] ~ "," ~ eta == 1)
run.iters.obj.fn$title[run.iters.obj.fn$L == "L_Welfare" & run.iters.obj.fn$flavor == 5][[1]] <-
  bquote(L[Welfare][eta] ~ "," ~ eta == 5)


distn.name.converter <- c(
  "f_D.gamma" = "Gamma",
  "f_D.lnorm" = "Log-normal",
  "f_D.f" = "F"
)

support.interval <- 1:10000

for (i in 1:nrow(run.iters.obj.fn)) {
  
  loss.fn.temp <- run.iters.obj.fn$L[i]
  flavor.temp <- run.iters.obj.fn$flavor[i]
  
  which.temp <- which(names(run.iters$L) == loss.fn.temp &
      run.iters$flavor == flavor.temp)
  
  mat.data <- c()
  legend.labels <- c()
  
  for (j in which.temp) {
    a_i <- run.iters$f_D[[j]](results[[j]]$par, support)
    mat.data <- c(mat.data, theta_i[support.interval] -  a_i[support.interval]/sum(a_i) )
    
    legend.labels <- c(legend.labels, 
      paste0(distn.name.converter[ names(run.iters$f_D)[j] ], " | L = ", prettyNum(results[[j]]$value))
    )
  }
  
  png(paste0("images/dry-run/", loss.fn.temp, "-flavor-", flavor.temp, ".png"), width = 1000, height = 800)
  
  par(mar = c(8, 4, 5, 2) + 0.1)
  
  matplot(matrix(mat.data, nrow = length(support.interval)), 
    main = run.iters.obj.fn$title[i][[1]], 
    ylab = bquote(f[S](x) - f[D](x)),
    xlab = "x",
    log = "x", xaxt = "n",
    type = "l", lty = 1, col = c("red", "blue", "green"))
  
  abline(h = 0, lty = 2)
  
  legend("top", 
    legend = legend.labels, 
    lty = c(1,1,1), col = c("red", "blue", "green"),
    inset = c(0, -0.04), xpd = NA, y.intersp = 0.5,
    bty = "n", horiz = TRUE)
  
  for ( j in support.interval) {
    rug(support.interval[j], ticksize = theta_i[j] * 50)
    rug(support.interval[j], ticksize = (-1) * sum(theta_i[1:j] * 0.2), col = "green")
  }
  
  axis(1, at = c(1, 10, 100, 1000, 10000))
  
  axis(1, at = 10000 * 0.8, 
    labels = paste0(round(sum(theta_i[1:max(support.interval)]) * 100), "% cumulative mass"), 
    tick =  FALSE, line = 2)
  
  dev.off()
  
}

