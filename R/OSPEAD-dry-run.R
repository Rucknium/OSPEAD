
# install.packages("readr")
# install.packages("actuar")
# install.packages("distributionsrd")
# install.packages("future.apply")
# install.packages("RColorBrewer")

xmr.Moser <- as.data.frame(readr::read_csv("data-raw/log_spend_times.csv.xz"))
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
spend_times_blocks <- xmr.Moser$spend_times_blocks
rm(xmr.Moser)


param.trans <- list()

f_D.lgamma <- function(param, support, return.log = FALSE) {
  actuar::dlgamma(support + 1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
}

param.trans$lgamma <- c(exp, exp)


f_D.f <- function(param, support, return.log = FALSE) {
  stats::df(support, df1 = exp(param[1]), df2 = exp(param[2]), ncp = exp(param[3]), log = return.log)
}

param.trans$f <- c(exp, exp, exp)

f_D.rpln <- function(param, support, return.log = FALSE) {
  distributionsrd::drightparetolognormal(support, 
    shape2 = exp(param[1]), meanlog = param[2], sdlog = exp(param[3]), log = return.log)
}

param.trans$rpln <- c(exp, I, exp)

atan.0.1 <- function(x) { atan(x)/pi + 0.5 }

f_D.lgamma.f.mix <- function(param, support, return.log = FALSE) {
  atan.0.1(param[1]) * f_D.lgamma(param[2:3], support, return.log = return.log) + 
    (1 - atan.0.1(param[1])) * f_D.f(param[4:6], support, return.log = return.log)
}

param.trans$lgamma.f.mix <- c(atan.0.1, exp, exp, exp, exp, exp)


f_D.lgamma.cosine <- function(param, support, return.log = FALSE) {
  density.x <- actuar::dlgamma(support + 1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
  
  density.x + atan.0.1(param[3]) * density.x * cos(support * pi * 2 / (24 * 30))
}

param.trans$lgamma.cosine <- c(exp, exp, atan.0.1)

distn.name.converter <- c(
  "f_D.lgamma" = "Log-gamma",
  "f_D.f" = "F",
  "f_D.rpln" = "Right-Pareto Log-normal",
  "f_D.lgamma.f.mix" = "Log-gamma + F Mixture",
  "f_D.lgamma.cosine" = "Log-gamma Cosine"
)


L_FGT <- function(param, f_D, support, flavor = 1) {
  alpha <- flavor
  a_i <- f_D(param, support)
  a_i <- a_i/sum(a_i)
  sum( (theta_i * ((theta_i - a_i)/theta_i)^alpha)[a_i < theta_i & theta_i > 0] , na.rm = FALSE )
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
  a_i <- f_D(param, support)
  a_i <- a_i/sum(a_i)
  (-1) * sum(  
    (theta_i * u_CRRA(ifelse(a_i/theta_i < 1, a_i/theta_i, 1), eta = eta ))[theta_i > 0]
    , na.rm = FALSE )
}

L_MLE <- function(param, f_D, ...) {
  (-1) * sum(
  f_D(param, spend_times_blocks, return.log = TRUE) 
    , na.rm = FALSE)
}


run.iters <- expand.grid(
  f_D = c(f_D.lgamma = f_D.lgamma, f_D.f = f_D.f, f_D.rpln = f_D.rpln), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare, L_MLE = L_MLE)
)

run.iters$flavor[names(run.iters$L) == "L_Welfare"] <- 
  rep(c(0.5, 1), each = sum(names(run.iters$L) == "L_Welfare")/2)

run.iters$flavor[names(run.iters$L) == "L_MLE"] <- 0
# MLE has only one flavor

run.iters <- unique(run.iters)


start.params <- list(
  f_D.lgamma = c(2, 0), 
  f_D.f = c(-1, 0, 2), 
  f_D.rpln = c(0.5, 2, 1))

future::plan(future::multiprocess())

results <- future.apply::future_lapply(1:nrow(run.iters), function(iter) {
  
  f_D.fun <- run.iters$f_D[[iter]]
  L.fun <- run.iters$L[[iter]]
  flavor <- run.iters$flavor[[iter]]
  start.params.optim <- start.params[[names(run.iters$f_D[iter])]]
  
  optim(
    start.params.optim,
    L.fun, 
    f_D = f_D.fun,
    support = as.numeric(support), 
    flavor = flavor,
    control = list(trace = 6, maxit = 10000))

})




run.iters.mix <- expand.grid(
  f_D = c(f_D.lgamma.f.mix = f_D.lgamma.f.mix), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare)
)

run.iters.mix$flavor[names(run.iters.mix$L) == "L_Welfare"] <- 
  rep(c(0.5, 1), each = sum(names(run.iters.mix$L) == "L_Welfare")/2)

# run.iters.mix$flavor[names(run.iters.mix$L) == "L_MLE"] <- 0


results.mix <- future.apply::future_lapply(1:nrow(run.iters.mix), function(iter) {
  
  f_D.fun <- run.iters.mix$f_D[[iter]]
  L.fun <- run.iters.mix$L[[iter]]
  flavor <- run.iters.mix$flavor[[iter]]
  start.params.optim <- c(0,
    results[[which(
      names(run.iters$f_D) == "f_D.lgamma" &
      names(run.iters$L) == names(run.iters.mix$L[iter]) & 
      run.iters$flavor == run.iters.mix$flavor[iter]
        )]]$par,
    results[[which(
      names(run.iters$f_D) == "f_D.f" &
        names(run.iters$L) == names(run.iters.mix$L[iter]) & 
        run.iters$flavor == run.iters.mix$flavor[iter]
    )]]$par
  )
  
  optim(
    start.params.optim,
    L.fun, 
    f_D = f_D.fun,
    support = as.numeric(support), 
    flavor = flavor,
    control = list(trace = 6, maxit = 10000))
  
})



run.iters.cosine <- expand.grid(
  f_D = c(f_D.lgamma.cosine = f_D.lgamma.cosine), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare)
)

run.iters.cosine$flavor[names(run.iters.cosine$L) == "L_Welfare"] <- 
  rep(c(0.5, 1), each = sum(names(run.iters.cosine$L) == "L_Welfare")/2)

# run.iters.cosine$flavor[names(run.iters.cosine$L) == "L_MLE"] <- 0



results.cosine <- future.apply::future_lapply(1:nrow(run.iters.cosine), function(iter) {
  
  f_D.fun <- run.iters.cosine$f_D[[iter]]
  L.fun <- run.iters.cosine$L[[iter]]
  flavor <- run.iters.cosine$flavor[[iter]]
  start.params.optim <- c(
    results[[which(
      names(run.iters$f_D) == "f_D.lgamma" &
        names(run.iters$L) == names(run.iters.cosine$L[iter]) & 
        run.iters$flavor == run.iters.cosine$flavor[iter]
    )]]$par,
    -25
  )
  
  optim(
    start.params.optim,
    L.fun, 
    f_D = f_D.fun,
    support = as.numeric(support), 
    flavor = flavor,
    control = list(trace = 6, maxit = 10000))
  
})





run.iters <- rbind(run.iters, run.iters.mix, run.iters.cosine)

results <- c(results, results.mix, results.cosine)

run.iters.results <- unique(run.iters[, c("L", "flavor")])
run.iters.results$L <- names(run.iters.results$L)

for (j in unique(names(run.iters$f_D))) {
  run.iters.results[, j] <- NA
}


for (i in 1:nrow(run.iters.results)) {
  for (j in unique(names(run.iters$f_D))) {
    
    if (is.na(which(names(run.iters$f_D) == j)[i])) {
      next
    }
    # For distributions that skip certain loss functions
    
    run.iters.results[i, j] <- results[[
      which(names(run.iters$f_D) == j)[i]
    ]]$value
  }
}

families <- colnames(run.iters.results)[grepl("f_D", colnames(run.iters.results))]

for (family in families) {
  run.iters.results[run.iters.results$L == "L_MLE", family] <- 
    2 * length(param.trans[[ gsub("f_D[.]", "", family) ]]) +
    2 * run.iters.results[run.iters.results$L == "L_MLE", family]
}
# AIC formula. Note that the value is the negative of the 
# log likelihood already, so just need to add

run.iters.results




write.csv(run.iters.results, "tables/dry-run/performance.csv", row.names = FALSE)



minimizer.params <- data.frame(f_D = names(run.iters[[1]]), L = names(run.iters[[3]]), flavor = run.iters[[2]],
  param_1 = NA, param_2 = NA, param_3 = NA, stringsAsFactors = FALSE)

dists.to.exclude <- c(
  "f_D.lgamma.f.mix",
  "f_D.lgamma.cosine"
)

minimizer.params <- minimizer.params[ ! minimizer.params$f_D %in% dists.to.exclude, ]

minimizer.params$f_D <- gsub("f_D[.]", "", minimizer.params$f_D)

for ( i in 1:nrow(minimizer.params)) {
  
  temp.param.value <- results[[i]]$par
  
  for ( j in 1:length(temp.param.value)) {
    temp.param.value[j] <- param.trans[[minimizer.params[i, "f_D"]]][[j]](temp.param.value[j])
  }
  
  minimizer.params[i, paste0("param_", seq_along(results[[i]]$par))] <- temp.param.value
    
}

minimizer.params <- minimizer.params[order(minimizer.params$f_D), ]

write.csv(minimizer.params, "tables/dry-run/minimizer-params.csv", row.names = FALSE)


run.iters.obj.fn <- unique(run.iters[, c("L", "flavor")])
run.iters.obj.fn$L <- names(run.iters.obj.fn$L)
run.iters.obj.fn$title <- vector(mode = "list", length = nrow(run.iters.obj.fn))

run.iters.obj.fn$title[run.iters.obj.fn$L == "L_FGT" & run.iters.obj.fn$flavor == 1][[1]] <-
  bquote(L[FGT][alpha] ~ "," ~ alpha == 1)
run.iters.obj.fn$title[run.iters.obj.fn$L == "L_FGT" & run.iters.obj.fn$flavor == 2][[1]] <-
  bquote(L[FGT][alpha] ~ "," ~ alpha == 2)

run.iters.obj.fn$title[run.iters.obj.fn$L == "L_Welfare" & run.iters.obj.fn$flavor == 1][[1]] <-
  bquote(L[Welfare][eta] ~ "," ~ eta == 1)
run.iters.obj.fn$title[run.iters.obj.fn$L == "L_Welfare" & run.iters.obj.fn$flavor == 0.5][[1]] <-
  bquote(L[Welfare][eta] ~ "," ~ eta == 0.5)

run.iters.obj.fn$title[run.iters.obj.fn$L == "L_MLE"][[1]] <-
  bquote(L[MLE])



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
    mat.data <- c(mat.data, a_i[support.interval]/sum(a_i) / theta_i[support.interval] )
    
    legend.labels <- c(legend.labels, 
      paste0(distn.name.converter[ names(run.iters$f_D)[j] ], " | L = ", prettyNum(results[[j]]$value))
    )
  }
  
  png(paste0("images/dry-run/estimate-div-target/estimate-div-target-", loss.fn.temp, "-flavor-", flavor.temp, ".png"), 
    width = 1000 * 1.5, height = 800 * 1.5)
  
  par(mar = c(8, 6, 8, 2) + 0.1, cex = 1.5)
  
  matplot(matrix(mat.data, nrow = length(support.interval)), 
    main = run.iters.obj.fn$title[i][[1]], 
    ylab = bquote(frac(f[D](x), f[S](x))),
    xlab = "x",
    log = "x", xaxt = "n",
    type = "l", lty = 1, 
    ylim = c(min(c(0, mat.data)), min(c(2, max(mat.data)))), 
    col = RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))), "Set1"))
  
  abline(h = 1, lty = 2)
  
  legend("top", 
    legend = legend.labels, 
    lty = rep(1, length(unique(names(run.iters$f_D)))), 
    lwd = rep(5, length(unique(names(run.iters$f_D)))), 
    col = RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))), "Set1"),
    inset = c(0, -0.08), xpd = NA, y.intersp = 1,
    bty = "n", ncol = 3)
  
  for ( j in support.interval) {
    rug(support.interval[j], ticksize = theta_i[j] * 50, col = "blue")
    rug(support.interval[j], ticksize = (-1) * sum(theta_i[1:j] * 0.2), col = "green")
  }
  
  axis(1, at = c(1, 10, 100, 1000, 10000))
  
  axis(1, at = 10000 * 0.8, 
    labels = paste0(round(sum(theta_i[1:max(support.interval)]) * 100), "% cumulative mass"), 
    tick =  FALSE, line = 2)
  
  dev.off()
  
}




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
    mat.data <- c(mat.data, a_i[support.interval]/sum(a_i) )
    
    legend.labels <- c(legend.labels, 
      paste0(distn.name.converter[ names(run.iters$f_D)[j] ], " | L = ", prettyNum(results[[j]]$value))
    )
  }
  
  png(paste0("images/dry-run/estimate/estimate-", loss.fn.temp, "-flavor-", flavor.temp, ".png"), 
    width = 1000 * 1.5, height = 800 * 1.5)
  
  par(mar = c(8, 4, 8, 2) + 0.1, cex = 1.5)
  
  matplot(matrix(mat.data, nrow = length(support.interval)), 
    main = run.iters.obj.fn$title[i][[1]], 
    ylab = bquote(f[D](x)),
    xlab = "x",
    log = "x", xaxt = "n",
    type = "l", lty = 1, 
    ylim = c(min(c(0, mat.data)), min(c(2, max(mat.data)))), 
    col = RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))), "Set1"))
  
  abline(h = 0, lty = 2)
  
  legend("top", 
    legend = legend.labels, 
    lty = rep(1, length(unique(names(run.iters$f_D)))), 
    lwd = rep(5, length(unique(names(run.iters$f_D)))), 
    col = RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))), "Set1"),
    inset = c(0, -0.08), xpd = NA, y.intersp = 1,
    bty = "n", ncol = 3)
  
  for ( j in support.interval) {
    rug(support.interval[j], ticksize = theta_i[j] * 50, col = "blue")
    rug(support.interval[j], ticksize = (-1) * sum(theta_i[1:j] * 0.2), col = "green")
  }
  
  axis(1, at = c(1, 10, 100, 1000, 10000))
  
  axis(1, at = 10000 * 0.8, 
    labels = paste0(round(sum(theta_i[1:max(support.interval)]) * 100), "% cumulative mass"), 
    tick =  FALSE, line = 2)
  
  dev.off()
  
}


