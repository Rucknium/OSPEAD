
# install.packages("readr")
# install.packages("actuar")
# install.packages("distributionsrd")
# install.packages("future.apply")
# install.packages("RColorBrewer")
# install.packages("VGAM")
# install.packages("huxtable")
library(huxtable)

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
  support.p1 <- support + 1
  actuar::dlgamma(support.p1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
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

f_D.gev <- function(param, support, return.log = FALSE) {
  VGAM::dgev(support, location = param[1], scale = exp(param[2]), shape = param[3], log = return.log)
}

param.trans$gev <- c(I, exp, I)


atan.0.1 <- function(x) { atan(x)/pi + 0.5 }


f_D.lgamma.f.mix <- function(param, support, return.log = FALSE) {

  atan.0.1(param[1]) * f_D.lgamma(param[2:3], support, return.log = return.log) + 
    (1 - atan.0.1(param[1])) * f_D.f(param[4:6], support, return.log = return.log)
}

param.trans$lgamma.f.mix <- c(atan.0.1, exp, exp, exp, exp, exp)

f_D.gev.f.mix <- function(param, support, return.log = FALSE) {

  atan.0.1(param[1]) * f_D.gev(param[2:4], support, return.log = return.log) + 
    (1 - atan.0.1(param[1])) * f_D.f(param[5:7], support, return.log = return.log)
}

param.trans$f_D.gev.f.mix <- c(atan.0.1, I, exp, I, exp, exp, exp)

# DONT USE COSINE:

# f_D.lgamma.cosine <- function(param, support, return.log = FALSE) {
#   density.x <- actuar::dlgamma(support + 1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
  
#   density.x + atan.0.1(param[3]) * density.x * cos(support * pi * 2 / (24 * 30))
# }

# param.trans$lgamma.cosine <- c(exp, exp, atan.0.1)




f_D.lgamma.periodic.laplace <- function(param, support, return.log = FALSE) {
  
  laplace.multiplier <- 
    (2 * 24 * 30) / sum(VGAM::dlaplace((support + 1) %% (2 *24 * 30), location = 24 * 30, scale = exp(param[4]), log = return.log)[1:(2 * 24 * 30)])
  
  density.x <- actuar::dlgamma(support + 1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
  
  (1 - atan.0.1(param[3])) * density.x + atan.0.1(param[3]) * density.x * laplace.multiplier * 
    VGAM::dlaplace((support + 1) %% (2 *24 * 30), location = 24 * 30, scale = exp(param[4]), log = return.log)
}

# TODO: The above actually has a period of two days since it fits better with the Moser data. 1-day period would be:

# f_D.lgamma.periodic.laplace <- function(param, support, return.log = FALSE) {
  
#   laplace.multiplier <- 
#     (24 * 30) / sum(VGAM::dlaplace((support + 1 + (24 * 30 / 2)) %% (24 * 30), location = 24 * 30 / 2, scale = param[4], log = return.log)[1:(24 * 30)])
  
#   density.x <- actuar::dlgamma(support + 1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
  
#   (1 - atan.0.1(param[3])) * density.x + atan.0.1(param[3]) * density.x * laplace.multiplier * 
#     VGAM::dlaplace((support + 1 + (24 * 30 / 2)) %% (24 * 30), location = 24 * 30 / 2, scale = param[4], log = return.log)
    # }

param.trans$f_D.lgamma.periodic.laplace <- c(exp, exp, atan.0.1, exp)

distn.name.converter <- c(
  "f_D.lgamma" = "Log-gamma",
  "f_D.f" = "F",
  "f_D.rpln" = "Right-Pareto Log-normal",
  "f_D.gev" = "Generalized Extreme Value",
  "f_D.lgamma.f.mix" = "Log-gamma + F Mixture",
  "f_D.gev.f.mix" = "Generalized Extreme Value + F Mixture",
  "f_D.lgamma.periodic.laplace" = "Log-gamma + Periodic Laplace"
  #"f_D.lgamma.cosine" = "Log-gamma Cosine"
)


L_FGT <- function(param, f_D, support, flavor = 1, theta_i, ...) {
  alpha <- flavor
  a_i <- f_D(param, support)
  a_i <- a_i/sum(a_i)
  sum( (theta_i * ((theta_i - a_i)/theta_i)^alpha)[a_i < theta_i & theta_i > 0] , na.rm = FALSE )
}



L_Welfare <- function(param, f_D, support, flavor = 1, theta_i, ...) {
  
  u_CRRA <- function(x, eta) {
    if (eta != 1) {
      return( (x^(1 - eta)) / (1 - eta)  )
    } else {
      return(log(x))
    }
  }
  
  eta <- flavor
  a_i <- f_D(param, support)
  a_i <- a_i/sum(a_i)
  (-1) * sum(  
    (theta_i * u_CRRA(ifelse(a_i/theta_i < 1, a_i/theta_i, 1), eta = eta ))[theta_i > 0]
    , na.rm = FALSE )
}

L_MLE <- function(param, f_D, spend_times_blocks = spend_times_blocks, ...) {
  (-1) * sum(
  f_D(param, spend_times_blocks, return.log = TRUE) 
    , na.rm = FALSE)
}


run.iters.simple <- expand.grid(
  f_D = c(f_D.lgamma = f_D.lgamma, f_D.f = f_D.f, f_D.rpln = f_D.rpln, f_D.gev = f_D.gev), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare, L_MLE = L_MLE)
)

run.iters.simple$flavor[names(run.iters.simple$L) == "L_Welfare"] <- 
  rep(c(0.5, 1), each = sum(names(run.iters.simple$L) == "L_Welfare")/2)

run.iters.simple$flavor[names(run.iters.simple$L) == "L_MLE"] <- 0
# MLE has only one flavor

run.iters.simple <- unique(run.iters.simple)


start.params <- list(
  f_D.lgamma = c(2, 0), 
  f_D.f = c(-1, 0, 2), 
  f_D.rpln = c(0.5, 2, 1),
  f_D.gev = c(100, log(100), 1))


future::plan(future::multiprocess())

results.simple <- future.apply::future_lapply(1:nrow(run.iters.simple), function(iter) {
  
  f_D.fun <- run.iters.simple$f_D[[iter]]
  L.fun <- run.iters.simple$L[[iter]]
  flavor <- run.iters.simple$flavor[[iter]]
  start.params.optim <- start.params[[names(run.iters.simple$f_D[iter])]]
  
  optim(
    start.params.optim,
    L.fun, 
    f_D = f_D.fun,
    support = as.numeric(support), 
    flavor = flavor,
    control = list(trace = 6, maxit = 10000),
    theta_i = theta_i,
    spend_times_blocks = spend_times_blocks)

})



run.iters.mix <- expand.grid(
  f_D = c(
    f_D.lgamma.f.mix = f_D.lgamma.f.mix,
    f_D.gev.f.mix = f_D.gev.f.mix), 
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
  
  mixture.name <- names(run.iters.mix$f_D[iter])
  mixture.name <- gsub("(f_D.)|(.mix)", "", mixture.name)
  mixture.name <- strsplit(mixture.name, "[.]")[[1]]
  # NOTE: Naming is crucial and order matters
  
  start.params.optim <- c(0,
    results.simple[[which(
      names(run.iters$f_D) == paste0("f_D.", mixture.name[1]) &
      names(run.iters$L) == names(run.iters.mix$L[iter]) & 
      run.iters$flavor == run.iters.mix$flavor[iter]
        )]]$par,
    results.simple[[which(
      names(run.iters$f_D) == paste0("f_D.", mixture.name[2]) &
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
    control = list(trace = 6, maxit = 10000),
    theta_i = theta_i,
    spend_times_blocks = spend_times_blocks,
    atan.0.1 = atan.0.1,
    f_D.lgamma = f_D.lgamma,
    f_D.gev = f_D.gev,
    f_D.f = f_D.f)
  
})



# run.iters.cosine <- expand.grid(
#   f_D = c(f_D.lgamma.cosine = f_D.lgamma.cosine), 
#   flavor = 1:2,
#   L = c(L_FGT = L_FGT, L_Welfare = L_Welfare)
# )

run.iters.periodic <- expand.grid(
  f_D = c(f_D.lgamma.periodic.laplace = f_D.lgamma.periodic.laplace), 
  flavor = 1:2,
  L = c(L_FGT = L_FGT, L_Welfare = L_Welfare)
)

run.iters.periodic$flavor[names(run.iters.periodic$L) == "L_Welfare"] <- 
  rep(c(0.5, 1), each = sum(names(run.iters.periodic$L) == "L_Welfare")/2)

# run.iters.periodic$flavor[names(run.iters.periodic$L) == "L_MLE"] <- 0

results.periodic <- future.apply::future_lapply(1:nrow(run.iters.periodic), function(iter) {
  
  f_D.fun <- run.iters.periodic$f_D[[iter]]
  L.fun <- run.iters.periodic$L[[iter]]
  flavor <- run.iters.periodic$flavor[[iter]]
  start.params.optim <- c(
    results.simple[[which(
      names(run.iters$f_D) == "f_D.lgamma" &
        names(run.iters$L) == names(run.iters.periodic$L[iter]) & 
        run.iters$flavor == run.iters.periodic$flavor[iter]
    )]]$par,
    -25, # atan.0.1-transformed mixing factor
    log(24 * 30 / 2) # Scale parameter of the Laplace distribution
  )
  
  optim(
    start.params.optim,
    L.fun, 
    f_D = f_D.fun,
    support = as.numeric(support), 
    flavor = flavor,
    control = list(trace = 6, maxit = 10000),
    theta_i = theta_i,
    spend_times_blocks = spend_times_blocks,
    atan.0.1 = atan.0.1)
  
})



run.iters <- rbind(run.iters.simple, run.iters.mix, run.iters.periodic)

results <- c(results.simple, results.mix, results.periodic)

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

hux.run.iters.results <- run.iters.results

hux.run.iters.results$flavor[hux.run.iters.results$flavor == 0] <- NA

colnames(hux.run.iters.results) <- c("Loss function", "Loss function parameter", "Log-gamma", "F", "Right-Pareto Log-normal", "Generalized Extreme Value", "Log-gamma + F mix", "Log-gamma + GEV mix", "Log-gamma + Laplace Periodic")

hux.run.iters.results <- huxtable::as_hux(hux.run.iters.results)
hux.run.iters.results <- t(hux.run.iters.results)
hux.run.iters.results <- huxtable::set_bottom_border(hux.run.iters.results, row = 2, col = huxtable::everywhere)
hux.run.iters.results <- huxtable::set_align(hux.run.iters.results, row = huxtable::everywhere, col = 1, value = "left")
hux.run.iters.results <- huxtable::set_number_format(hux.run.iters.results, row = 3:9, col = 2:5, value = 4)

YlGn.colors <- rev(RColorBrewer::brewer.pal(7, "YlGn"))

for (iter in 1:5) {
  if (iter == 5) {
    hux.run.iters.results <- huxtable::set_background_color(hux.run.iters.results, 
      col = 1 + iter, row = order(unlist(run.iters.results[iter, 3:6])) + 2, value = YlGn.colors[1:4])
  } else {
    hux.run.iters.results <- huxtable::set_background_color(hux.run.iters.results, 
      col = 1 + iter, row = order(unlist(run.iters.results[iter, 3:9])) + 2, value = YlGn.colors)
  }
}

hux.run.iters.results <- huxtable::add_footnote(hux.run.iters.results, 
  text = "Note: Values should be compared down columns. Lower values (darker green) indicate better performance. MLE value is Akaike Information Criterion (AIC).")

hux.run.iters.results <- huxtable::set_wrap(hux.run.iters.results, col = 1, row = 1:(nrow(hux.run.iters.results) - 1), value = FALSE)

width(hux.run.iters.results) <- 1
# Makes the cells wrap

hux.run.iters.results <- huxtable::set_col_width(hux.run.iters.results, col = 2:6, value = (1/ncol(hux.run.iters.results)) * 0.8)

cat(huxtable::to_latex(hux.run.iters.results, tabular_only = TRUE), file = "tables/dry-run/performance.tex")



minimizer.params <- data.frame(f_D = names(run.iters[[1]]), L = names(run.iters[[3]]), flavor = run.iters[[2]],
  param_1 = NA, param_2 = NA, param_3 = NA, stringsAsFactors = FALSE)

dists.to.exclude <- c(
  "f_D.lgamma.f.mix",
  "f_D.gev.f.mix",
  "f_D.lgamma.periodic.laplace"
  #"f_D.lgamma.cosine"
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

minimizer.params

write.csv(minimizer.params, "tables/dry-run/minimizer-params.csv", row.names = FALSE)


hux.minimizer.params <- minimizer.params

hux.minimizer.params$f_D <- distn.name.converter[match(hux.minimizer.params$f_D, 
  gsub("f_D.", "", names(distn.name.converter)))]

colnames(hux.minimizer.params) <- c("Distribution", "Loss fn", "Loss fn param", "param_1", "param_2", "param_3")

hux.minimizer.params <- huxtable::as_hux(hux.minimizer.params)

hux.minimizer.params <- huxtable::set_wrap(hux.minimizer.params, col = 1, row = 1:(nrow(hux.minimizer.params)), value = FALSE)

width(hux.minimizer.params) <- 0.95

hux.minimizer.params <- huxtable::add_footnote(hux.minimizer.params, 
  text = "F Distribution: param_1 is first degree of freedom parameter; param_2 is second degree of freedom parameter; param_3 is non-centrality parameter.")

hux.minimizer.params <- huxtable::add_footnote(hux.minimizer.params, 
  text = "Generalized Extreme Value Distribution: param_1 is location parameter; param_2 is scale parameter; param_3 is shape parameter.")

hux.minimizer.params <- huxtable::add_footnote(hux.minimizer.params, 
  text = "Log-gamma Distribution: param_1 is shape parameter; param_2 is rate parameter.")

hux.minimizer.params <- huxtable::add_footnote(hux.minimizer.params, 
  text = "Right-Pareto Log-normal Distribution: param_1 is shape parameter; param_2 is mean parameter; param_3 is variance parameter.")

hux.minimizer.params <- huxtable::add_footnote(hux.minimizer.params, 
  text = "Mixture distributions are omitted from this table.")


cat(huxtable::to_latex(hux.minimizer.params, tabular_only = TRUE), file = "tables/dry-run/minimizer-params.tex")



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
    width = 1000 * 1.5, height = 1200 * 1.5)
  
  par(mar = c(10, 6, 10, 2) + 0.1, cex = 1.5)
  
  matplot(matrix(mat.data, nrow = length(support.interval)), 
    main = run.iters.obj.fn$title[i][[1]], 
    ylab = bquote(frac(f[D](x), f[S](x))),
    xlab = "Age of spent outputs in terms of number of blocks. Log scale.                                                            ",
    log = "x", xaxt = "n", 
    type = "l", lty = 1, 
    ylim = c(min(c(0, mat.data)), min(c(2, max(mat.data)))), 
    col = (RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))) + 1, "Set1")[-6]))
  
  abline(h = 1, lty = 2)
  
  legend("top", 
    legend = legend.labels, 
    lty = rep(1, length(unique(names(run.iters$f_D)))), 
    lwd = rep(10, length(unique(names(run.iters$f_D)))), 
    col = (RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))) + 1, "Set1")[-6]),
    inset = c(0, -0.07), xpd = NA, y.intersp = 1,
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
    width = 1000 * 1.5, height = 1200 * 1.5)
  
  par(mar = c(10, 6, 10, 2) + 0.1, cex = 1.5)
  
  matplot(matrix(mat.data, nrow = length(support.interval)), 
    main = run.iters.obj.fn$title[i][[1]], 
    ylab = bquote(f[D](x)),
    xlab = "Age of spent outputs in terms of number of blocks. Log scale.                                                            ",
    log = "x", xaxt = "n",
    type = "l", lty = 1, 
    ylim = c(min(c(0, mat.data)), min(c(2, max(mat.data)))), 
    col = (RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))) + 1, "Set1")[-6]))
  
  abline(h = 0, lty = 2)
  
  legend("top", 
    legend = legend.labels, 
    lty = rep(1, length(unique(names(run.iters$f_D)))), 
    lwd = rep(10, length(unique(names(run.iters$f_D)))), 
    col = (RColorBrewer::brewer.pal( length(unique(names(run.iters$f_D))) + 1, "Set1")[-6]),
    inset = c(0, -0.07), xpd = NA, y.intersp = 1,
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


