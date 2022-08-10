
# install.packages("collapse")
# install.packages("MARSS")
# install.packages("doParallel")
# install.packages("actuar")
# install.packages("moments")
# install.packages(data.table)

library(doParallel)
library(MARSS)
library(data.table)

n.threads <- 50
# Number of CPU threads for parallel processing

blockchain.args <- list(
  BTC  = list(data.dir = "", target.block.time = 600),
  BCH  = list(data.dir = "", target.block.time = 600),
  LTC  = list(data.dir = "", target.block.time = 150),
  DOGE = list(data.dir = "", target.block.time = 60)
)
# Must add path to master_edgelist_output_spent-uncompressed.rds file
# for each blockchain, with trailing slash.


forecast.horizon <- 8
# Number of weeks

years.covered <- "(201[5-9])|(202[0-2])"
# Since 2015

blockchain.results <- list()


for (blockchain in names(blockchain.args)) {
  
  data.dir <- blockchain.args[[blockchain]]$data.dir
  target.block.time <- blockchain.args[[blockchain]]$target.block.time
  
  master.edgelist.output.spent <- readRDS(paste0(data.dir, "master_edgelist_output_spent-uncompressed.rds"))
  
  
  week.set <- levels(master.edgelist.output.spent$output.spent.block_time.week)
  week.set <- week.set[grepl(years.covered, week.set)]
  
  master.edgelist.output.spent <- master.edgelist.output.spent[output.spent.block_time.week %in% week.set, ]
  
  master.edgelist.output.spent[, 
    spend_times_blocks := ifelse(output.spend.age > 0, 
      round(output.spend.age / target.block.time), 0), by = "output.spent.block_time.week"]
  # by = week to do multi-threaded processing
  
  master.edgelist.output.spent[, output.spent.block_time := NULL]
  master.edgelist.output.spent[, output.created.block_height := NULL]
  master.edgelist.output.spent[, output_index := NULL]
  master.edgelist.output.spent[, output.created.block_time := NULL]
  master.edgelist.output.spent[, output.spent.block_height := NULL]
  master.edgelist.output.spent[, output.spend.age := NULL]
  
  
  master.edgelist.output.spent[, output.spent.block_time.week := factor(output.spent.block_time.week)]
  # Remove unused levels
  
  data.table::setDF(master.edgelist.output.spent)
  master.edgelist.output.spent <- split(master.edgelist.output.spent$spend_times_blocks,
    f = master.edgelist.output.spent$output.spent.block_time.week, drop = TRUE)
  
  
  fit.age.distribution <- function(x, f_D, par = NULL, only.eval = FALSE) {
    
    theta_i <- as.data.frame(collapse::qtab(Var1 = x))
    theta_i$Freq <- theta_i$Freq/sum(theta_i$Freq)
    theta_i$Var1 <- as.numeric(as.character(theta_i$Var1))
    theta_i <- merge(data.frame(Var1 = 1:max(theta_i)), theta_i, all.x = TRUE)
    theta_i <- theta_i$Freq
    theta_i[is.na(theta_i)] <- 0
    support <- 1:max(x)
    
    drightparetolognormal.stable <- function(x, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, log = FALSE)  {
      d <- log(shape2 * x^(-shape2 - 1)) +  
        (shape2 * meanlog + (shape2^2 * sdlog^2)/2) + 
        pnorm((log(x) - meanlog - shape2 * sdlog^2)/sdlog, log.p = TRUE)
      if (!log) {
        d <- exp(d)
      }
      return(d)
    }
    # Based on distributionsrd::drightparetolognormal, but 
    # log-linearized to be more numerically stable
    
    f_D.rpln <- function(param, support, return.log = FALSE) {
      drightparetolognormal.stable(support, 
        shape2 = exp(param[1]), meanlog = param[2], sdlog = exp(param[3]), log = return.log)
    }
    
    f_D.lgamma <- function(param, support, return.log = FALSE) {
      support.p1 <- support + 0.001
      actuar::dlgamma(support.p1, shapelog = exp(param[1]), ratelog = exp(param[2]), log = return.log)
    }
    
    L_FGT <- function(param, f_D, support, flavor = 1) {
      alpha <- flavor
      a_i <- f_D(param, support)
      a_i <- a_i/sum(a_i)
      sum( (theta_i * ((theta_i - a_i)/theta_i)^alpha)[a_i < theta_i & theta_i > 0] , na.rm = FALSE )
    }
    
    start.params <- list(
      f_D.lgamma = c(2, 0), 
      f_D.rpln = c(0.5, 2, 1))
    
    stopifnot(f_D %in% c("rpln", "lgamma"))
    
    if (f_D == "rpln") {
      f_D.fun <- f_D.rpln
      if (is.null(par)) {
        start.params.optim  <- distributionsrd::rightparetolognormal.mle(x + 1)$coefficients
        start.params.optim  <- c(log(start.params.optim["shape2"]), 
          start.params.optim["meanlog"], log(start.params.optim["sdlog"]))
        #start.params.optim <- start.params$f_D.rpln
      } else {
        start.params.optim <- par
      }
    }
    
    if (f_D == "lgamma") {
      f_D.fun <- f_D.lgamma
      if (is.null(par)) {
        start.params.optim <- univariateML::mllgamma(x + 1.001)
        start.params.optim <- c(log(start.params.optim["shapelog"]), log(start.params.optim["ratelog"])) 
        # start.params.optim <- start.params$f_D.lgamma
      } else {
        start.params.optim <- par
      }
    }
    
    L.fun <- L_FGT
    flavor <- 1
    
    if (only.eval) {
      return(L.fun(start.params.optim, f_D = f_D.fun, support = as.numeric(support), flavor = flavor))
    }
    
    results <- optim(
      start.params.optim,
      L.fun, 
      method = "BFGS",
      f_D = f_D.fun,
      support = as.numeric(support), 
      flavor = flavor,
      control = list(trace = 6, maxit = 10000))
    
    results
  }
  
  
  cl <- makeCluster(n.threads)
  registerDoParallel(cl)
  
  summary.stats <- data.table::rbindlist( foreach(x = master.edgelist.output.spent) %dopar% { 
    data.table::data.table(mean = mean(x), median = median(x), sd = sd(x),
      skewness = moments::skewness(x), kurtosis = moments::kurtosis(x))
  })
  
  weekly.age <- list()
  
  weekly.age[["rpln"]] <- foreach(x = master.edgelist.output.spent) %dopar% { 
    fit.age.distribution(x, "rpln")
  }
  
  weekly.age[["lgamma"]] <- foreach(x = master.edgelist.output.spent) %dopar% { 
    fit.age.distribution(x, "lgamma")
  }
  
  end.forecast.elements <- (length(master.edgelist.output.spent) - forecast.horizon + 1):length(master.edgelist.output.spent)
  
  naive.forecast.period <- list()
  
  naive.forecast.period[["rpln"]] <- fit.age.distribution(
    unlist(master.edgelist.output.spent[end.forecast.elements]), "rpln")
  
  naive.forecast.period[["lgamma"]] <- fit.age.distribution(
    unlist(master.edgelist.output.spent[end.forecast.elements]), "lgamma")
  
  
  stopCluster(cl)
  
  
  ## Forecast results
  
  forecast.results <- list()
  
  for (fitted.distribution in names(weekly.age)) {
    
    fitted.par <- t(sapply(weekly.age[[fitted.distribution]], function(x) {
      x$par
    }))
    
    rownames(fitted.par) <- names(master.edgelist.output.spent)
    
    fitted.par.trim <- apply(fitted.par, 2, FUN = function(x) {
      lower <- quantile(x, prob = 0.05)
      upper <- quantile(x, prob = 0.95)
      x > (median(x) - 5 * abs(median(x) - lower)) & x < median(x) + 5 * abs(median(x) + upper)
    })
    
    fitted.par <- fitted.par[apply(fitted.par.trim, 1, FUN = all), ]
    
    fitted.par.training <- fitted.par[1:(nrow(fitted.par) - forecast.horizon), ]
    fitted.par.holdout <-  fitted.par[(nrow(fitted.par) - forecast.horizon + 1):nrow(fitted.par), ]
    
    fit <- MARSS(t(fitted.par.training), method = "BFGS", 
      silent = 2, control = list(maxit = 10000))
    # BFGS get convergence faster than EM in this case
    
    
    forecasted.distribution <- forecast.marssMLE(fit, type = "ytT", h = floor(forecast.horizon/2), interval = "none")
    forecasted.distribution <- forecasted.distribution$pred
    forecast.par <- forecasted.distribution[forecasted.distribution$t == max(forecasted.distribution$t), "estimate"]
    forecast.par.naive.final.week <- fitted.par.training[nrow(fitted.par.training), ]
    forecast.par.naive.horizon.period <- naive.forecast.period[[fitted.distribution]]$par
    
    forecast.accuracy <- sapply(master.edgelist.output.spent[rownames(fitted.par.holdout)], 
      FUN = fit.age.distribution, f_D = fitted.distribution, par = forecast.par, only.eval = TRUE)
    # WARNING: This requires fitted.par.holdout to keep its rownames attribute
    
    forecast.accuracy.naive.final.week <- sapply(master.edgelist.output.spent[rownames(fitted.par.holdout)], 
      FUN = fit.age.distribution, f_D = fitted.distribution, par = forecast.par.naive.final.week, only.eval = TRUE)
    
    forecast.accuracy.naive.horizon.period <- sapply(master.edgelist.output.spent[rownames(fitted.par.holdout)], 
      FUN = fit.age.distribution, f_D = fitted.distribution, par = forecast.par.naive.horizon.period, only.eval = TRUE)
    
    forecast.results[[fitted.distribution]] <- 
      list(forecast.accuracy = forecast.accuracy,
        forecast.accuracy.naive.final.week = forecast.accuracy.naive.final.week,
        forecast.accuracy.naive.horizon.period = forecast.accuracy.naive.horizon.period)
    
  }
  
  
  cl <- makeCluster(n.threads)
  registerDoParallel(cl)
  
  empirical.pmf <-  foreach(x = master.edgelist.output.spent) %dopar% { 
    support.viz <- 0:max(x)
    theta_i <- as.data.frame(collapse::qtab(Var1 = x))
    theta_i$Freq <- theta_i$Freq/sum(theta_i$Freq)
    theta_i$Var1 <- as.numeric(as.character(theta_i$Var1))
    theta_i <- merge(data.frame(Var1 = support.viz), theta_i, all.x = TRUE)
    theta_i$Freq[support.viz]
  }
  
  names(empirical.pmf) <- names(master.edgelist.output.spent)
  
  stopCluster(cl)
  
  blockchain.results[[blockchain]] <- list(
    blockchain.args = blockchain.args[[blockchain]],
    week.set = week.set,
    summary.stats = summary.stats,
    weekly.age = weekly.age,
    forecast.results = forecast.results,
    empirical.pmf = empirical.pmf
  )
  
  
}


saveRDS(blockchain.results, file = paste0(blockchain.args$BTC$data.dir, "blockchain-results.rds"))


