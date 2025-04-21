

testthat::test_that("gen.old.new.dsa.rings() produces correct distributions", {
  
  
  dbetabinomial <- function(x, n, alpha, beta) {
    choose(n, x) * beta(x + alpha, n - x + beta) / beta(alpha, beta)
    # Note: beta() is the beta function and beta alone is the parameter value
    # From https://en.wikipedia.org/wiki/Beta-binomial_distribution
  }
  
  
  end.support <- 201
  n <- 16
  C <- 0.4
  beta <- 0.10
  
  real_spend.pmf <- dbetabinomial(0:(end.support - 1), end.support - 1, 0.5, 0.5)
  new_decoy.pmf  <- dbetabinomial(0:(end.support - 1), end.support - 1, 1.5,   2)
  G_star.pmf     <- dbetabinomial(0:(end.support - 1), end.support - 1,   2, 1.5)
  # Note that support of beta-binomial includes 0,
  # but we want support to start at 1. See the draw() functions:
  
  draw.real.spend <- function(n) {
    sample(1:end.support, size = n, replace = TRUE, prob = real_spend.pmf)
  }
  
  draw.decoy.new.dsa <- function(n) {
    sample(1:end.support, size = n, replace = TRUE, prob = new_decoy.pmf)
  }
  
  draw.decoy.old.dsa <- function(n) {
    sample(1:end.support, size = n, replace = TRUE, prob = G_star.pmf)
  }
  
  set.seed(314)
  
  ring.data <- testthat::expect_no_error(
    gen.old.new.dsa.rings(n.rings = 100000, beta = beta,
      old.dsa = draw.decoy.old.dsa,
      new.dsa = draw.decoy.new.dsa,
      ring.size = n,
      change.probability = C,
      real.spend = draw.real.spend, check.data.sorted = TRUE)
  )
  
  
  
  create.f <- function(f_S, f_AD, f_BD, beta, C, n) {
    f_A = (1/n) * f_S + ((n - 1)/n) * f_AD
    f_B = (1/n) * f_S + ((n - 1)/n) * f_BD
    lambda_B = beta + (1/n) * C * (1 - beta)
    f.alpha_B = f_B * lambda_B + f_A * (1 - lambda_B)
    lambda_A = (1 - beta) + (1/n) * C * (beta)
    f.alpha_A = f_A * lambda_A + f_B * (1 - lambda_A)
    f.alpha = (1 - beta) * f.alpha_A + beta * f.alpha_B
    lambda.real.spend = C + (1 - C) * beta
    f.alpha_B.real.spend = f_B * lambda.real.spend + f_A * (1 - lambda.real.spend)
    list(f_A = f_A, f_B = f_B, f.alpha_B = f.alpha_B,
      f.alpha_B.real.spend = f.alpha_B.real.spend,
      f.alpha = f.alpha)
  }
  
  f <- create.f(f_S = real_spend.pmf, f_AD = G_star.pmf,
    f_BD = new_decoy.pmf, beta = beta, C = C, n = n)
  
  
  
  # First, checking if the ring distribution of the tx of interest
  # is the correct distribution
  
  # Chi-sq test is good here. With very large discrete support, KS test may
  # be appropriate.
  
  # "Don't make Type I error", i.e. _should not_ reject null hypothesis
  f_A.type.I.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 0, substr(colnames(X), 1, 41) == "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f_A)
  })
  
  testthat::expect_gt(f_A.type.I.test.results$p.value, sqrt(.Machine$double.eps))
  rm(f_A.type.I.test.results)
  
  # "Don't make Type II error", i.e. _should_ reject null hypothesis
  f_A.type.II.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 0, substr(colnames(X), 1, 41) == "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f_B)
  })
  
  testthat::expect_equal(f_A.type.II.test.results$p.value, 0)
  rm(f_A.type.II.test.results)
  
  # "Don't make Type I error", i.e. _should not_ reject null hypothesis
  f_B.type.I.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 1, substr(colnames(X), 1, 41) == "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f_B)
  })
  
  testthat::expect_gt(f_B.type.I.test.results$p.value, sqrt(.Machine$double.eps))
  rm(f_B.type.I.test.results)
  
  # "Don't make Type II error", i.e. _should_ reject null hypothesis
  f_B.type.II.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 1, substr(colnames(X), 1, 41) == "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f_A)
  })
  
  testthat::expect_equal(f_B.type.II.test.results$p.value, 0)
  rm(f_B.type.II.test.results)
  
  
  
  
  # Get the rings of the antecedent txs that are the real spend
  f_A.type.II.real.spend.antecedents.test.results <- with(ring.data, {
    X <- cbind(y = y, X)
    data <- apply(X[z == 0, ], 1, FUN = function(x) {
      x[substr(colnames(X), 1, 38 + nchar(x["y"])) == paste0("Inputs.0.Previous_Tx_Block_Num_Delta.", x["y"], ".")]
    })
    chisq.test(table(factor(c(data), levels = 1:end.support)),
      p = f$f.alpha_B.real.spend)
  })
  
  
  testthat::expect_equal(f_A.type.II.real.spend.antecedents.test.results$p.value, 0)
  rm(f_A.type.II.real.spend.antecedents.test.results)
  
  
  f_B.type.I.real.spend.antecedents.test.results <- with(ring.data, {
    X <- cbind(y = y, X)
    data <- apply(X[z == 1, ], 1, FUN = function(x) {
      x[substr(colnames(X), 1, 38 + nchar(x["y"])) == paste0("Inputs.0.Previous_Tx_Block_Num_Delta.", x["y"], ".")]
    })
    chisq.test(table(factor(c(data), levels = 1:end.support)),
      p = f$f.alpha_B.real.spend)
  })
  
  testthat::expect_gt(f_B.type.I.real.spend.antecedents.test.results$p.value, sqrt(.Machine$double.eps))
  rm(f_B.type.I.real.spend.antecedents.test.results)
  
  
  
  
  
  # Get the rings of all the antecedent txs
  f_A.type.II.all.antecedents.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 0, substr(colnames(X), 1, 41) != "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f.alpha_B)
  })
  
  testthat::expect_equal(f_A.type.II.all.antecedents.test.results$p.value, 0)
  rm(f_A.type.II.all.antecedents.test.results)
  
  
  f_B.type.I.all.antecedents.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[z == 1, substr(colnames(X), 1, 41) != "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f.alpha_B)
  })
  
  testthat::expect_gt(f_B.type.I.all.antecedents.test.results$p.value, sqrt(.Machine$double.eps))
  rm(f_B.type.I.all.antecedents.test.results)
  
  
  
  
  # Test all antecedent rings
  f.type.I.all.antecedents.test.results <- with(ring.data, {
    chisq.test(table(factor(c(X[, substr(colnames(X), 1, 41) != "Inputs.0.Previous_Tx_Block_Num_Delta.self"]), levels = 1:end.support)),
      p = f$f.alpha)
  })
  
  testthat::expect_gt(f.type.I.all.antecedents.test.results$p.value, sqrt(.Machine$double.eps))
  rm(f.type.I.all.antecedents.test.results)
  
  
  
})

