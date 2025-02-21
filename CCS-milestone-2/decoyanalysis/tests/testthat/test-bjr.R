
skip.R.tests <- FALSE



testthat::test_that("C implementation of est.AA.inner() matches R implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip("C implementation of est.AA.inner() not ready")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation.R <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    debug.return = "after est.AA.inner")
  
  bjr.validation.C <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    use.C = TRUE, debug.return = "after est.AA.inner")
  
  testthat::expect_equal(bjr.validation.C$B4, bjr.validation.R$B4,
    tolerance = 1e-14)
  
})





testthat::test_that("C implementation of est.cdf.and.mixing.prop.inner() matches R implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip("C implementation of est.cdf.and.mixing.prop.inner() not ready")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())

  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation.R <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    debug.return = "after est.cdf.and.mixing.prop.inner")
  
  bjr.validation.C <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    use.C = TRUE, debug.return = "after est.cdf.and.mixing.prop.inner")
  
  testthat::expect_equal(bjr.validation.C$B4, bjr.validation.R$B4,
    tolerance = 1e-14)
  
})



testthat::test_that("Overall C implementation of parts of bjr() matches R implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip("Overall C implementation of parts of bjr() not ready")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation.R <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  bjr.validation.C <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    use.C = TRUE)
  
  testthat::expect_equal(bjr.validation.C$B4, bjr.validation.R$B4,
    tolerance = 1e-14)
  
})






testthat::test_that("bjr() matches Octave implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  testthat::expect_equal(bjr.validation$cdf$CDF, ValidationResultsOctaveCDFChebychev500N,
    tolerance = 1e-14)
  
})


testthat::test_that("Two-thread bjr() matches Octave implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  stopifnot(parallelly::availableCores() >= 2)
  future::plan(future::multisession(workers = 2))
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  testthat::expect_equal(bjr.validation$cdf$CDF, ValidationResultsOctaveCDFChebychev500N,
    tolerance = 1e-14)
  
})


testthat::test_that("Nested-thread bjr() matches Octave implementation (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  
  testthat::skip_if_not(parallelly::availableCores() >= 4,
    "Need at least 4 CPU threads for this test")
 
  future::plan(list(
    future::tweak(future::multisession, workers = 2),
    future::tweak(future::multisession, workers = 2)
    ))
  # The main purpose of the nested threads is to use cluster computing
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    control = list(cluster.threads = c(2, 2)))
  
  testthat::expect_equal(bjr.validation$cdf$CDF, ValidationResultsOctaveCDFChebychev500N,
    tolerance = 1e-14)
  
})




testthat::test_that("bjr() matches Octave implementation (Monte Carlo, indicator as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "indicator",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  testthat::expect_equal(bjr.validation$cdf$CDF, ValidationResultsOctaveCDFindicator500N,
    tolerance = 1e-14)

})


testthat::test_that("bjr() matches Octave implementation (Waterdata, indicator as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  y <- Waterdata[, c(10, 5, 6, 7)]
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  ymin = min(y) - diff(range(y)) * 0.001
  ymax = max(y) + diff(range(y)) * 0.001
  ya = (ymin+ymax)/2
  yb = (ymax-ymin)/2

  L = 121
  XIk = matrix(0:(L-1), ncol = 1)
  XI = cos(XIk*pi/(L-1))
  XIXI = c(yb*XI+ya)
  cdf.points = XIXI

  bjr.validation <- bjr(y, II = 7, K = 3, basis = "indicator",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  testthat::expect_equal(bjr.validation$cdf$CDF, ValidationResultsOctaveCDFindicatorWaterdata,
    tolerance = 1e-09)
  
})




testthat::test_that("bjr() correctly calculates mixing proportions (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  mixing.proportions.validation <- c(0.700996345913642, 0.296700701485614)
  
  testthat::expect_equal(bjr.validation$mixing.proportions, mixing.proportions.validation,
    tolerance = 1e-14)
  
})





testthat::test_that("Two-thread bjr() correctly calculates mixing proportions (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  stopifnot(parallelly::availableCores() > 1)
  future::plan(future::multisession(workers = 2))
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  mixing.proportions.validation <- c(0.700996345913642, 0.296700701485614)
  
  testthat::expect_equal(bjr.validation$mixing.proportions, mixing.proportions.validation,
    tolerance = 1e-14)
  
})




testthat::test_that("Nested-thread bjr() correctly calculates mixing proportions (Monte Carlo, Chebychev as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  
  testthat::skip_if_not(parallelly::availableCores() >= 4,
    "Need at least 4 CPU threads for this test")
  
  future::plan(list(
    future::tweak(future::multisession, workers = 2),
    future::tweak(future::multisession, workers = 2)
  ))
  # The main purpose of the nested threads is to use cluster computing
  
  cdf.points = seq(-4, 7, by = .25)
  
  withr::with_preserve_seed({ y <- gen.standard.bjr.test.dataset() })
  # Don't affect the random seed in the R global environment
  K <- 2
  
  bjr.validation <- bjr(y, II = 10, K = K, basis = "Chebychev",
    cdf.points = cdf.points, estimate.mean.sd = FALSE,
    control = list(cluster.threads = c(2, 2)))
  
  mixing.proportions.validation <- c(0.700996345913642, 0.296700701485614)
  
  testthat::expect_equal(bjr.validation$mixing.proportions, mixing.proportions.validation,
    tolerance = 1e-06)
  
})





testthat::test_that("bjr() correctly calculates mixing proportions (Waterdata, indicator as basis)", {
  
  testthat::skip_if_not( ! skip.R.tests, "R-specific tests skipped")
  
  y <- Waterdata[, c(10, 5, 6, 7)]
  
  restore.plan <- future::plan()
  on.exit(future::plan(restore.plan))
  # Save the multi-threaded future plan that the user has set so it can be restored later
  future::plan(future::sequential())
  
  ymin = min(y) - diff(range(y)) * 0.001
  ymax = max(y) + diff(range(y)) * 0.001
  ya = (ymin+ymax)/2
  yb = (ymax-ymin)/2
  
  L = 121
  XIk = matrix(0:(L-1), ncol = 1)
  XI = cos(XIk*pi/(L-1))
  XIXI = c(yb*XI+ya)
  cdf.points = XIXI
  
  bjr.validation <- bjr(y, II = 7, K = 3, basis = "indicator",
    cdf.points = cdf.points, estimate.mean.sd = FALSE)
  
  mixing.proportions.validation <- c(0.575646080680339, 0.233299211784955, 0.210189455424116)
  
  testthat::expect_equal(bjr.validation$mixing.proportions, mixing.proportions.validation,
    tolerance = 1e-13)
  
})



