
#' Bonhomme-Jochmans-Robin Estimator
#'
#' @param y Data
#' @param II Number of evaluation points for basis function
#' @param K Number of distribution components
#' @param basis Choice of basis function
#' @param cdf.points Points to evaluate the estimated CDF
#' @param estimate.mean.sd Estimate mean and sd?
#' @param estimate.cdf Estimate CDF?
#' @param estimate.mixing.prop Estimate mixing proportions?
#' @param use.C Use C implementation of certain sub-functions?
#' @param debug.return Return function workspace objects before end of function.
#' Supported values include "before est.AA.inner", "after est.AA.inner",
#' "before est.cdf.and.mixing.prop.inner", and "after est.cdf.and.mixing.prop.inner"
#' @param debug.output.to.csv TRUE/FALSE to output objects to csv when
#' debug.return is enabled. A directory above the current working directory
#' called "decoyanalysis-debug-output" will be created to contain the csv files.
#' @param control control list
#'
#' @return The results
#' @export
#'
#' @examples
#' \dontrun{
#' N = 500 # Sample size
#' M = 3   # Number of variables (repeated measurements) 
#' K = 2   # Number of distribution components
#' theta1 = c(0.0, 2.0) # Mean parameters
#' theta2 = c(1.0, 1.0) # Standard deviation parameters
#' omega  = c(0.3, 0.7) # Mixing proportions
#' 
#' set.seed(314)
#' y <- gen.bjr.normal.data(N = N, M = M, K = K, theta1 = theta1,
#'   theta2 = theta2, omega = omega)
#' bjr.results <- bjr(y, II = 10, K = 2)
#' }
bjr <- function(
  y,
  II,
  K,
  basis = "Chebychev",
  cdf.points = NULL,
  estimate.mean.sd = TRUE,
  estimate.cdf = TRUE,
  estimate.mixing.prop = TRUE,
  use.C = FALSE,
  debug.return = "never",
  debug.output.to.csv = FALSE,
  control = list(cluster.threads = NULL)
  ) {
  
  cat(base::date(), "Begin bjr()\n")

  N <- nrow(y)
  M <- ncol(y)
  
  if (is.null(cdf.points)) {
    cdf.points <- quantile(c(y), probs = (0:20)/20, names = FALSE)
  }
  
  
  if (! (basis %in% c("Chebychev", "indicator")) ) {
    stop("basis must be 'Chebychev' or 'indicator'")
  }
  
  ymin = min(y) - diff(range(y)) * 0.001
  ymax = max(y) + diff(range(y)) * 0.001
  ya = (ymin+ymax)/2
  yb = (ymax-ymin)/2
  
  XI = 0:(II+1)
  XI = cos(XI*pi/(II+1))
  XI = yb*XI+ya
  basis.indic = XI[2:(length(XI)-1)]
  rm(XI)
  
  if (basis == "Chebychev") {
    # Chebychev polynomials
    
    g = function(i,y) {
      (2/pi) * cos((i-1) * acos( (y-ya)/yb)) / (2^((i-1)==0) * sqrt(1-((y-ya)/yb)^2))
    }
    
    g.m = function(i,y) {
      (2/pi) * cos((i-1) %*% acos( t((y-ya)/yb))) / (2^((i-1)==0) %*% sqrt(1-t((y-ya)/yb)^2))
    }
    # Matrix version of g()
    
  }
  
  if (basis == "indicator") {
    # indicator functions
    
    g = function(i,y)  {
      y <= basis.indic[i]
    }
    
    g.m = function(i,y) {
      t(y %*% pracma::ones(1, NROW(i)) <= pracma::ones(NROW(y),1) %*% basis.indic[i])
    }
    # Matrix version of g()
    
  }
  
  A <- est.A(y, II, M, N, g.m, basis, ya, yb, basis.indic, control)

  cat(base::date(), "est.A() finished\n")

  # Eigendecomposition of A = VLV' and construction of W
  
  svd.A <- svd(A)
  V <- svd.A$u
  Lambda <- diag(svd.A$d)

  G = sqrt(abs(Lambda))
  G = Re(diag(1 / diag(G[1:K,1:K])))
  W = G %*% t(V[, 1:K])
  
  AA <- est.AA(II, M, g, y, basis, use.C, debug.return, debug.output.to.csv,
    ya, yb, basis.indic, control)
  
  if (debug.return %in% c("before est.AA.inner", "after est.AA.inner")) {
    return(AA)
  }
  
  cat(base::date(), "est.AA() finished\n")
  
  C <- array(dim = c(K, K, II))
  
  # Estimation of joint eigenvectors
  for (j in 1:II) {
    C[, , j] = W %*% AA[, , j] %*% t(W)
  }
  
  frjd.C <- JADE::frjd(C)
  U <- frjd.C$V
  
  if (debug.return == "before est.cdf") {
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv,
      csv.objects = c("cdf.points", "M", "K", "II", "N", "U", "W", "y")))
  }
  
  
  if ((estimate.cdf) & (estimate.mixing.prop)) {
  
    cdf.and.mixing.prop <- est.cdf.and.mixing.prop(cdf.points,
      M, K, II, N, U, W, g, y, basis, use.C, debug.return, debug.output.to.csv,
      ya, yb, basis.indic, control)
    
    if (debug.return %in% c("before est.cdf.and.mixing.prop.inner",
        "after est.cdf.and.mixing.prop.inner")) {
      return(cdf.and.mixing.prop) 
    }
    
    cdf <-  cdf.and.mixing.prop$cdf
    mixing.proportions <- cdf.and.mixing.prop$mixing.proportions
    
    if (!estimate.mean.sd) {
      cat(base::date(), "est.cdf.and.mixing.prop() finished\n")
      cat(base::date(), "End bjr()\n")
      return(list(mu.sigma = NULL, cdf = cdf, mixing.proportions = mixing.proportions))
    }
    
  }
  
  
  if (estimate.cdf) {
    cdf <- est.cdf(cdf.points, M, K, II, N, U, W, g, y)
  } else {
    cdf <- NULL
  }
  
  if (debug.return == "after est.cdf") {
    CDF <- cdf$CDF
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv,
      csv.objects = c("CDF", "cdf.points", "M", "K", "II", "N", "U", "W", "y")))
  }
  
  if ((!estimate.mean.sd) & (!estimate.mixing.prop)) {
    return(list(mu.sigma = NULL, cdf = cdf, mixing.proportions = NULL))
  }
  # If neither estimate.mean.sd nor estimate.mixing.prop, then can stop computations
  
  
  if ((!estimate.mean.sd) & estimate.mixing.prop) {
    mixing.proportions <- est.mixing.prop(M, K, II, N, U, W, g, y)
    return(list(mu.sigma = NULL, cdf = cdf, mixing.proportions = mixing.proportions))
  }
  
  
  if (debug.return == "before est.upsilon") {
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv, csv.objects = c("A", "y", "II", "M", "N")))
  }
  
  upsilon <- est.upsilon(A, y, II, M, N, g)
  
  if (debug.return == "after est.upsilon") {
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv, csv.objects = c("A", "y", "II", "M", "N", "upsilon")))
  }
  
  
  if (debug.return == "before est.Upsilon") {
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv, csv.objects = c("II", "M", "N", "y")))
  }
  
  Upsilon <- est.Upsilon(II, M, N, AA, g, y)
  
  if (debug.return == "after est.Upsilon") {
    return(return.bjr.env(objects = ls(), debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv, csv.objects = c("Upsilon", "II", "M", "N", "y")))
  }
  
  D <- frjd.C$D
  
  for (j in 1:II) {
    D[, , j] = diag(diag(D[, , j]))
  }
  
  LK = Lambda[1:K,1:K]
  
  LL = (kronecker(Lambda, diag(1, K)) - kronecker(diag(1, II), LK))
  
  Jw1 = -kronecker(V,eye(K)) %*% pinv(LL) %*% kronecker(t(V),W) - 
    0.5*KhatriRao(t(W),diag(1,K)) %*% t(KhatriRao(t(W),t(W)))
  LL = (kronecker(LK,diag(1,II))-kronecker(diag(1,K),Lambda ))
  Jw2 =  kronecker(diag(1, K),V) %*% pinv(LL) %*% kronecker(W,t(V))-
    0.5*KhatriRao(diag(1, K),t(W)) %*% t(KhatriRao(t(W),t(W))) 
  
  psi.Ju <- est.psi.Ju(II, K, W, AA, upsilon, Upsilon, Jw1, Jw2, D, U, N)
  psi <- psi.Ju$psi
  Ju <- psi.Ju$Ju
  rm(psi.Ju)
  
  if (estimate.mean.sd & (!estimate.mixing.prop)) {
    mu.sigma <- est.mu.sigma(M, K, II, N, U, Jw2, upsilon, W, Ju, psi, g, y)
    return(list(mu.sigma = mu.sigma, cdf = cdf, mixing.proportions = NULL))
  }
  

  
  if (estimate.mean.sd & estimate.mixing.prop) {
    # est.mu.sigma() and est.mixing.prop() share some inner loop calculations,
    # so computation can be more efficient when combining them if the user chooses
    # to do both
    mu.sigma.and.mixing.prop <- est.mu.sigma.and.mixing.prop(M, K, II, N, U, Jw2, upsilon, W, Ju, psi, g, y)
    mu.sigma <- mu.sigma.and.mixing.prop$mu.sigma
    mixing.proportions <- mu.sigma.and.mixing.prop$mixing.proportions
  } 
  
  return(list(mu.sigma = mu.sigma, cdf = cdf, mixing.proportions = mixing.proportions))
  
}



return.bjr.env <- function(objects, debug.return, debug.output.to.csv, csv.objects) {
  all.objects <- vector("list", length(objects))
  names(all.objects) <- objects
  for (i in objects) {
    all.objects[[i]] <- dynGet(i)
  }
  
  if (debug.output.to.csv) {
    if (! dir.exists("../decoyanalysis-debug-output")) {
      dir.create("../decoyanalysis-debug-output")
    }
    debug.return.name <- gsub("( )|([.])", "-", debug.return)
    debug.dir <- paste0("../decoyanalysis-debug-output/", debug.return.name, "_",
      format(Sys.time(), "%Y-%m-%d-%Hh_%Mm_%Ss"))
    dir.create(debug.dir)
    for (i in csv.objects) {
      x <- dynGet(i)
      if (length(dim(x)) > 0 && length(dim(x)) == 3) {
        x <- apply(x, MARGIN = 2, function(y) {as.matrix(y)})
        # Flatten a 3-dimensional array
      }
      if ( (!is.vector(x) && (!is.matrix(x))) ) {
        stop("objects in csv.objects must be vectors or matrices")
      }
      if (is.vector(x)) {
        cat(x, file = paste0(debug.dir, "/", i, ".txt"), sep = "\n")
      }
      if (is.matrix(x)) {
        write.table(x, file = paste0(debug.dir, "/", i, ".csv"), sep = ",",
          row.names = FALSE, col.names = FALSE)
      }
    }
  }
  
  return(all.objects)
}



nested.future <- function(X, FUN, future.globals = TRUE, control) {
  
  if ( is.null(control$cluster.threads) ) {
    
    n.future.workers <- future::nbrOfWorkers()
    
    if (n.future.workers > 1) {
      X.split <- split(as.data.frame(X), cut(seq_len(nrow(X)), breaks = n.future.workers))
      X.split <- lapply(X.split, function(x) {unname(as.matrix(x))})
      # Convert back to matrix
    } else {
      X.split <- list(X)
    }
    
    rm(X)
    
    return(future.apply::future_lapply(X.split, FUN, future.globals = future.globals))
    
  } else {
    
    X.split.outer <- split(as.data.frame(X), cut(seq_len(nrow(X)),
      breaks = sum(control$cluster.threads)))
    X.split.outer <- lapply(X.split.outer, as.matrix)
    
    X.split.inner <- vector("list", length(control$cluster.threads))
    for (i in seq_along(X.split.inner)) {
      X.split.inner[[i]] <- X.split.outer[ seq_len(control$cluster.threads[i]) ]
      X.split.outer[ seq_len(control$cluster.threads[i]) ] <- NULL
      # Setting elements to NULL (deleting them) resets the index of the next 
      # set of cores to "1". seq_len(control$cluster.threads[i]), which
      # starts from "1", can be used again. It also deletes unneeded data
    }
    
    rm(X.split.outer)
    
    FUN.inner <- FUN
    
    inner.return <- future.apply::future_lapply(X.split.inner, FUN = function(X) {
      return(
        future.apply::future_lapply(X, FUN.inner,
          future.globals = future.globals)
      )
    }, future.globals = c(FUN.inner = FUN.inner, future.globals))
    # Having "FUN.inner" be included in the global list is the same FUN
    # as is being passed to the overall nested.future() function.
    
    inner.return <- unlist(inner.return, recursive = FALSE)
    
    return(inner.return)
    
  }
  
}



est.A <- function(y, II, M, N, g.m, basis, ya, yb, basis.indic, control) {

  ## ESTIMATION OF A and TRANSFORMATION MATRIX W
  # Estimation of two-way contingency table

  i1 <- matrix(rep(1:II, each = II), nrow = II, byrow = TRUE)
  
  if (basis == "Chebychev") {
    mode(i1) <- "numeric"
  }
  # If Chebychev, then g.m() will need to do floating point operations with i1,
  # so convert now. If indicator function as basis, keep as integer
  
  i2 <- t(i1)
  i1 <- c(i1)
  i2 <- c(i2)
  i1.len <- length(i1)
  
  ii <- sort(unique(i1))
  # Will give 1:II vector, but numeric or integer depending
  # on the basis used

  A.future <- function(y) {
    
    N.subset <- nrow(y)
    
    g.m.precompute <- array(0, dim = c(II, N.subset, M))
    
    for (j in 1:M) {
      g.m.precompute[, , j] <- g.m(ii, y[, j])
    }
    # 
    A <- vector("numeric", II^2)
    
    for (j1 in 1:(M-1)) {
      for (j2 in (j1+1):M){
        x.rowSums <- g.m.precompute[i1, , j1] * g.m.precompute[i2, , j2]
        A = A + .Internal(rowSums(x.rowSums, i1.len, N.subset, FALSE))
        # A = A + rowSums(g.m(i1, y[, j1]) * g.m(i2, y[, j2]) )
        # This is sum instead of mean since A is divided by N in
        # this line later:
        # A <- Reduce(`+`, A) / N
      }
    }

    return(A)

  }
  
  A <- nested.future(y, A.future,
    future.globals = list(II = II, M = M, g.m = g.m, i1 = i1, i2 = i2,
      i1.len = i1.len, ii = ii, ya = ya, yb = yb,
      basis.indic = basis.indic), control)

  A <- Reduce(`+`, A) / N
  # Works OK if single-threaded (i.e. number of elements in A is 1)

  A = A*2/M/(M-1)
  A = matrix(A, II, II)
  A = (A+t(A))/2
  A
}




est.AA <- function(II, M, g, y, basis, use.C, debug.return, debug.output.to.csv,
  ya, yb, basis.indic, control) {
  
  ## ESTIMATION OF A1,A2,...,AI AND JOINT EIGENVECTORS U
  # Estimation of three-way contingency tables
  
  Combn.M <- t(combn(M, 3))
  # nrow of Combn.M is number of combinations, not permutations
  # nchoosek(4, 3) = 4
  # nchoosek(11, 3) = 165
  # nchoosek(16, 3) = 560
  Combn.M.nrow <- nrow(Combn.M)
  
  precompute.ind <- unname(as.matrix(expand.grid(1:II, 1:II, 1:II)))
  # Number of rows is II^3
  
  ii <- 1:II
  # Create a vector (ordered set) of integers between 1 and II
  
  if (basis == "Chebychev") {
    mode(ii) <- "numeric"
  }
  # If Chebychev, then g() will need to do floating point operations with ii,
  # so convert now. If indicator function as basis, keep as integer
  
  if (use.C) {
    est.AA.inner <- est.AA.inner.C
  } else {
    est.AA.inner <- est.AA.inner.R
  }
  
  AA.future <- function(y) {
    
    N.subset <- nrow(y)
    
    g.precompute <- array(0, dim = c(II, M, N.subset))
    
    for (i in 1:II) {
      for (j in 1:M) {
        g.precompute[i, j, ] <- g(ii[i], y[, j])
      }
    }
    
    
    if (debug.return == "before est.AA.inner") {
      objects <- c("II", "Combn.M", "Combn.M.nrow", "g.precompute",
        "precompute.ind")
      return(return.bjr.env(objects = objects, debug.return = debug.return,
        debug.output.to.csv = debug.output.to.csv,
        csv.objects = setdiff(objects, "g")))
    }
    
    
    B3 <- est.AA.inner(II, Combn.M, Combn.M.nrow, g.precompute, precompute.ind)
    
    
    if (debug.return == "after est.AA.inner") {
      objects <- c("B3", "II", "Combn.M", "Combn.M.nrow", "g.precompute",
        "precompute.ind")
      return(return.bjr.env(objects = objects, debug.return = debug.return,
        debug.output.to.csv = debug.output.to.csv,
        csv.objects = setdiff(objects, "g")))
    }
    
    aperm.index <- gtools::permutations(3, 3)
    
    AA <- array(0, dim = c(II, II, II))
    
    for (i in 1:nrow(B3)) {
      
      B4 <- array(B3[i, ], dim = c(II, II, II))
      
      B4 <- Reduce(`+`, apply(aperm.index, 1, function(x) {
          aperm(B4, x)
        }, simplify = FALSE))
      
      AA <- AA + B4/nrow(aperm.index)
      
      #B4 = B4 + aperm(B4, c(1, 3, 2)) + aperm(B4, c(2, 1, 3)) +
      #  aperm(B4, c(2, 3, 1)) + aperm(B4, c(3, 1, 2)) + aperm(B4, c(3, 2, 1))
      
      # AA = AA + B4/6
      
    }
    
    return(AA)
    
  }
  
  AA <- nested.future(y, AA.future,
    future.globals = list(est.AA.inner = est.AA.inner, II = II, M = M,
      ii = ii, Combn.M = Combn.M, Combn.M.nrow = Combn.M.nrow,
      precompute.ind = precompute.ind, g = g,
      basis = basis, ya = ya, yb = yb,
      basis.indic = basis.indic, debug.return = debug.return, 
      debug.output.to.csv = debug.output.to.csv), control)
  
  if (debug.return %in% c("before est.AA.inner", "after est.AA.inner")) {
    return(AA[[1]])
  }
  
  AA <- Reduce(`+`, AA)
  
  AA <- AA/nrow(Combn.M)
  
  return(AA)
  
}



est.AA.inner.R <- function(II, Combn.M, Combn.M.nrow, g.precompute, precompute.ind) {
  
  # II is an integer. Will probably be between 10 and 30 in the final version.
  
  # Combn.M is a Combn.M.nrow x 3 matrix of integers. Combn.M.nrow is 560 when
  # ring size is 16.
  
  # Combn.M.nrow is an integer scalar.
  
  # g.precompute is a three-dimensional array of floats with
  # dimensions II x M x N.subset. M is ring size. N.subset is the number of
  # rows (rings) in the subset of the data that has been sent to this
  # particular CPU thread.
  
  # precompute.ind is a matrix of integers with dimensions II^3 x 3. It contains
  # integers between 1 and II.

  B3 <- matrix(0, nrow = Combn.M.nrow, ncol = II^3)
  # Create an empty matrix of floats with dimensions Combn.M.nrow x II^3.
  
  for (q in 1:Combn.M.nrow) {
    
    for (j in 1:(II^3)) {
      B3[q, j] <- .Internal(mean(
        g.precompute[precompute.ind[j, 1], Combn.M[q, 1], ] *
        g.precompute[precompute.ind[j, 2], Combn.M[q, 2], ] *
        g.precompute[precompute.ind[j, 3], Combn.M[q, 3], ]
      ))
    }
    # Get the vector from the g.precompute matrix that is in the position of the
    # _row_ given by the index number in the j'th row, first column of the
    # precompute.ind matrix and the _column_ given by the index number in the
    # q'th row and first column of the Combn.M matrix. The length of
    # this vector will be N.subset. Do the same for the second and third
    # columns of precompute.ind and Combn.M. Then compute the element-wise
    # (Hadamard) product of these three vectors. That will give you a vector
    # with length N.subset. Compute the mean of that vector and assign it to
    # the j'th element of the B3 matrix in its q'th row, j'th column.
    # .Internal() is used for speed boost:
    # https://adv-r.hadley.nz/perf-improve.html#mean
    
  }
  
  return(B3)
  # Return B3
  
}


est.AA.inner.C <- function(II, Combn.M, Combn.M.nrow, g.precompute, precompute.ind) {
  
  # Implement est.AA.inner.R() in C
  
}





est.cdf.and.mixing.prop <- function(cdf.points, M, K, II, N, U, W, g, y, basis,
  use.C, debug.return, debug.output.to.csv, ya, yb, basis.indic, control) {
  # cdf.points is an ordered set of floats. In final version, the length will be 20-100
  # M is an integer. It is the "ring size". So it will be 11 for pre-hardfork and 16 for post-hardfork data
  # K is an integer. It is the number of distribution components to be estimated. It is supposed to 
  # be the number of distinct decoy selection algorithms being used in the wild
  # II is an integer. Will be between 9 and 24 in the final version.
  # N is an integer. It is the number of rings in the dataset. Each week of data has about 200,000 rings.
  # U is a K x K matrix of floats
  # W is a K x II matrix of floats
  # g is the chosen basis function
  # y is a N x M matrix of floats containing the ring member age data
  
  Perms <- gtools::permutations(M, 3)
  # Create a matrix from all three-element permutations of the integers from 1 to M.
  P <- nrow(Perms)
  # P is the number of permutations. 990 when M = 11. 3360 when M = 16.
  
  L = length(cdf.points)
  # L is the number of elements of the set
  
  ii <- 1:II
  # Create a vector (ordered set) of integers between 1 and II
  
  if (basis == "Chebychev") {
    ii <- as.numeric(ii)
  }
  # If Chebychev, then g() will need to do floating point operations with ii, so convert now.
  # If indicator function as basis, keep as integer
  
  if (use.C) {
    est.cdf.and.mixing.prop.inner <- est.cdf.and.mixing.prop.inner.C
  } else {
    est.cdf.and.mixing.prop.inner <- est.cdf.and.mixing.prop.inner.R
  }
  
  CDF.a.hat.B.hat.future <- function(y) {
    
    N.subset <- nrow(y)
    
    W.t <- t(W)
    U.t <- t(U)
    
    g.ind <- unique(Perms[, 1:2])
    ind.3rd <- Perms[, 3]
    
    g.ind.match <- apply(g.ind, 1, paste, collapse = "-")
    Perms.match <- apply(Perms[, 1:2], 1, paste, collapse = "-")
    precompute.ind <- match(Perms.match, g.ind.match)
    
    
    if (debug.return == "before est.cdf.and.mixing.prop.inner") {
      objects <- c("y", "II", "K", "P", "L", "M", "N.subset", "g.ind",
        "precompute.ind", "ind.3rd", "g", "basis", "ii", "W", "W.t", "U", "U.t",
        "cdf.points", "ya", "yb", "basis.indic")
      return(return.bjr.env(objects = objects, debug.return = debug.return,
        debug.output.to.csv = debug.output.to.csv,
        csv.objects = setdiff(objects, "g")))
    }
    
    
    a.hat.B.hat.CDF <- est.cdf.and.mixing.prop.inner(y, II, K, P, L, M, N.subset,
      g.ind, precompute.ind, ind.3rd, g, basis, ii, W, W.t, U, U.t, 
      cdf.points, ya, yb, basis.indic)
    
    a.hat <- a.hat.B.hat.CDF$a.hat
    B.hat <- a.hat.B.hat.CDF$B.hat
    CDF <- a.hat.B.hat.CDF$CDF
    rm(a.hat.B.hat.CDF)
    
    if (debug.return == "after est.cdf.and.mixing.prop.inner") {
      objects <- c("a.hat", "B.hat", "CDF", "y", "II", "K", "P", "L", "M",
        "N.subset", "g.ind", "precompute.ind", "ind.3rd", "g", "basis", "ii",
        "W", "W.t", "U", "U.t", "cdf.points", "ya", "yb", "basis.indic")
      return(return.bjr.env(objects = objects, debug.return = debug.return,
        debug.output.to.csv = debug.output.to.csv,
        csv.objects = setdiff(objects, "g")))
    }
    
    B.hat <- rowSums(B.hat, dims = 2)
    # N.subset x K
    
    return(list(CDF = CDF, B.hat = B.hat, a.hat = a.hat))
    
  }
  
  CDF.a.hat.B.hat <- nested.future(y, CDF.a.hat.B.hat.future,
    future.globals = list(
      est.cdf.and.mixing.prop.inner = est.cdf.and.mixing.prop.inner,
      K = K, M = M, II = II, P = P, L = L, Perms = Perms, g = g,
      cdf.points = cdf.points, basis = basis, ii = ii, W = W, U = U,
      ya = ya, yb = yb, basis.indic = basis.indic, debug.return = debug.return,
      debug.output.to.csv = debug.output.to.csv), control)
  
  if (debug.return %in% c("before est.cdf.and.mixing.prop.inner",
    "after est.cdf.and.mixing.prop.inner")) {
   return(CDF.a.hat.B.hat[[1]]) 
  }
  
  
  CDF <- lapply(CDF.a.hat.B.hat, FUN = function(x) x$CDF)
  CDF <- Reduce(`+`, CDF)
  CDF <- CDF/(as.numeric(P) * as.numeric(N))
  # as.numeric() to prevent possible integer overflow
  
  a.hat <- lapply(CDF.a.hat.B.hat, FUN = function(x) x$a.hat)
  a.hat <- Reduce(`+`, a.hat)
  a.hat <- a.hat/(as.numeric(M) * as.numeric(N))
  a.hat <- matrix(a.hat, nrow = II, ncol = 1)
  
  B.hat <- lapply(CDF.a.hat.B.hat, FUN = function(x) x$B.hat)
  rm(CDF.a.hat.B.hat)
  B.hat <- Reduce(`+`, B.hat)
  B.hat <- B.hat/(as.numeric(P) * as.numeric(N))
  
  mixing.proportions <- solve(t(B.hat) %*% B.hat) %*% t(B.hat) %*% a.hat
  
  list(cdf = list(CDF = t(CDF), cdf.points = cdf.points),
    mixing.proportions = c(mixing.proportions))
  # Return the CDF matrix, cdf.points vector, and mixing proportions vector,
  # organized in a list
  
}



est.cdf.and.mixing.prop.inner.R <- function(
    y, II, K, P, L, M, N.subset, g.ind, precompute.ind, ind.3rd, g,
    basis, ii, W, W.t, U, U.t, cdf.points, ya, yb, basis.indic) {
  
  # y is a N.subset x M matrix of floats containing the ring member age data.
  
  # II is an integer. Will probably be between 10 and 30 in the final version.
  
  # K is an integer. It is the number of distribution components to be
  # estimated. It is supposed to be the number of distinct decoy selection
  # algorithms being used in the wild.
  
  # P is an integer. It is the number of three-element permutations of the
  # integers from 1 to M. 990 when M = 11. 3360 when M = 16.
  
  # L an integer. It is the length of the cdf.points vector. It will probably
  # be between 100 and 1,000 in practice. It could go higher.
  
  # M is an integer. It is the number of repeated measurements. In our setting,
  # it is the ring size.
  
  # N.subset is an integer. N is the number of rings in the dataset (200,000 in 
  # a typical week). N.subset is the size of the subset of the data that has
  # been sent to this particular CPU thread.
  
  # g.ind is a two-column matrix of integers containing every two-element
  # permutation of the integers from 1 to M. When M = 11, the number of
  # rows is 110. When M = 16, the number of rows is 240
  
  # precompute.ind is a one-dimensional vector of integers. Its length
  # is P. It contains every integer between 1 and the number of rows of g.ind
  
  # ind.3rd is a one-dimensional vector of integers. Its length
  # is P. It contains every integer between 1 and M.
  
  # g is the basis function implemented in R. It won't be used by C.
  
  # basis is a character object equal to "Chebychev" or "indicator". This can
  # be used by the C code in an if statement to switch between basis functions.
  
  # ii is a one-dimensional vector with length II. The elements of the vector
  # will be floats if basis == "Chebychev" and integers if basis == "indicator".
  
  # U is a K x K matrix of floats. U.t is its matrix transpose.
  
  # W is a K x II matrix of floats. W.t is its matrix transpose.
  
  # cdf.points is a one-dimensional vector of floats with length L.
  
  # ya and yb are float scalars used by the Chebychev basis function. Because
  # of R's scoping rules, they do not need to be explicitly passed to g()
  
  # basis.indic is a one-dimensional vector of floats with length II used
  # by the indicator basis function.
  
  
  a.hat <- vector("numeric", II)
  # Create an empty one-dimensional vector of floats with length II
  
  B.hat <- array(0, dim = c(II, K, P))
  # Create an empty three-dimensional array of floats with dimensions II, K, and P
  
  CDF <- matrix(0, nrow = K, ncol = L)
  # Create an empty matrix of floats with K rows and L columns
  
  for (n in 1:N.subset) {
    
    g.precompute <- matrix(0, nrow = II, ncol = M)
    # Create an empty matrix of floats with II rows and M columns
    
    for (m in 1:M) {
      g.precompute[, m] <- g(ii, y[n, m])
    }
    # Compute the value of the basis function with ii as its first argument
    # and the element og y in its n'th row and m'th column. Assign this
    # value to the m'th column of g.precompute.
    # ii is a one-dimensional vector of integers or floats rather than a scalar. 
    # R deals with this case seamlessly, but C may need a loop.
    
    a.hat <- a.hat + rowSums(g.precompute)
    # Start with the II x M matrix g.precompute. Compute the sum of each of its
    # rows to create a one-dimensional vector with length II. Add the vector
    # to a.hat
    
    for (k in 1:K) {
      
      k.U   <-   U[, k, drop = FALSE]
      k.U.t <- U.t[k, , drop = FALSE]
      # Get the k'th column of the U matrix and k'th row of the U.t matrix.
      # "drop = FALSE" keeps the object as a two-dimensional matrix instead of
      # turning it into a one-dimensional vector so that matrix multiplication
      # can be performed later
      
      tau.precompute <- vector("numeric", nrow(g.ind))
      # Create an empty one-dimensional vector of floats with length equal
      # to the number of rows of g.ind
      
      for (mm in 1:nrow(g.ind)) {
        
        g.center <- outer(
          g.precompute[, g.ind[mm, 1] ],
          g.precompute[, g.ind[mm, 2] ]
        )
        # Get the element of g.ind in its mm'th row, 1st column to get the
        # desired column index of g.precompute. Do the same for the mm'th row,
        # 2nd column of g.ind. Then compute the matrix outer product of these
        # two vector. Since the first and second object both have length II, the
        # dimensions of the resulting matrix are II x II
        
        tau.precompute[mm] <- k.U.t %*% W %*% g.center %*% W.t %*% k.U
        # This is a series of matrix multiplications (not element-wise
        # multiplications). The %*% is the matrix multiplication operator
        # The result of these multiplications is a scalar, which is
        # assigned to the mm'th element of tau.precompute.
        
      }
      
      tauP <- tau.precompute[precompute.ind]
      # precompute.ind is a vector of integers intended to be an indexer. It
      # contains some duplicate index number in a specific order. when it
      # is used as an index for tau.precompute, the elements of tau.precompute
      # that correspond to the integers in precompute.ind are assigned to
      # tauP in a specific order and with (appropriately) duplicated elements.
      # The tauP result is a one-dimensional vector with length P.
      
      basis.3rd <- g.precompute[, ind.3rd]
      # Like the line above, ind.3rd serves as an indexer. In this case, it
      # operates on a matrix instead of a vector. It is only operating on
      # the column indices of g.precompute. It takes the whole column for every
      # index number in ind.3rd. The matrix produced by this indexing
      # selection, basis.3rd, is a II x P matrix
      
      tauP.expanded <- matrix(tauP, nrow = II, ncol = P, byrow = TRUE)
      # Produces a II x P matrix whose rows are each a duplicate of
      # the tauP vector. 
      
      B.hat[, k, ] <- B.hat[, k, ] + tauP.expanded * basis.3rd
      # tauP.expanded and basis.3rd are both II x P matrices. The k'th "column"
      # of the B.hat three-dimensional array is a II x P matrix. This operation
      # computes the element-wise (Hadamard) product of the tauP.expanded and
      # basis.3rd matrices and adds them to the matrix defined by the
      # k'th "column" of B.hat.
      
      for (l in 1:L) {
        CDF[k, l] <- CDF[k, l] + sum(tauP * (y[n, ind.3rd] <= cdf.points[l]))
      }
      # y[n, ind.3rd] gets the n'th row of y and the columns of y corresponding
      # to the elements of the ind.3rd indexer. ind.3rd has length P.
      # cdf.points[l] produces a sclar by getting the l'th element of the
      # cdf.points vector. The "<=" operator between them produces a logical
      # vector with length P. Each element of the vector is TRUE when the left
      # side is less than or equal to the right side. FALSE otherwise.
      # The logical vector is multiplied (element-wise) by tauP. The logical
      # vector is converted to a numeric vector when multiplied. TRUE = 1
      # and FALSE = 0. The vector result of this multiplication is summed
      # to produce a scalar. It is added to the element in the CDF matrix
      # in the k'th row and l'th column.
      
    }
  }
  
  return(list(a.hat = a.hat, B.hat = B.hat, CDF = CDF))
  # Return the a.hat, B.hat, and CDF objects as three objects in a list
  
}



est.cdf.and.mixing.prop.inner.C <- function(
  y, II, K, P, L, M, N.subset, g.ind, precompute.ind, ind.3rd, g,
  basis, ii, W, W.t, U, U.t, cdf.points, ya, yb, basis.indic) {
  
  # Implement est.cdf.and.mixing.prop.inner.R() in C
    
    
}





######################################################
######################################################
# IGNORE BELOW THIS LINE
######################################################
######################################################





est.upsilon <- function(A, y, II, M, N, g) {
  # A is a II x II matrix of floats
  # y is a N x M matrix of floats containing the ring member age data
  # II is an integer. Will be between 9 and 24 in the final version
  # M is an integer. It is the "ring size". So it will be 11 for pre-hardfork and 16 for post-hardfork data
  # N is an integer. It is the number of rings in the dataset. Each week of data has about 200,000 rings.
  # g is the chosen basis function
  
  X <- array(0, dim = c(II, II, N))
  # Create an empty II x II x N array
  
  xi <- matrix(0, nrow = II^2, ncol = N)
  # Create an exmpty II^2 x N matrix
  
  perm_expr <- gtools::permutations(M, 2)
  # Create a matrix from all two-element permutations of the integers from 1 to M.
  
  perm_expr_expanded <- str2expression(paste0("g(i1,y[n,", perm_expr[, 1], "]) * g(i2,y[n,",
    perm_expr[, 2], "])", collapse = " + "))
  # Create the expression to be evaluated
  
  # Number of inner iterations: N*II^2
  for (n in 1:N) {
    for (i1 in 1:II) {
      for (i2 in 1:II) {
        X[i1,i2,n] = eval(perm_expr_expanded)
        # Assign the result of eval(perm_expr_expanded) to the element of X in
        # the i1th, i2th, nth position of the array
      }
    }
    xi[, n] = c(X[, , n])/(M^2 - M)
    # Take the elements of X in the nth position of its 3rd dimension, divide by (M^2 - M),
    # and assign to the nth row of xi
  }
  upsilon = xi - matrix(A, nrow = II^2, ncol = 1) %*% matrix(1, nrow = 1, ncol = N)
  # This can be done outside of the C code
  # Form a II^2 x 1 matrix from the contents of the A matrix. Take its outer
  # product with a 1 x N matrix, then subtract it from the xi matrix to get upsilon
  
  upsilon
  # Return upsilon
  
}




est.Upsilon <- function(II, M, N, AA, g, y) {
  # II is an integer. Will probably be between 9 and 24 in the final version.
  # M is an integer. It is the "ring size". It will be 11 for pre-hardfork and 16 for post-hardfork data.
  # N is an integer. It is the number of rings in the dataset. Each week of data has about 200,000 rings.
  # AA is an 3-dimensional array of floats with dimensions II x II x II
  # g is the chosen basis function
  # y is a N x M matrix of floats containing the ring member age data
  
  perm_expr <- gtools::permutations(M, 3)
  # Create a matrix from all three-element permutations of the integers from 1 to M.
  
  perm_expr_expanded <- str2expression(paste0("g(i1,y[n,", perm_expr[, 1], "]) * g(i2,y[n,",
    perm_expr[, 2], "]) * g(i3,y[n,", perm_expr[, 3], "])", collapse = " + "))
  # Concatenate the expression (take the product of each triple and then sum the products)
  
  Xi <- matrix(0, nrow = II^3, ncol = N)
  # Create an empty matrix with II^3 rows and N columns
  
  # There will be N*II^3 inner iterations
  for (n in 1:N) {
    X = array(0, dim = c(II,II,II))
    # Create an empty 3-dimensional array with dimensions II x II x II 
    for (i1 in 1:II) {
      for (i2 in 1:II) {
        for (i3 in 1:II) {
          X[i1,i2,i3] = eval(perm_expr_expanded)
          # Assign the result of the evaluated expression to the i1'th, i2'th, i3'th element of the X matrix
        }
      }
    }
    Xi[, n] = c(X)/nrow(perm_expr)
    # Fill the n'th row of Xi with the "flattened" X array, divided by the 
    # total number of permutations in perm_expr
  }
  
  Upsilon = Xi - matrix(AA, II^3, 1) %*% matrix(1, ncol = N)
  # This can be done outside of the C code
  # Form a II^3 x 1 matrix from the contents of the AA matrix. Take its outer
  # product with a 1 x N matrix, then subtract it from the Xi matrix to get Upsilon
  
  Upsilon
  # Return Upsilon
  
}





est.psi.Ju <- function(II, K, W, AA, upsilon, Upsilon, Jw1, Jw2, D, U, N) {
  
  BWL = array(dim = c(K*K, K*II, II))
  BBL = matrix(0, K*K,K*II)
  for (i in 1:II) {
    BWL[, , i] = kronecker(W %*% AA[, , i], diag(1, K))
    BBL = rbind(BBL, BWL[ , , i])
  }
  BWL = BBL
  BWL <- BWL[(-1) * 1:(K*K), ]
  
  BWR = array(dim = c(K*K, K*II, II))
  BBR = matrix(0, K*K,K*II)
  for (i in 1:II) {
    BWR[ , , i] = kronecker(diag(1, K), W %*% AA[, , i])
    BBR = rbind(BBR, BWR[ , , i])
  }
  BWR = BBR
  BWR <- BWR[(-1) * 1:(K*K), ]
  
  
  psi = BWL %*% Jw1 %*% upsilon+BWR %*% Jw2 %*% upsilon+kronecker(diag(1, II),kronecker(W,W)) %*% Upsilon
  
  DD = matrix(0,K*K,K*K)
  RR = matrix(0,K*K,K*K) 
  for (i in 1:II) {
    DD = DD + (kronecker(D[ , , i],diag(1,K))-kronecker(diag(1,K),D[ , , i]))^2
    RR = cbind(RR, (kronecker(D[ , , i],diag(1,K))-kronecker(diag(1,K),D[ , , i]))) 
  }
  
  RR <- RR[, (-1) * 1:(K*K)]
  
  Ju = -kronecker(diag(1, K),U) %*% pinv(DD) %*% RR %*% kronecker(diag(1, II),kronecker(t(U),t(U)))
  
  list(psi = psi, Ju = Ju)
  
}




est.mu.sigma <- function(M, K, II, N, U, Jw2, upsilon, W, Ju, psi, g, y) {
  
  # ESTIMATION OF WEIGHT FUNCTION AND MOMENTS
  Perms = gtools::permutations(M, 3)
  
  mu = matrix(0,K,1)
  dmu = matrix(0,K,II)
  sq = matrix(0,K,1)
  dsq = matrix(0,K,II)
  iota = kronecker(t(U),diag(1,II)) %*% Jw2 %*% upsilon+kronecker(diag(1,K),t(W)) %*% Ju %*% psi
  EYE = diag(1,K)
  imu = array(dim = c(N, nrow(Perms), K))
  isq = array(dim = c(N, nrow(Perms), K))
  sd = vector(length = K)
  
  inf_mu <- matrix(nrow = N, ncol = K)
  se_mu <- vector(length = K)
  inf_sq <- matrix(nrow = N, ncol = K)
  se_sq <- vector(length = K)
  inf_va <- matrix(nrow = N, ncol = K)
  se_va <- vector(length = K)
  
  ii <- 1:II
  
  for (k in 1:K) {
    for (n in 1:N) {
      for (m in 1:nrow(Perms)) {
        i1 = Perms[m,1]
        i2 = Perms[m,2]
        i3 = Perms[m,3]
        Xn = g(ii,y[n,i1]) %*% t(g(ii,y[n,i2])) 
        
        meann = y[n,i3]
        squan = meann^2
        taun = t(U[, k]) %*% W %*% Xn %*% t(W) %*% U[, k]
        dtaun = 2*t(U[, k]) %*% W %*% Xn
        
        mu[k] = mu[k] + taun %*% meann/(nrow(Perms)*N)
        dmu[k, ] = dmu[k, ] + dtaun*meann/(nrow(Perms)*N) # first  moment 
        
        sq[k] = sq[k] + taun %*% squan/(nrow(Perms)*N)
        dsq[k, ] = dsq[k, ] + dtaun*squan/(nrow(Perms)*N) # second moment
        
        imu[n,m,k] = taun %*% meann
        isq[n,m,k] = taun %*% squan
        
      } 
    } 
    
    sd[k] = sqrt(sq[k]-mu[k]^2) # standard deviation
    
    inf_mu[, k] = as.matrix((rowMeans(imu[, , k]) - mu[k])+ t(dmu[k, ] %*% kronecker(EYE[k, , drop = FALSE],diag(1, II)) %*% iota))
    se_mu[k] = sqrt(diag(t(inf_mu[, k]) %*% inf_mu[, k])/N^2)
    inf_sq[, k] = as.matrix((rowMeans(isq[, , k]) - sq[k])+ t(dsq[k, ] %*% kronecker(EYE[k, , drop = FALSE],eye(II)) %*% iota))
    se_sq[k] = sqrt(diag(t(inf_sq[, k]) %*% inf_sq[, k])/N^2)
    inf_va[, k] = inf_sq[, k]-.5*mu[k] %*% inf_mu[, k]
    se_va[k] = sqrt(diag(t(inf_va[, k]) %*% inf_va[, k])/N^2) # influence function for variance 
  }
  
  se_sd = (.5/sd)*se_va
  
  list(mu = mu, se_mu = se_mu, sd = sd, se_sd = se_sd)
  
}







est.cdf <- function(cdf.points, M, K, II, N, U, W, g, y) {
  # cdf.points is an ordered set of floats. In final version, the length will be 20-100
  # M is an integer. It is the "ring size". So it will be 11 for pre-hardfork and 16 for post-hardfork data
  # K is an integer. It is the number of distribution components to be estimated. It is supposed to 
  # be the number of distinct decoy selection algorithms being used in the wild
  # II is an integer. Will be between 9 and 24 in the final version.
  # N is an integer. It is the number of rings in the dataset. Each week of data has about 200,000 rings.
  # U is a K x K matrix of floats
  # W is a K x II matrix of floats
  # g is the chosen basis function
  # y is a N x M matrix of floats containing the ring member age data
  
  Perms <- gtools::permutations(M, 3)
  # Create a matrix from all three-element permutations of the integers from 1 to M.
  M <- nrow(Perms)
  # M is the number of permutations. 990 when M = 11. 3360 when M = 16.
  
  L = length(cdf.points)
  # L is the number of elements of the set
  
  CDF = matrix(0, nrow = K, ncol = L)
  # Create a K x L matrix
  
  ii <- 1:II
  # Create an ordered set of integers between 1 and II
  
  # Total number of inner iterations: K*L*N*M
  for (k in 1:K) {
    for (l in 1:L) {
      for (n in 1:N) {
        for (m in 1:M) {
          
          i3 = Perms[m, 3]
          # Get the 3rd index number from the Perms matrix
          
          if (! y[n,i3] <= cdf.points[l] ) {
            next
          }
          # If the i3'th ring member of the n'th ring is not below the
          # l'th cdf.support.point, the loop doesn't need to calculate anything
          
          i1 = Perms[m,1]
          i2 = Perms[m,2]
          # Get 1st and 2nd index number from the Perms matrix
          
          Xn = outer(g(ii,y[n,i1]), g(ii,y[n,i2])) 
          # ii is an ordered set of integers of length II rather than a scalar. 
          # R deals with this case seamlessly, but C may need a loop.
          # outer() computers the matrix outer product of these two objects
          # Since the first and second object both have length II, the
          # dimensions of the resulting matrix are II x II
          
          
          taun = t(U[, k]) %*% W %*% Xn %*% t(W) %*% U[, k]
          # This is a series of matrix multiplications (not element-wise multiplications)
          # %*% is the matrix multiplication operator
          # t() creates the transpose of the matrix.
          # U[, k] gets the k'th column of the U matrix
          
          CDF[k,l] = CDF[k,l] + taun/(M*N)
          # Add taun divided by M*N to the element in the k'th row, l'th column
          # of the CDF matrix
          
        } 
      }
      cat(base::date(), "    CDF:  K: ", k, "/", K, "  L: ", l, "/", L, "\n", sep = "")
    }
  }
  
  list(CDF = t(CDF), cdf.points = cdf.points)
  # Return the CDF matrix and cdf.points ordered set,
  # organized in a list
  
}






est.mixing.prop <- function(M, K, II, N, U, W, g, y) {
  
  Perms = gtools::permutations(M, 3)
  
  B.hat = matrix(0, II, K)
  
  ii <- 1:II
  
  for (k in 1:K) {
    for (n in 1:N) {
      for (m in 1:nrow(Perms)) {
        i1 = Perms[m,1]
        i2 = Perms[m,2]
        i3 = Perms[m,3]
        Xn = g(ii, y[n, i1]) %*% t(g(ii, y[n, i2])) 
        
        basis.3rd = g(ii, y[n, i3])
        
        taun = t(U[, k]) %*% W %*% Xn %*% t(W) %*% U[, k]
        
        B.hat[,k] = B.hat[,k] + c(taun) * basis.3rd
        
      } 
    } 
  }
  
  B.hat = B.hat/(nrow(Perms)*N)
  
  a.hat <- matrix(0,II,1)
  
  for (n in 1:N) {
    for (m in 1:M) {
      a.hat = a.hat + g(ii, y[n, m])
    }
  }
  
  a.hat = a.hat/(M*N)
  
  mixing.proportions <- solve(t(B.hat) %*% B.hat) %*% t(B.hat) %*% a.hat
  
  mixing.proportions
  
}












est.mu.sigma.and.mixing.prop <- function(M, K, II, N, U, Jw2, upsilon, W, Ju, psi, g, y) {
  
  # ESTIMATION OF WEIGHT FUNCTION AND MOMENTS
  Perms = gtools::permutations(M, 3)
  
  mu = matrix(0,K,1)
  dmu = matrix(0,K,II)
  sq = matrix(0,K,1)
  dsq = matrix(0,K,II)
  iota = kronecker(t(U),diag(1,II)) %*% Jw2 %*% upsilon+kronecker(diag(1,K),t(W)) %*% Ju %*% psi
  EYE = diag(1,K)
  imu = array(dim = c(N, nrow(Perms), K))
  isq = array(dim = c(N, nrow(Perms), K))
  sd = vector(length = K)
  
  inf_mu <- matrix(nrow = N, ncol = K)
  se_mu <- vector(length = K)
  inf_sq <- matrix(nrow = N, ncol = K)
  se_sq <- vector(length = K)
  inf_va <- matrix(nrow = N, ncol = K)
  se_va <- vector(length = K)
  
  B.hat = matrix(0, II, K)
  
  ii <- 1:II
  
  for (k in 1:K) {
    for (n in 1:N) {
      for (m in 1:nrow(Perms)) {
        i1 = Perms[m,1]
        i2 = Perms[m,2]
        i3 = Perms[m,3]
        Xn = g(ii,y[n,i1]) %*% t(g(ii,y[n,i2])) 
        
        meann = y[n,i3]
        squan = meann^2
        taun = t(U[, k]) %*% W %*% Xn %*% t(W) %*% U[, k]
        dtaun = 2*t(U[, k]) %*% W %*% Xn
        
        mu[k] = mu[k] + taun %*% meann/(nrow(Perms)*N)
        dmu[k, ] = dmu[k, ] + dtaun*meann/(nrow(Perms)*N) # first  moment 
        
        sq[k] = sq[k] + taun %*% squan/(nrow(Perms)*N)
        dsq[k, ] = dsq[k, ] + dtaun*squan/(nrow(Perms)*N) # second moment
        
        imu[n,m,k] = taun %*% meann
        isq[n,m,k] = taun %*% squan
        
        basis.3rd = g(ii, y[n, i3])
        
        B.hat[,k] = B.hat[,k] + c(taun) * basis.3rd
        
      } 
    } 
    
    sd[k] = sqrt(sq[k]-mu[k]^2) # standard deviation
    
    inf_mu[, k] = as.matrix((rowMeans(imu[, , k]) - mu[k])+ t(dmu[k, ] %*% kronecker(EYE[k, , drop = FALSE],diag(1, II)) %*% iota))
    se_mu[k] = sqrt(diag(t(inf_mu[, k]) %*% inf_mu[, k])/N^2)
    inf_sq[, k] = as.matrix((rowMeans(isq[, , k]) - sq[k])+ t(dsq[k, ] %*% kronecker(EYE[k, , drop = FALSE],eye(II)) %*% iota))
    se_sq[k] = sqrt(diag(t(inf_sq[, k]) %*% inf_sq[, k])/N^2)
    inf_va[, k] = inf_sq[, k]-.5*mu[k] %*% inf_mu[, k]
    se_va[k] = sqrt(diag(t(inf_va[, k]) %*% inf_va[, k])/N^2) # influence function for variance 
  }
  
  se_sd = (.5/sd)*se_va
  
  B.hat = B.hat/(nrow(Perms)*N)
  
  a.hat <- matrix(0,II,1)
  
  for (n in 1:N) {
    for (m in 1:M) {
      a.hat = a.hat + g(ii, y[n, m])
    }
  }
  
  a.hat = a.hat/(M*N)
  
  mixing.proportions <- solve(t(B.hat) %*% B.hat) %*% t(B.hat) %*% a.hat
  
  list(mu.sigma = list(mu = mu, se_mu = se_mu, sd = sd, se_sd = se_sd),
    mixing.proportions = mixing.proportions)
  
}







