ptycho.all <- function(data, across=c("none","traits","sites"),
                       doGrpIndicator, dir.out,
                       nreplicates=NULL, ncpu.replicates=1,
                       replicateIterator=ifelse(ncpu.replicates==1,
                                             "replicateLoop","replicateLoopMC"),
                       doSetSeed=TRUE, ncolumns=NULL, ...) {
  across <- match.arg(across)
  if (!file.exists(dir.out)) {
    cat(gettextf("Creating directory %s\n", dir.out))
    dir.create(dir.out, recursive=TRUE)
  }
  if (is.null(nreplicates)) nreplicates <- seq_along(data$replicates)
  do.call(replicateIterator, list(data, across, doGrpIndicator, dir.out,
                                  nreplicates, ncpu.replicates, doSetSeed,
                                  ncolumns, ...))
}

# ncpu.replicates is ignored, but having it eliminates annoyance in ptycho.all
replicateLoop <- function(data, across, doGrpIndicator, dir.out,
                          nreplicates, ncpu.replicates=1,
                          doSetSeed, ncolumns, ...) {
  for (nn in nreplicates) {
    ptycho.replicate(data$X, data$replicates[[nn]], across, doGrpIndicator,
                     dir.out, nn, doSetSeed, ncolumns, ...)
  }
}

# CRAN says I must eliminate the NOTE "no visible binding for global variable"
utils::globalVariables(c("nn"))

replicateLoopMC <- function(data, across, doGrpIndicator, dir.out, nreplicates,
                            ncpu.replicates=1, doSetSeed, ncolumns, ...) {
  # If CRAN were not a bureaucracy, I'd use the following two lines and only
  # Suggest foreach, doMC, and doRNG.
  #require(doMC)
  #doMC::registerDoMC(ncpu.replicates)
  registerDoMC(ncpu.replicates)
  # Assignment to foo prevents samples from being written to stdout
  foo <- foreach (nn=nreplicates) %dopar% {
    ptycho.replicate(data$X, data$replicates[[nn]], across, doGrpIndicator,
                     dir.out, nn, doSetSeed, ncolumns, ...)
  }
}

# ncolumns ignored if across equals "traits"
ptycho.replicate <- function(X, repl, across, doGrpIndicator, dir.out, n.repl,
                             doSetSeed, ncolumns, ...) {
  Y <- repl$y
  q.in <- ncol(Y)
  omtrue <- repl$omega
  if (doSetSeed) set.seed(n.repl)
  if (across == "traits") {
    if (q.in == 1) {
      stop("ptycho.replicates: cannot pool phenotypes when only one input")
    }
    smpl <- ptycho.wrap(X, Y, omega.true=omtrue, doGrpIndicator=doGrpIndicator,
                        dir.out=dir.out, n.repl=n.repl, n.col=1, ...)
  } else {
    q.in <- min(ncolumns, ncol(Y))
    # If I use a*ply or apply here, the object passed into ptycho.wrap is not a
    # matrix and consequently has lost its column name.
    for (nn in seq_len(q.in)) {
      smpl <- ptycho.wrap(X, Y[,nn,drop=FALSE],
                          omega.true=omtrue[,nn,drop=FALSE],
                          doGrpIndicator=doGrpIndicator, dir.out=dir.out,
                          n.repl=n.repl, n.col=nn, ...)
    }
  }
}

# n.repl, n.col are used only for creating output file name
# omega.true is used only for initializing chains
ptycho.wrap <- function(X, y, omega.true=NULL, doGrpIndicator, dir.out, n.repl,
                        n.col, tau.min=0.01, tau.max=10, groups=NULL,
                        ncpu.chain=1, ...) {
  ### I needed to fill in gaps in output.
  #if (file.exists(sprintf("%s/rpl%dcol%d.Rdata", dir.out, n.repl, n.col))) {
  #  cat(paste("Skipping replicate", n.repl, "column", n.col, "\n"))
  #  return(NULL)
  #}
  st0 <- startStates(X, y, doGrpIndicator, omega.true, tau.min, tau.max, groups)
  smpl <- ptycho(X, y, st0, groups, tau.min, tau.max, ncpu=ncpu.chain, ...)
  save(smpl, file=sprintf("%s/rpl%dcol%d.Rdata", dir.out, n.repl, n.col))
  smpl
}

# omega is actual probabilities, not prior on omega
startStates <- function(X, y, doGrpIndicator, omega=NULL, tau.min, tau.max,
                        groups=NULL) {
  nchains <- 4
  tau <- 0.5 * (tau.min + tau.max)
  p <- ncol(X)
  q <- ncol(y)
  # xy[i,j] corresponds to x[,i] and y[,j]
  xy <- -abs(cor(X, y))
  # In each column, largest absolute correlation has rank 1.
  xy <- aaply(xy, 2, rank)
  if (q > 1) xy <- t(xy)
  z <- vector("list", nchains)
  for (nn in seq_len(nchains)) {
    # Initialize indic.var; initialize indic.grp from indic.var at end of
    # function
    if (nn==1) {
      indic.var <- matrix(FALSE, nrow=p, ncol=q)
    } else if (nn==2) {
      indic.var <- matrix(TRUE, nrow=p, ncol=q)
    } else {
      jin <- (if (is.null(omega)) c(10,20)[nn-2]
              else (nn-2)*pmax(1,colSums(omega)))
      if (length(jin) < q) jin <- rep(jin, q)
      indic.var <- (scale(xy, center=jin, scale=FALSE) <= 0)
    }
    if (doGrpIndicator) {
      if (q > 1) {
        indic.grp <- (rowSums(indic.var) > 0)
      } else {
        indic.grp <- laply(groups$group2var, function(x) { sum(indic.var[x]) })
        indic.grp <- (indic.grp > 0)
      }
    } else {
      indic.grp <- NULL
    }
    z[[nn]] <- list(indic.var=indic.var, indic.grp=indic.grp, tau=tau)
  }
  z
}
