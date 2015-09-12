\name{mcmcLoop}
\alias{mcmcLoop}
\title{MCMC Sampling Loop}
\description{Loop to perform MCMC sampling for one chain.}
\usage{
mcmcLoop(data, params, state)
}
\arguments{
  \item{data}{List containing the observed data}
  \item{params}{List containing model and computational parameters}
  \item{state}{List containing initial state of the MCMC iteration}
}
\details{
  The only reason this is a user level function is so that users can use
  different packages for parallelization by writing their own version of the
  function that calls \code{mcmcLoop}.  The name of the custom function should
  be the value of argument \code{chainIterator} to \code{\link{ptycho}}.
  Such functions should be similar to the following function that is included in
  the package:
  \preformatted{
    chainLoopRNG <- function(data, params, initStates, ncpu) {
      require(doMC)
      require(doRNG)
      registerDoMC(ncpu)
      nchains <- length(initStates)
      z <- foreach (nn=1:nchains) \%dorng\% {
        state <- initStates[[nn]]
        mcmcLoop(data, params, state)
      }
      z
    }
  }
  The package \code{doRNG} is used in this code to enable reproducibility.
}
\value{
  An \code{\link{mcmc}} object if \code{params$only.means} is \code{NULL}.
  Otherwise, a matrix.
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{\code{\link{ptycho}}}
