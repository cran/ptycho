\name{PosteriorStatistics}
\alias{PosteriorStatistics}
\alias{meanTau}
\alias{varTau}
\alias{meanIndicators}
\alias{meanVarIndicators}
\alias{meanGrpIndicators}
\title{Extract Posterior Statistics}
\description{
  Extract posterior statistics from a \code{\link{ptycho}} object.
}
\usage{
meanTau(obj)
varTau(obj)
meanIndicators(obj)
meanVarIndicators(obj)
meanGrpIndicators(obj)
}
\arguments{
  \item{obj}{A \code{\link{ptycho}} object}
}
\details{
  A \code{\link{ptycho}} object contains means for many different variables.  If
  multiple chains were run, it has separate means for each, and it may have
  running means from different points within each chain.  The functions
  described here simplify extracting from the input object certain statistics of
  the posterior distribution sampled by \code{\link{ptycho}}.
  
  The function \code{meanTau} identifies the last iteration saved in the input
  \code{\link{ptycho}} object and computes the mean of \eqn{\tau} at that
  iteration across all chains.  The function \code{varTau} is analogous,
  computing var\eqn{(\tau)}.

  Similarly, \code{meanIndicators} returns the mean across all chains of each
  indicator variable.  The functions \code{meanVarIndicators} and
  \code{meanGrpIndicators} compute the means only of the indicators of variants
  or only of second-level indicator variables, respectively.
}
\value{
  Both \code{meanTau} and \code{varTau} return a scalar.

  The other functions, which extract means of indicator variables, return
  vectors with names copied from the column names of the input \code{obj}.
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{
  \code{\link{ptycho}}, \link{WhichCols}; also \link{Data} describes
  \code{ptychoIn} and \code{ptychoOut} in example below
}
\examples{
data(ptychoIn)
data(ptychoOut)
# Compare averages of sampled group indicator variables to truth.
cbind(ptychoIn$replicates[[1]]$indic.grp,
      meanGrpIndicators(ptychoOut))
# Compare averages of sampled covariate indicator variables to truth.
cbind(ptychoIn$replicates[[1]]$indic.var,
      meanVarIndicators(ptychoOut))
# Compare average of sampled values of tau to truth.
ptychoIn$replicates[[1]]$tau
meanTau(ptychoOut)
# Variance of sampled values of tau is reasonable because sampled model
# is usually NOT empty.
varTau(ptychoOut)
}
