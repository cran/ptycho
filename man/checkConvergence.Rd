\name{checkConvergence}
\alias{checkConvergence}
\title{Compute Differences Between MCMC Chains}
\description{
  Compute the differences between the chains in a \code{\link{ptycho}} object
  of the means of \eqn{\tau} and the indicator variables.
}
\usage{
checkConvergence(obj, doLastIterOnly=TRUE)
}
\arguments{
  \item{obj}{A \code{\link{ptycho}} object}
  \item{doLastIterOnly}{Logical specifying whether to compute differences only
    for the last MCMC iteration in the input object or for all iterations}
}
\details{
  For \eqn{\tau} and each indicator variable in the input \code{\link{ptycho}}
  object, compute the difference between the maximum and the minimum mean for
  each chain.  If \code{doLastIterOnly} is \code{TRUE}, then the differences are
  only computed for the last iteration in each chain; otherwise, the differences
  are computed at each iteration in the input object.
}
\value{
  A data frame with the following columns:
  \describe{
    \item{\code{iter}}{MCMC iteration number}
    \item{\code{type}}{Factor specifying the type of the variable; one of
      \dQuote{tau}, \dQuote{var} for variant indicator variable, or \dQuote{grp}
      for second-level indicator variable}
    \item{\code{index}}{Number specifying the pertinent column in the design
      matrix (for \code{type} equal to \dQuote{var} or for \code{type} equal to
      \dQuote{grp} when \emph{Across Traits} prior was used) or the variant
      group index (for \code{type} equal to \dQuote{grp} when \emph{Across
        Sites} prior was used);
      equal to 1 for \code{type} equal to \dQuote{tau}}
    \item{\code{y}}{Factor specifying the name of the response; empty for
      \code{type} equal to \dQuote{tau} or for \code{type} equal to \dQuote{grp}
      when \emph{Across Traits} prior was used}
    \item{\code{range}}{Difference between maximum and minimum across chains}
  }
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{
  \code{\link{ptycho}}; also \code{\link{ptychoOut}} for example below
}
\examples{
data(ptychoOut)
cvg <- checkConvergence(ptychoOut, doLastIterOnly=FALSE)
reshape2::dcast(cvg, ... ~ iter, value.var="range")
}
