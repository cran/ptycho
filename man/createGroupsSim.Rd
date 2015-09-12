\name{createGroupsSim}
\alias{createGroupsSim}
\title{Create Groups of Covariates}
\description{
  Create an object specifying groups of covariates as needed for some of the
  simulations.
}
\usage{
createGroupsSim(G, p)
}
\arguments{
  \item{G}{Number of groups}
  \item{p}{Number of covariates}
}
\details{
  If \code{G} divides \code{p}, then each group will have \code{p/G} consecutive
  covariates.  If \code{G} does not divide \code{p}, then the last group will
  have fewer covariates.
}
\value{
  List containing the following components:
  \describe{
    \item{\code{var2group}}{Integer vector of length \code{p}, with entry
      \eqn{j} being the index of the group containing covariate \eqn{j}}
    \item{\code{group2var}}{List of length \code{G}, each entry of which is an
      integer vector containing the indices of the covariates belonging to that
      group}
    \item{\code{sizes}}{Vector of length \code{G} containing the number of
      covariates in each group}
  }
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{\code{\link{createData}}}
\examples{
grp <- createGroupsSim(G=3, p=15)
# Which covariates are in group 2?  Two ways to find out:
which(grp$var2group == 2)
grp$group2var[[2]]
}
