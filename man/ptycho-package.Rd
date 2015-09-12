\name{ptycho-package}
\alias{ptycho-package}
\docType{package}
\title{Bayesian Variable Selection with Hierarchical Priors}
\description{
  Bayesian variable selection for linear regression models using hierarchical
  priors. There is a prior that combines information across responses and one
  that combines information across covariates, as well as a standard spike and
  slab prior for comparison. An MCMC samples from the marginal posterior
  distribution for the 0-1 variables indicating if each covariate belongs to the
  model for each response.
}
\details{
This package provides functions to carry out Bayesian model selection combining
different layers of information: across multiple traits or across multiple
variants in the same gene.  The priors are described by Stell and Sabatti
(2015).  To sample the posterior distribution for specified genotype and
phenotype matrices, use \code{\link{ptycho}}.

This package also provides functions to generate simulated data as in Stell and
Sabatti (2015); see \code{\link{createData}} and
\url{web.stanford.edu/~lstell/ptycho/}.  Those datasets are not included in this
package because they have images about 20 MB or larger.  Instead small data
objects are included for examples; see \link{Data}.

Functions for post-processing \code{ptycho} objects are described at
\code{\link{checkConvergence}} and \link{PosteriorStatistics}.
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\references{
  Stell, L. and Sabatti, C. (2015) Genetic variant selection: learning across
  traits and sites, arXiv:1504.00946.
}
\keyword{ package }
\keyword{models}
\keyword{regression}
