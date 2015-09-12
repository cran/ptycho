\name{createData}
\alias{createData}
\alias{createDataBayesModel}
\alias{createPubData}
\title{Simulate Data}
\description{
  Create a data object suitable for \code{\link{ptycho.all}}.
}
\usage{
createData(X, y, omega = NULL, beta = NULL)
createDataBayesModel(mode = c("exchange","pleiotropy","gene"), n, p, q,
                     nreps, tau.min, tau.max, G)
createPubData(mode = c("tinysim","ptychoIn",
                       "exchange","pleiotropy","gene",
                       "actualGeno","actualPheno",
                       "fixedOmega","uniformEffects"),
              X=NULL, y=NULL, var.detail=NULL)
}
\arguments{
  \item{X}{Design matrix or alist specifying how to generate such a matrix.  If
    a list, the first entry is a function name and the second is a list of
    arguments to the function.  In \code{createPubData}, \code{X} is ignored
    unless \code{mode} is \dQuote{actualGeno} or \dQuote{actualPheno}.}
  \item{y}{
    Numeric vector or matrix or list with the following components:
    \describe{
      \item{\code{nreps}}{Number of replicates to simulate}
      \item{\code{q}}{Number of responses to generate for each replicate}
      \item{\code{sd}}{Standard deviation of the simulated noise}
    }
    In \code{createPubData}, \code{y} is ignored unless \code{mode} is
    \dQuote{actualPheno}.}
  \item{omega}{Numeric vector or matrix or list specifying how to generate a
    list with the component \code{omega}; see Details for its meaning.  If this
    is a list, the first entry is a function name and the second is a list of
    arguments to the function, which will be prepended by the number of rows in
    output \code{X} and the number of columns in output \code{y}.  Only used if
    \code{y} is a list.}
  \item{beta}{List specifying how to generate a matrix of effect sizes.  The
    first entry of the list is a function name and the second is a list of
    arguments to the function, which will be prepended by a matrix specifying
    the variables selected and \code{y$sd}.  Only used if \code{y} is a list.}
  \item{n}{Number of observations to simulate}
  \item{p}{Number of covariates to simulate}
  \item{q}{Number of responses to simulate for each replicate}
  \item{nreps}{Number of replicates to simulate}
  \item{mode}{String specifying type of dataset to create:
    \describe{
      \item{\code{tinysim}}{Simulated data included with this package;
        equivalent to mode \code{pleiotropy} except that the dataset is tiny,
        with \code{n=100}, \code{p=10}, \code{q=5}, and \code{nreps=10}}
      \item{\code{ptychoIn}}{Simulated data included with this package;
        equivalent to mode \code{gene} except that the dataset is tiny,
        with \code{n=3000}, \code{p=10}, \code{q=1}, and \code{nreps=1}}
      \item{\code{exchange}}{Create orthogonal \code{X} and exchangeable
        variants; \code{n=5000}, \code{p=50}, \code{q=5}, and \code{nreps=100}}
      \item{\code{pleiotropy}}{Create orthogonal \code{X}, and several variants
        have nonzero effects on multiple responses; \code{n=5000}, \code{p=50},
        \code{q=5}, and \code{nreps=100}}
      \item{\code{gene}}{Create orthogonal \code{X}, and each group of variants
        typically has either several or no variants that effect a response;
        \code{n=5000}, \code{p=50}, \code{q=5}, and \code{nreps=100}}
      \item{\code{actualGeno}}{Simulate responses for input \code{X}}
      \item{\code{actualPheno}}{Put input \code{X} and \code{y} into data
        object}
      \item{\code{fixedOmega}}{Create orthogonal \code{X}, and each variant has
        a certain probability of a nonzero effect size}
      \item{\code{uniformEffects}}{Same as mode \code{fixedOmega} except that
        effect sizes are uniformly rather than normally distributed}
    }
    For \code{createDataBayesModel}, \code{mode} must be one of
    \dQuote{exchange}, \dQuote{pleiotropy}, or \dQuote{gene}.}
  \item{tau.min, tau.max}{Endpoints of uniform distribution from which to draw
  \code{tau}}
  \item{G}{Number of groups of covariates; unused if \code{mode} is not
    \dQuote{gene}}
  \item{var.detail}{Data frame with row names same as column names of \code{X};
    must have columns \dQuote{MAF} and \dQuote{GENE}.  Ignored unless
    \code{mode} is \dQuote{actualGeno}.}
}
\details{
  We describe \code{createData} and then describe its wrappers
  \code{createDataBayesModel} and \code{createPubData}.

  Although \code{createData} can form the data object required by
  \code{ptycho.all} when \code{X} and \code{y} are input, it primarily exists to
  simplify simulating data from \eqn{Y=X\beta+\epsilon}{Y=X*beta+epsilon}, where
  \eqn{\epsilon} is normal with mean zero and specified standard deviation and
  \eqn{\beta} is sparse with entries simulated as specified.

  The function generates a specified number of replicates, all of which use the
  same design matrix \eqn{X}.  If this matrix is not input, then its argument
  must specify a function call to generate it.  In either case, suppose \eqn{X}
  has \eqn{n} rows and \eqn{p} columns.

  If the input \code{y} is numeric, then it will be used for the lone replicate.
  If it is a matrix, it must have \eqn{n} rows; let \eqn{q} be its number of
  columns.  If input \code{y} is a numeric vector, it must have \eqn{n} entries
  and will be cast as a matrix with \eqn{q=1} column.  Otherwise, input \code{y}
  is a list specifying, along with the arguments \code{omega} and \code{beta},
  how to simulate the response(s).  Because it is useful in analysis of the
  estimation of the marginal posterior distribution, the returned object always
  contains, regardless of how \code{X} and \code{y} are specified, a matrix
  \code{eta2} with \eqn{(j,k)} entry equal to
  \eqn{\mathbf{x}_j^T \mathbf{y}_k / (n \mathbf{y}_k^T \mathbf{y}_k)}{%
           <x_j,y_k> / (n <y_k,y_k>).}

  If \code{y} is to be simulated, the first step is to choose the probability
  that each covariate is associated with each reponse as specified by the input
  argument \code{omega}.  If this argument is a matrix, it must have size
  \eqn{p}-by-\eqn{q}.  If it is not a matrix but is numeric, it will be passed
  to \code{\link{matrix}} to create a matrix of the correct size.  Otherwise,
  the matrix for each replicate will be generated by calling the function whose
  name is given by \code{omega[[1]]} with argument list
  \code{(p, q, omega[[2]])}.  This function must return a list with component
  \code{omega} set to a \eqn{p}-by-\eqn{q} matrix; the list may also contain
  additional components.  The package contains several functions whose names
  start with \dQuote{createOmega} that might guide users in writing their own
  functions.

  The next step is to draw a \eqn{p}-by-\eqn{q} matrix \code{indic.var} whose
  \eqn{(j,k)} entry is equal to one with probability \code{omega[j,k]} and zero
  otherwise.  This matrix will be drawn until all column sums are positive.

  For each entry in \code{indic.var} that is equal to one, the effect size must
  be drawn.  This is done by calling the function whose name is given by
  \code{beta[[1]]} with argument list \code{(indic.var, y$sd, beta[[2]])}.  This
  function must return a list with component \code{beta} set to a
  \eqn{p}-by-\eqn{q} matrix; the list may also contain additional components.
  If \code{indic.var[j,k]} is zero, then \code{beta[j,k]} should be zero.  The
  package contains functions whose names start with \dQuote{createBeta} that
  might guide users in writing their own functions.

  Finally, an \eqn{n}-by-\eqn{q} matrix of noise is drawn from
  \eqn{N(0,\sigma^2)}, where \eqn{\sigma} is the input \code{noise.sd}, and
  added to \eqn{X\beta}{X*beta} to obtain \code{y}.  The column names of each
  response matrix generated will be \code{y1}, \code{y2}, and so forth.

  The function \code{createPubData} generates the data sets used in Stell and
  Sabatti (2015).  For \code{mode} equal to \dQuote{exchange},
  \dQuote{pleiotropy}, or \dQuote{geno}, it calls \code{createData} via
  \code{createDataBayesModel}; otherwise, it calls \code{createData} directly.
  These functions also serve as additional examples of the use of
  \code{createData}.  For reproducibility, \code{createPubData} first sets the
  random seed to 1234, except that it is set to 4 when \code{mode} equals
  \dQuote{ptychoIn}.

  In \code{createDataBayesModel}, if \code{mode} is \dQuote{exchange}, then
  one \eqn{\omega \sim \mbox{Beta}(12,48)}{omega ~ Beta(12,48)} is drawn
  independently for each trait.  If \code{mode} is \dQuote{pleiotropy}, then one
  probability of association for a trait is drawn from Beta(16,55) for each data
  set, that probability is used to draw \code{indic.grp} for each variant, and
  then the probability of nonzero \code{indic.var[j,k]} is drawn from
  Beta(48,12) for each nonzero \code{indic.grp[j]}.  Finally, if \code{mode} is
  \dQuote{gene}, the process is analogous to pleiotropy except that each trait
  is simulated independently.
}
\value{
  List containing:
  \item{X}{Design matrix}
  \item{q}{Number of columns in each response}
  \item{noise.sd}{Standard deviation of the simulated noise; \code{NULL} if
    input \code{y} is numeric}
  \item{omega}{Input \code{omega}}
  \item{beta}{Input \code{beta}}
  \item{replicates}{List of length \code{y$nreps} (length 1 if \code{y} is
    numeric), each entry of which is a list with the following components:
    \describe{
      \item{\code{omega}}{Matrix containing probabilities of association between
        covariates and responses; row names are \code{colnames(X)} and column
        names are \code{colnames(y)}; \code{NULL} if input \code{y} is numeric}
      \item{\code{indic.var}}{Matrix containing ones for associations and zeros
        otherwise; row and column names are same as for \code{omega};
        \code{NULL} if input \code{y} is numeric}
      \item{\code{beta}}{Matrix of effect sizes; row and column names are same
        as for \code{omega}; \code{NULL} if input \code{y} is numeric}
      \item{\code{y}}{Response matrix}
      \item{\code{eta2}}{Matrix with row names equal to \code{colnames(X)} and
        column names equal to \code{colnames(y)}}
    }
    For \code{createDataBayesModel} with \code{mode} that uses a second level of
    indicator variables, each entry in the \code{replicate} list also has
    components \code{omega.grp} and \code{indic.grp} containing the intermediate
    steps of drawing the second-level indicator variable before drawing
    \code{omega}.  If the argument \code{beta} to \code{createData} is
    \dQuote{createBetaNormal} (which it is when called by
    \code{createDataBayesModel}), then each replicate will also have a component
    \code{tau} giving the value drawn by a call to
    \code{runif(1, tau.min, tau.max)}.
  }
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\references{
  Stell, L. and Sabatti, C. (2015) Genetic variant selection: learning across
  traits and sites, arXiv:1504.00946.
}
\seealso{
  \code{\link{createOrthogonalX}}, \code{\link{createGroupsSim}};
  also \link{Data} describes \code{tinysim} in example below as well as another
  object output by \code{createData}
}
\examples{
### EXAMPLE 1
data(tinysim)
# Data generated with mode equal to pleiotropy, so indic.grp exists and
# has an entry for each column in X.
colnames(tinysim$X)
tinysim$replicates[[5]]$indic.grp
# X4, X6, and X9 are associated with some responses.
tinysim$replicates[[5]]$indic.var

### EXAMPLE 2
# Generate miniature data set with information shared across covariates.
set.seed(1234)
tiny1 <- createDataBayesModel(mode="gene", n=100, p=10, q=5, nreps=10,
                              tau.min=0.045, tau.max=0.063, G=2)
# A covariate can only have indic.var=1 if the group it belongs to has
# indic.grp=1.  For example,indic.grp[1,4]=0 implies
# indic.var[groups$group2var[1],4]=0.
tiny1$replicates[[1]]$indic.grp
tiny1$omega[[2]]$groups$group2var[1]
tiny1$replicates[[1]]$indic.var

### EXAMPLE 3
# Alternatively, call createData directly
groups <- createGroupsSim(G=2, p=10)
omegaargs <- list(indic.grp.shape1=16, indic.grp.shape2=55,
                  shape1=48, shape2=12, groups=groups)
betaargs <- list(tau.min=0.045, tau.max=0.063)
set.seed(1234)
tiny2 <- createData(X=list("createOrthogonalX", list(n=100, p=10)),
                    y=list(nreps=10, q=5, sd=1),
                    omega=list("createOmegaCrossVars", omegaargs),
                    beta=list("createBetaNormal", betaargs))
identical(tiny1, tiny2)
### SEE THE CODE FOR createPubData FOR MORE EXAMPLES.
}
