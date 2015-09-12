\name{createOrthogonalX}
\alias{createOrthogonalX}
\title{Create Design Matrix With Orthogonal Columns}
\description{
  Create a design matrix whose columns are orthogonal to each other.
}
\usage{
createOrthogonalX(n, p)
}
\arguments{
  \item{n}{Number of rows in \eqn{X}}
  \item{p}{Number of columns in \eqn{X}}
}
\details{
  First create
  \eqn{\hat{X} = (I_p\,\,\,I_p\,\,\,\cdots\,\,\,I_p)^T,}{X = (I_p  I_p  ...  I_p)^T,}
  where \eqn{I_p}, the identity matrix of size \eqn{p}, is repeated
  \code{ceiling(n/p)} times.  If \code{p} does not divide \code{n}, remove
  rows at the bottom so that \eqn{\hat{X}}{X} has \code{n} rows.  Divide by the
  root mean square of the columns of \eqn{\hat{X}}{X}.
}
\value{
  Matrix with \code{n} rows, \code{p} columns, and column names \code{X1},
  \code{X2}, and so forth.
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{\code{\link{createData}}}
\examples{
n <- 50; p <- 5
X <- createOrthogonalX(n, p)
XtX <- t(X) \%*\% X
D <- diag(n-1, nrow=p)
# XtX and D are not quite equal due to roundoff error
range(XtX - D)
}
