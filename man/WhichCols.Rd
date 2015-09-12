\name{WhichCols}
\alias{WhichCols}
\alias{indicVarCols}
\alias{indicGrpCols}
\title{Identify Columns Containing Indicator Variables}
\description{
  Determine which columns contain \code{indic.var} or \code{indic.grp} in an
  object returned by \code{\link{ptycho}}.
}
\usage{
indicVarCols(obj)
indicGrpCols(obj)
}
\arguments{
  \item{obj}{Object output by \code{\link{ptycho}} or any numeric object}
}
\value{
  If the input object is a numeric vector, returns the indices of its entries
  that have names starting with \dQuote{indic.var} or \dQuote{indic.grp},
  respectively.

  Otherwise, it returns the indices of \code{colnames(obj)} that start with
  \dQuote{indic.var} or \dQuote{indic.grp}, respectively.
}
\author{
  Laurel Stell and Chiara Sabatti\cr
  Maintainer: Laurel Stell <lstell@stanford.edu>
}
\seealso{
  \code{\link{ptycho}}; also \code{\link{ptychoOut}} and
  \link{PosteriorStatistics} for example below
}
\examples{
data(ptychoOut)
colnames(ptychoOut)[indicVarCols(ptychoOut)]
# Can also apply these functions to output of meanIndicators ...
mi <- meanIndicators(ptychoOut)
mi[indicGrpCols(mi)]
# ... instead of using meanGrpIndicators or meanVarIndicators
meanGrpIndicators(ptychoOut)
}
