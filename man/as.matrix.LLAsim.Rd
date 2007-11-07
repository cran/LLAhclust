\name{as.matrix.LLAsim}
\alias{as.matrix.LLAsim}
\alias{format.LLAsim}
\alias{names.LLAsim}
\alias{print.LLAsim}

\title{Useful functions for dealing with LLAsim objects}

\description{
 The function \code{as.matrix.LLAsim} converts a \code{LLAsim} object into a square
 symmetrical matrix. The usual R functions \code{format}, \code{print} and
\code{names} have also been extended to deal with \code{LLAsim} objects.}

\usage{
as.matrix.LLAsim(x, ...)
}

\arguments{
  \item{x}{the \code{LLAsim} object to be converted.}
  \item{\dots}{nothing so far.}
}

\value{An object of class \code{matrix}.}

\references{
  I.C. Lerman (1981), \emph{Classification et analyse ordinale de donnés},
  Dunod, Paris.
  
  I.C. Lerman (1991), \emph{Foundations of the likelihood linkage analysis
    classification method}, Applied Stochastic Models and Data
  Analysis, 7, pages 63--76.

  I.C. Lerman (1993), \emph{Likelihood linkage analysis
    classification method: An example treated by hand},
  Biochimie, 75, pages 379--397.
  
  I.C. Lerman, Ph. Peter and H. Leredde (1993), \emph{Principes et calculs de
    la méthode implantée dans le programme CHAVL (Classification
    Hiérarchique par Analyse de la Vraisemblance des Liens)}, Modulad, 12,
  pages 33-101.
}

\seealso{\code{\link{LLAsimvar}}, \cr
  \code{\link{LLAsimobj}}, \cr
  \code{\link{as.LLAsim}}.}

\examples{
data(USArrests)

## Compute similarities between objects based on
## a local Euclidean distance (see references above):
s <- LLAsimobj(USArrests)

## Convert to a matrix object:
as.matrix(s)

## Other useful functions:
print(s, upper=TRUE)
names(s)
## For the format function, see the R help.
}

\keyword{cluster}