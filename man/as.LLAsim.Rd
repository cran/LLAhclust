\name{as.LLAsim}
\alias{as.LLAsim}

\title{Converts the lower triangle of a square matrix into a LLAsim
  object}

\description{
Converts the lower triangle of a square matrix into a \code{LLAsim} object. The
\code{LLAsim} object contains similarity coefficients among objects or
variables of interest.}

\usage{
as.LLAsim(m, upper = FALSE, probabilistic = FALSE)
}

\arguments{
  \item{m}{input square similarity matrix.}
  \item{upper}{logical value indicating whether the upper triangle of the
    similarity matrix should be printed by \code{print.LLAsim}.}
  \item{probabilistic}{logical value indicating whether the coefficients
    in the input similarity matrix should be treated as probabilistic
    similarity values. If set to \code{FALSE}, the input similarity
    coefficients are scaled. See examples below.}
}

\value{Returns an object of class \code{LLAsim} whose attributes are
  very similar to those of objects of class \code{dist}. See
  \code{\link[stats]{dist}} for more details.}

\details{The following functions are also defined for objects of class
  \code{LLAsim}: \code{names.LLAsim}, \code{format.LLAsim},
  \code{as.matrix.LLAsim} and \code{print.LLAsim}.}

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
  \code{\link{as.matrix.LLAsim}}, \cr
  \code{\link[stats]{dist}}.}

\examples{
## Assume that we have at hand a probabilistic similarity matrix
## between 5 objects (lower triangle only):
m <- matrix(runif(25), 5, 5)

## The corresponding LLAsim object is obtained as follows:
s <- as.LLAsim(m, probabilistic=TRUE)

## Display the initial matrix and the LLAsim object:
m
s

## Assume now that we have at hand a non-probabilistic similiarty
## matrix:
m <- matrix(rnorm(25), 5, 5)

## The corresponding LLAsim object is obtained as follows:
s <- as.LLAsim(m, probabilistic=FALSE)

## Display the initial matrix and the LLAsim object:
m
s
## Notice that the coefficients in s are scaled:
mean(s)
sd(s)
}

\keyword{cluster}