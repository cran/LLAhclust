\name{LLAsimobj}
\alias{LLAsimobj}

\title{Computes similarities among objects}

\description{Computes similarities among objects using the likelihood
  linkage analysis approach proposed by Lerman. The likelihood
  linkage analysis method mainly consists in replacing the value
  of the similarity coefficient between two objects
  by the probability of finding a lower value under the
  hypothesis of \emph{absence of link}. See the references below for
  more details.}

\usage{
LLAsimobj(x, method = "LLAeuclidean", upper = FALSE)
}

\arguments{
  \item{x}{a numeric matrix or data frame.}
  \item{method}{Can be one of \code{LLAeuclidean},
    \code{LLAcosinus}, \code{LLAcategorical}, \code{LLAordinal},
    or \code{LLAboolean}. The two first methods can be used to compute
    similiarty coefficients between objects described by numerical variables.}
  \item{upper}{logical value indicating whether the upper triangle of the
    similarity matrix should be printed by \code{print.LLAsim}.}
}
      
\details{The following functions are also defined for objects of class
  \code{LLAsim}: \code{names.LLAsim}, \code{format.LLAsim},
  \code{as.matrix.LLAsim} and \code{print.LLAsim}.}
      
\value{Returns an object of class \code{LLAsim} whose attributes are
  very similar to those of objects of class \code{dist}. See
  \code{\link[stats]{dist}} for more details.}

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
  \code{\link{as.LLAsim}}, \cr
  \code{\link{LLAhclust}}, \cr
  \code{\link{LLAparteval}}, \cr
  \code{\link[stats]{dist}}.}


\examples{
data(USArrests)

## Compute similarities between objects based on
## a local Euclidean distance (see references above):
s <- LLAsimobj(USArrests)
s
}

\keyword{cluster}