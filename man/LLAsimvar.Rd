\name{LLAsimvar}
\alias{LLAsimvar}

\title{Computes similarities among variables using the likelihood
  linkage analysis approach}

\description{
  Computes similarities among variables using the likelihood
  linkage analysis approach proposed by Lerman. The likelihood
  linkage analysis method mainly consists in replacing the value
  of the estimated similarity coefficient between two variables
  by the probability of finding a lower value under the
  hypothesis of stochastic independence, called
  \emph{absence of link} in that context. Nine similarity coefficients
  can be computed using the \code{LLAsimvar} function.}

\usage{
LLAsimvar(x, method = "LLAnumerical", upper = FALSE,
                 simulated.distribution = NULL)
}

\arguments{
  \item{x}{a numeric matrix or data frame.}
  \item{method}{Can be one of \code{LLAnumerical},
    \code{LLAcategorical}, \code{LLAordinal}, \code{LLAboolean},
    \code{chi.square}, \code{pearson.abs},
    \code{spearman.abs}, \code{kendall.abs} or
    \code{empirical.copula}. The methods \code{LLA*} were initially
    defined by Lerman (see references below). The four remaining methods
  compute the similarity between two variables as one minus the p-value
  obtained from a test of independence. See the last reference and the
  example section below for more details.}
  \item{upper}{logical value indicating whether the upper triangle of the
    similarity matrix should be printed by \code{print.LLAsim}.}
  \item{simulated.distribution}{Object of class
    \code{empcopula.simulation}.
    Should be set only if the method \code{empirical.copula} is selected.
    See function \code{\link{empcopula.simulate}} and the example
    section below for more details. 
  }
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

  P. Deheuvels (1979), La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, Acad. Roy.
  Belg. Bull. Cl. Sci. 5th Ser. 65, 274-292.
   
  C. Genest and B. Rémillard (2004). \emph{Tests of independence and
    randomness based on the empirical copula process.} Test, 13, 335-369.

  I. Kojadinovic (2007), \emph{Hierarchical clustering of continuous variables
  based on the empirical copula process}, submitted.

}

\seealso{
  \code{\link{as.LLAsim}}, \cr
  \code{\link{empcopula.simulate}}, \cr
  \code{\link{LLAsimobj}}, \cr
  \code{\link{LLAhclust}}, \cr
  \code{\link{LLAparteval}}, \cr
  \code{\link[stats]{dist}}.} 


\examples{
data(USArrests)

## Compute similarities between variables using the
## LLAnumerical method:
s <- LLAsimvar(USArrests)
s

## Compute similarities between variables using the classical
## bilateral test of independence based on Spearman's rho:
s <- LLAsimvar(USArrests, method = "spearman.abs")
s

## Compute similarities between variables using the classical
## bilateral test of independence based on Kendall's tau:
s <- LLAsimvar(USArrests, method = "kendall.abs")
s

## Compute similarities between variables using the test of
## independence à la Deheuvels based on the empirical copula
## process recently studied by Genest and Rémillard:
s <- LLAsimvar(USArrests, method = "empirical.copula")
s

## The previous computation could have been done in two steps:
d <- empcopula.simulate(n=50,N=2000)
s <- LLAsimvar(USArrests, method = "empirical.copula",
                       simulated.distribution = d)
s
}


\keyword{cluster}
