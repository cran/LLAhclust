\name{LLAparteval}
\alias{LLAparteval}

\title{Evalutes the quality of each partition compatible with a
  hierarchy in terms of several indices}

\description{
Evalutes the quality of each partition compatible with the hierarchy
returned by \code{LLAhclust} . If the hierarchy is obtained from
similarity coefficients computed using \code{LLA*} methods, the global
and local statistics proposed by Lerman are calculated. Otherwise, for
similarity coefficients obtained from independence tests
(see \code{\link{LLAsimvar}}),  for each partition, the inter-class
p-values are combined using Tippett's and Fisher's rules. Furthermore,
the minimum inter-class p-value and the maximum intra-class p-value are
given. See the last reference and the examples below for more details.
}

\usage{
LLAparteval(tree, s, m=NULL)
}

\arguments{
  \item{tree}{An object of class \code{hclust} as returned by \code{LLAhclust}.}
  \item{s}{An object of class \code{LLAsim} as returned by
    \code{LLAsimvar}, \code{LLAsimobj} or \code{as.LLAsim}.}
  \item{m}{Integer. If set, the quality of the m coarsest partitions
    only is evaluated.}
}


\value{Returns a \code{data.frame} whose columns are: \code{global.stat}
  and \code{local.stat} if the hierarchy is obtained from
  similarity coefficients computed using \code{LLA*} methods, and
  \code{tippett.inter}, \code{fisher.inter}, \code{min.inter} and
  \code{max.intra} in case of similarity coefficients obtained from
  independence tests.
}

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

  I. Kojadinovic (2007), \emph{Hierarchical clustering of continuous variables
  based on the empirical copula process}, submitted.
}

\seealso{\code{\link{LLAsimvar}}, \cr
  \code{\link{LLAsimobj}}, \cr
  \code{\link{as.LLAsim}}, \cr
  \code{\link{LLAhclust}}.}

  \examples{

data(USArrests)

## Compute similarities between variables based on
## the LLAnumerical method:
s <- LLAsimvar(USArrests)
s

## Perform the hierarchical clustering of the variables:
h <- LLAhclust(s)
plot(h)

## Compute the quality of the partitions compatible
## with the hierarchy in terms of the statistics defined by Lerman:
LLAparteval(h,s)

## Compute similarities between variables using the classical
## bilateral test of independence based on Spearman's rho:
s <- LLAsimvar(USArrests, method = "spearman.abs")
s

## Perform the hierarchical clustering of the variables
## using Fisher's p-value combination method:
h <- LLAhclust(s,method="fisher")
plot(h)
## NB: the height in the dendrogram is a p-value
## and can be used to identify mutually independent classes of
## variables, if any.

## Compute the quality of the partitions compatible
## with the hierarchy in terms of the indices defined in the
## last reference:
LLAparteval(h,s)
}

\keyword{cluster}
