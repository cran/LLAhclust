\name{LLAhclust}
\alias{LLAhclust}

\title{Likelihood linkage analysis hiearichal clustering}

\description{Builds a hierarchy from similarity coefficients among objects
  or variables  as returned by \code{LLAsimvar}, \code{LLAsimobj}
  or \code{as.LLAsim}. The default aggregation criteria, called
  \code{lla}, can be regarded as a probabilistic version of the single linkage.
}
\usage{
LLAhclust(s, method = "lla", epsilon = 1, members = NULL)
}
\arguments{
  \item{s}{Similarity coefficients as returned by \code{LLAsimvar},
     \code{LLAsimobj} or \code{as.LLAsim}.}
  \item{method}{Linkage method (i.e. aggregation criterion). Can be one
    of \code{lla} (default), \code{tippett} (Tippett's p-value
    combination method), \code{average}, \code{complete}, \code{fisher}
    (Fisher's p-value combination method), \code{uniform}
    (uniform p-value combination method; can be regarded as a
    probabilistic version of the average linkage), \code{normal}
    (normal p-value combination method) or \code{maximum} (maximum p-value combination
    method; can be regarded as a probabilistic version of the complete
    linkage). See the last reference for more details.
  }
  \item{epsilon}{Coefficient used in the \code{lla} linkage. Should lie
    in [0,1]: \code{epsilon=0}
    corresponds to the single linkage, \code{epsilon=1} (default) yields
  a probabilistic version of the single linkage.}
  \item{members}{"Weights" of the objects to be clustered if not of
    equal "weight". See \code{\link[stats]{hclust}} for more details.}
}

\value{An object of class \code{hclust} with the corresponding
  attributes. See \code{\link[stats]{hclust}} for more details.}

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
  \code{\link{LLAparteval}}, \cr
  \code{\link[stats]{hclust}}.}


\examples{
data(USArrests)

## Compute similarities between variables based on
## the LLAnumerical method:
s <- LLAsimvar(USArrests)
s

## Perform the hierarchical clustering of the variables
## using the default aggregation criterion (lla):
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
