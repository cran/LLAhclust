\name{empcopula.simulate}
\alias{empcopula.simulate}

\title{Simulation step used in the independence test  based on the empirical 
  copula process implemented in the LLAsimvar function}

\description{
  Simulation step used in the independence test  based on the empirical 
  copula process as proposed by Christian Genest and Bruno 
  Rémillard. To be used in conjunction with the \code{LLAsimvar}
  function (\code{method="empirical.copula"}). The simulation step
  consists in simulating the
  distribution of the test statistic under independence for the sample
  size under consideration. More details can be found in the articles
  cited in the reference section.
}

\usage{
empcopula.simulate(n, N = 2000)
}

\arguments{
  \item{n}{Sample size when simulating the distribution of the test
    statistic under independence.}
  \item{N}{Number of repetitions when simulating under independence.}
}

\details{See the references below for more details, especially the third
  one.
}

\value{The function \code{empcopula.simulate} returns an object of class
  \code{empcop.simulation} whose attributes are: \code{sample.size},
  \code{number.repetitons} and
  \code{dist.independence} (a vector of length \code{N} containing
  the values of the test statistic for each each repetition).
}

\references{
  P. Deheuvels (1979), La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, Acad. Roy.
  Belg. Bull. Cl. Sci. 5th Ser. 65, 274-292.
  
  P. Deheuvels (1981), \emph{A non parametric test for independence},
  Publ. Inst. Statist. Univ. Paris 26, 29-50.
  
  C. Genest and B. Rémillard (2004). \emph{Tests of independence and
    randomness based on the empirical copula process.} Test, 13, 335-369.
  
  C. Genest, J.-F. Quessy and B. Rémillard (2006).  \emph{Local
    efficiency of a Cramer-von Mises test of independence.}
  Journal of Multivariate Analysis, 97, 274-294.
  
  C. Genest, J.-F. Quessy and B. Rémillard (2007). \emph{Asymptotic
    local efficiency of Cramer-von Mises tests for multivariate
    independence.} The Annals of Statistics, 35, in press.
  
  I. Kojadinovic (2007), \emph{Hierarchical clustering of continuous variables
    based on the empirical copula process}, submitted.
}

\seealso{\code{\link{LLAsimvar}}, \cr \code{\link{LLAhclust}}. }

\examples{
data(USArrests)

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
