##############################################################################
#
# Copyright © 2007 Ivan Kojadinovic, Israël-César Lerman and Philippe Peter 
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

##############################################################################
##
## Hierarchical clustering. Initial version from R package "stats".
## Credit to the R people! 
## A range of criteria are supported; also there is a
## storage-economic option.
##
## This uses the very efficient nearest neighbor chain algorithm,
## which makes this algorithm of O(n^2) computational time,
## and differentiates it from the less efficient -- i.e. O(n^3) -- implementations
##
## Clustering (linkage) methods: LLA linkages
##
## Original author: F. Murtagh, May 1992
## R Modifications: Ross Ihaka, Dec 1996
##	               Friedrich Leisch, Apr 1998, Jun 2000
## LLA linkage methods: Ivan Kojadinovic, Apr 2007
##
##############################################################################

LLAhclust <- function(s, method="lla", epsilon=1, members=NULL)
{
    METHODS <- c("lla", "tippett",  "average", "complete", "fisher",
                 "uniform",  "normal", "maximum")
    method <-  pmatch(method, METHODS)
    if(is.na(method))
	stop("invalid clustering method")
    if(method == -1)
	stop("ambiguous clustering method")
    
    n <- as.integer(attr(s, "Size"))
    probabilistic <- attr(s, "Probabilistic")
    if(is.null(n) || is.null(probabilistic))
	stop("invalid similarities")
    
    if(n < 2)
        stop("must have n >= 2 objects to cluster")

    len <- as.integer(n*(n-1)/2)
    if(length(s) != len)
        (if (length(s) < len)stop else warning)("similarities of improper length")

    if(epsilon > 1 || epsilon < 0)
        stop("epsilon should lie in the interval [0,1]")
    
    if(is.null(members))
        members <- rep(1, n)
    else if(length(members) != n)
        stop("invalid length of members")

    if (probabilistic == FALSE)
        s <- pnorm(s)
    if (any(s == 0) || any(s == 1))
        (if (method == 1)stop else warning)("some similarities are 0 or 1")
    
    if (method == 1)
        s <- -log(-log(s))
        
    if (method == 6)
        s2 <- s
    else
        s2 <- NULL

    if (method <= 4)
        hcl <- .Fortran("llahclustnn",
                        n = n,
                        len = len,
                        iopt = as.integer(method),
                        ia = integer(n),
                        ib = integer(n),
                        crit = double(n),
                        membr = as.double(members),
                        nn = integer(n),
                        simnn = double(n),
                        flag = logical(n),
                        simi = as.double(s),
                        epsilon = as.double(epsilon),
                        PACKAGE="LLAhclust")
    else
        hcl <- .Fortran("llahclustn3",
                        n = n,
                        len = len,
                        iopt = as.integer(method),
                        ia = integer(n),
                        ib = integer(n),
                        crit = double(n),
                        membr = as.double(members),
                        flag = logical(n),
                        simi = as.double(s),
                        simi2 = as.double(s2),
                        PACKAGE="LLAhclust")
   
    ## 2nd step: interpret the information that we now have
    ## as merge, height, and order lists.

    hcass <- .Fortran("hcass2",
                      n = as.integer(n),
                      ia = as.integer(hcl$ia),
                      ib = as.integer(hcl$ib),
                      order = integer(n),
                      iia = integer(n),
                      iib = integer(n),
                      PACKAGE="LLAhclust")

    tree <- list(merge = cbind(hcass$iia[1:(n-1)], hcass$iib[1:(n-1)]),
                 height= hcl$crit[1:(n-1)],
                 order = hcass$order,
                 labels=attr(s, "Labels"),
                 method=METHODS[method],
                 call = match.call(),
                 LLAsim.method = attr(s, "method"))
    class(tree) <- "hclust"

    return(tree)
}

##############################################################################


