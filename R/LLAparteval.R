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
##
## LLApareval: Ivan Kojadinovic, May 2007
##
##############################################################################

## CDF of the sum of uniform independent [0,1] random variables

psumunif <- function(y,n)
{
    .C("cdf_sum_uniform",
       as.double(y),
       as.integer(n),
       F = double(1),
       PACKAGE="LLAhclust")$F
}

##############################################################################

## Computes the quality of the partitions compatible with the hierarchy

LLAparteval <- function(tree, s, m=NULL)
{
    if (class(tree) != "hclust")
        stop("invalid tree")
    
    if (class(s) != "LLAsim")
        stop("invalid similarities")
    
    n <- as.integer(attr(s, "Size"))
    probabilistic <- attr(s, "Probabilistic")
    if(is.null(n) || is.null(probabilistic))
	stop("invalid similarities")
    
    len <- as.integer(n*(n-1)/2)
    if(length(s) != len)
        (if (length(s) < len) stop else warning)("similarities of improper length")
    
    if (length(tree$height) != n -1)
        stop("the tree was not obtained from these similarlities")

    if (!is.null(m))
    {
        if(!is.numeric(m <- as.integer(m)) || m < 2 || m > n)
            stop(paste("m should be an integer greater than 2 and lower than",n))
    }
    else
        m <- n
    
    if (probabilistic == FALSE) ## LLA partition index
    {   
        s <- as.matrix(s)

        ## computations necessary for the variance
        A1 <- sum(s)^2
        A2 <- sum(margin.table(s,1)^2)
        A3 <- sum(s^2)
        d1 <- n * (n - 1)
        d2 <- n * (n - 1) * (n - 2)
        d3 <- n * (n - 1) * (n - 2) * (n - 3)
        
        global <- numeric(m)
        local <- numeric(m)
        local[m] <- NA

        ## for each partition 
        for (i in m:2)
        {
            f <- matrix(0,n,n)
            part <- cutree(tree,i)
            for (j in 1:i)
            {
                p <- part
                p[p!=j] <- 0 
                p[p==j] <- 1
                f <- f + p %*% t(p) 
            }
            diag(f) <- numeric(n)
            
            B1 <- sum(f)^2
            B2 <- sum(margin.table(f,1)^2)
            B3 <- sum(f)
            
            v <- - A1 * B1 / d1^2 + 2 * A3 * B3 / d1 + 4 * (A2 - A3) * (B2 - B3) / d2 +
                (A1 - 4 * A2 + 2 * A3) * (B1 - 4 * B2 + 2 * B3) / d3

            if (v <= 0)
                global[i] <- NA
            else
                global[i] <- (sum(s * f) - sqrt(A1 * B1) / d1)/sqrt(v)
            if (i < m)
                local[i] <- global[i] - global[i+1]
        }

        res <- data.frame(global.stat = global[m:2],local.stat = local[m:2])
        row.names(res) <- m:2
        return(res)
    }
    else ## p-value combination methods
    {
        fisher.inter <- rep(NA,m)
        tippett.inter <- rep(NA,m)
        min.inter <- rep(NA,m)
        max.intra <- rep(NA,m)
        
        ## for each partition 
        for (i in m:1)
        {
            f <- matrix(0,n,n)
            part <- cutree(tree,i)
            for (j in 1:i)
            {
                p <- part
                p[p!=j] <- 0 
                p[p==j] <- 1
                f <- f + p %*% t(p) 
            }
            f <- f[row(f) > col(f)]

            if (i > 1)
            {
                len.inter <- length(f[f == 0])
                tippett.inter[i] <- 1 - max(s[f == 0])^len.inter
                fisher.inter[i] <- pchisq(sum(-2 * log(1 - s[f == 0])), 2 * len.inter, lower.tail = FALSE)
                min.inter[i] <- min(1 - s[f == 0])
            }

            if (i < n)
                max.intra[i] <- max(1 - s[f == 1])
        }
        
        res <- data.frame(tippett.inter = tippett.inter[m:1], fisher.inter = fisher.inter[m:1], 
			min.inter = min.inter[m:1], max.intra = max.intra[m:1])

        row.names(res) <- m:1
        return(res)
    }
}

#############################################################################
