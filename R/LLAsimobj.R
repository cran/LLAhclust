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

## Modification of the original R code "dist.R". Credit to the R people.
## LLA similarities among objects
## Ivan Kojadinovic, April 2007

LLAsimobj <- function(x, method="LLAeuclidean", upper=FALSE)
{
    METHODS <- c("LLAeuclidean","LLAcosinus","LLAcategorical","LLAordinal",
                 "LLAboolean")
    method <- pmatch(method, METHODS)
    if(is.na(method))
	stop("invalid similarity method")
    if(method == -1)
	stop("ambiguous similarity method")
    n <- nrow(x <- as.matrix(x))
    if (!is.numeric(x))
        stop("objects should be described by numerical values")
    if (sum(is.na(x)) > 0)
        stop("objects descriptions contain missing values")
    if (sum(sd(x,na.rm=TRUE) < 1e-16) > 0)
        stop("some variables seem to always take the same value")
    if (method %in% c(3,4) && !is.integer(x))
        stop("for the considered similarity method, objects should be described by integer values")
    if (method == 5 && (sum(x %in% c(0,1)) != length(x)))
        stop("for the considered similarity method, objects should be described by 0,1 values")
    
    s <- .C("similarity_objects",
            x = as.double(x),
            nr = n,
            nc = ncol(x),
            s = double(n*(n - 1)/2),
            method = as.integer(method),
            PACKAGE="LLAhclust")$s
    
    attr(s, "Size") <- n
    attr(s, "Probabilistic") <- FALSE
    attr(s, "Labels") <- dimnames(x)[[1]]
    attr(s, "Upper") <- upper
    attr(s, "method") <- paste("LLAsimobj:",METHODS[method])
    attr(s, "call") <- match.call()
    class(s) <- "LLAsim"
    return(s)
}

names.LLAsim <- function(x) attr(x, "Labels")

"names<-.LLAsim" <- function(x, value)
{
    if(length(value) != attr(x, "Size"))
	stop("invalid names for LLAsim object")
    attr(x, "Labels") <- value
    x
}

## Because names(o) != length(o) for "LLAsim"-object o, we need
format.LLAsim <- function(x, ...) format(as.vector(x), ...)

as.matrix.LLAsim <- function(x, ...)
{
    size <- attr(x, "Size")
    df <- matrix(0, size, size)
    df[row(df) > col(df)] <- x
    df <- df + t(df)
    labels <- attr(x, "Labels")
    dimnames(df) <-
	if(is.null(labels)) list(1:size,1:size) else list(labels,labels)
    df
}

as.LLAsim <- function(m, upper = FALSE, probabilistic = FALSE)
{
    if (inherits(m,"LLAsim"))
	ans <- m
    else { ## matrix |-> LLAsim
	m <- as.matrix(m)
        if(!is.numeric(m)) # coerce w/o losing attributes
            storage.mode(m) <- "numeric"
        p <- nrow(m)
        if(ncol(m) != p) warning("non-square matrix")
	ans <- m[row(m) > col(m)]
	attributes(ans) <- NULL
	if(!is.null(rownames(m)))
	    attr(ans,"Labels") <- rownames(m)
	else if(!is.null(colnames(m)))
	    attr(ans,"Labels") <- colnames(m)
	attr(ans,"Size") <- p
	attr(ans, "call") <- match.call()
	class(ans) <- "LLAsim"
    }
    if(is.null(attr(ans,"Upper")) || !missing(upper))
	attr(ans,"Upper") <- upper
    if(is.null(attr(ans,"Probabilistic")) || !missing(probabilistic))
	attr(ans,"Probabilistic") <- probabilistic
    if (probabilistic == TRUE)
    {
        if (any(ans > 1) || any(ans < 0))
            stop("some coefficients are lower than 0 or larger than 1")
    }
    else
        ans <- (ans - mean(ans))/sd(ans)  
    ans
}


print.LLAsim <- function(x, upper = NULL,
                         digits = getOption("digits"), justify = "none", right = TRUE, ...)
{
    if(is.null(upper))
	upper <- if(is.null(a <- attr(x,"Upper"))) FALSE else a

    m <- as.matrix(x)
    cf <- format(m, digits = digits, justify = justify)
    if(!upper)
	cf[row(cf) < col(cf)] <- ""
    ## No diag
    cf[row(cf) == col(cf)] <- ""

### Better: use an improved prettyNum() function -> ../../base/R/format.R

##-     if(any((i <- m == floor(m))))
##-         cf[i] <- sub("0+$", "", cf[i])
    print(if(upper) cf else cf[-1, -attr(x, "Size")],
	  quote = FALSE, right = right, ...)
    invisible(x)
}
