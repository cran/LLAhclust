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

## LLA similarities among variables
## Ivan Kojadinovic, April 2007

##############################################################################

similarity.numerical <- function(x)
{
    cat("LLAsim: the variables are supposed to be numerical\n")
    s <- cor(x)
    s <- s[row(s) > col(s)]
    attr(s, "Probabilistic") <- FALSE
    s <- (s - mean(s)) / sd(s)
    return(s)
}

##############################################################################

similarity.categorical <- function(x)
{
    cat("LLAsim: the variables are supposed to be categorical\n")

    p <- ncol(x)
    n <- nrow(x)
    s <- numeric(p * (p - 1) / 2)

    rac1 <- sqrt(2 * n * (n - 1))
    rac2 <- sqrt(n * (n - 1) * (n - 2))
    rac3 <- sqrt(4 * n * (n - 1) * (n - 2) * (n - 3))

    lambda <- numeric(p)
    rho <- numeric(p)
    theta <- numeric(p)
    
    for (i in 1:p)
    {
        m <- summary(as.factor(x[,i]))
        lambda[i] <- sum(m * (m - 1)) / rac1
        rho[i] <- sum(m * (m - 1) * (m - 2)) / rac2
        theta[i] <- (sum(m * (m - 1))^2 -
                          2 * sum(m * (m - 1) * (2 * m - 3))) / rac3
    }
    
    l <- 1
    for (i in 1:p)
        for (j in 1:p)
            if (j > i)
            {
                ct <- table(x[,i],x[,j])
                s[l] <- (sum(ct * (ct - 1) / 2) - lambda[i] * lambda[j]) /
                    sqrt(lambda[i] * lambda[j] + rho[i] * rho[j])# + theta[i] * theta[j] - lambda[i]^2 * lambda[j]^2)
                l <- l + 1
            }

    warning("calculs faits sans 'theta' etc")
    
    attr(s, "Probabilistic") <- FALSE
    s <- (s - mean(s)) / sd(s)
    return(s)
}

##############################################################################

similarity.ordinal <- function(x)
{
    cat("LLAsim: the variables are supposed to be ordinal\n")
    p <- ncol(x)
    n <- nrow(x)
    s <- numeric(p * (p - 1) / 2)

    rac1 <- sqrt(n * (n - 1))
    rac2 <- sqrt(n * (n - 1) * (n - 2))
    rac3 <- sqrt(n * (n - 1) * (n - 2) * (n - 3))
    
    lambda <- numeric(p)
    rhocc <- numeric(p)
    rhoff <- numeric(p)
    rhocf <- numeric(p)
    theta <- numeric(p)

    for (i in 1:p)
    {
        m <- summary(as.factor(x[,i]))
        d <- length(m)
        mm <- m %*% t(m)
        lambda[i] <- sum(mm[lower.tri(mm)]) / rac1

        mc <- cumsum(m)[1:(d-1)]
        mf <- cumsum(m[d:1])[(d-1):1]
        rhocc[i] <- sum(m[2:d] * mc * (mc - 1)) / rac2
        rhoff[i] <- sum(m[1:(d-1)] * mf * (mf - 1)) / rac2
        if (d > 2)
            rhocf[i] <- sum(m[2:(d-1)] * mc[1:(d-2)] * mf[2:(d-1)]) / rac2

        mpm <- mm * (matrix(m,d,d) + matrix(m,d,d,byrow=TRUE)
                     + sum(mm[lower.tri(mm)]) - 2*n + 1)
        theta[i] <- sum( mpm[lower.tri(mpm)]) / rac3
    }
    
    l <- 1
    for (i in 1:p)
        for (j in 1:p)
            if (j > i)
            {
                ct <- table(x[,i],x[,j])

                h <- dim(ct)[1] 
                k <- dim(ct)[2] 
                s[l] <- 0
                for (a in 1:(h-1))
                    for (b in 1:(k-1))
                        s[l] <- s[l] + sum(ct[a,b] * ct[(a+1):h,(b+1):k])

                s[l] <- (s[l] - lambda[i] * lambda[j]) /
                    sqrt(lambda[i] * lambda[j] + rhocc[i] * rhocc[j]
                        + rhoff[i] * rhoff[j] + 2 * rhocf[i] * rhocf[j])
                                        # + theta[i]* theta[j] - lambda[i]^2 * lambda[j]^2 )
                l <- l + 1
            }

    warning("calculs faits sans 'theta' etc")
    
    attr(s, "Probabilistic") <- FALSE
    s <- (s - mean(s)) / sd(s)
    return(s)
}

##############################################################################

similarity.boolean <- function(x)
{
    if (sum(x %in% c(0,1)) != length(x))
        stop("for the considered similarity method, the variables should be boolean, coded by the values 0,1")
 
    cat("LLAsim: the variables are supposed to be boolean, coded by the values 0,1\n")
    n <- margin.table(x,2)
    lambda <- n %*% t(n) / nrow(x)
    s <- (t(x) %*% x - lambda) / sqrt(lambda)
    s <- s[row(s) > col(s)]
    attr(s, "Probabilistic") <- FALSE
    s <- (s - mean(s)) / sd(s)
    return(s)
}

##############################################################################

similarity.chi.square <- function(x)
{
    cat("LLAsim: the variables are supposed to be categorical\n")
    p <- ncol(x)
    s <- numeric(p * (p - 1) / 2)
    l <- 1
    for (i in 1:p)
        for (j in 1:p)
            if (j > i)
            {
                ct <- table(x[,i],x[,j])
                s[l] <- 1 - chisq.test(ct)$p.value
                l <- l + 1
            }
    attr(s, "Probabilistic") <- TRUE
    return(s)        
}

##############################################################################

similarity.cor.abs <- function(x,method)
{
    cat("LLAsim: the variables are supposed to be numerical or ordinal\n")
    p <- ncol(x)
    s <- numeric(p * (p - 1) / 2)
    l <- 1
    for (i in 1:p)
        for (j in 1:p)
            if (j > i)
            {
                s[l] <- 1 - cor.test(x[,i],x[,j],method=method)$p.value
                l <- l + 1
            }
    attr(s, "Probabilistic") <- TRUE
    return(s)        
}

##############################################################################

# simulate the distribution of the empirical coopula statistic under independence

empcopula.simulate <- function(n,N=2000)
{
    if (!is.numeric(n) || as.integer(n) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(N) || as.integer(N) < 100)
        stop("N should be an integer greater than 100")

    Bn <- .C("simulate_empirical_copula",
             n = as.integer(n),
             Bn = double(N),
             N = as.integer(N),
             PACKAGE="LLAhclust")$Bn

    simulated.dist <- list(sample.size = n,
                           number.repetitions = N,
                           dist.independence = Bn)
    
    class(simulated.dist) <- "empcopula.simulation"
           
    return(simulated.dist)
}

##############################################################################

similarity.empirical.copula <- function(x, simulated.distribution=NULL)
{
    #cat("LLAsim: the variables are supposed to be numerical\n")

    p <- ncol(x)
    n <- nrow(x)
    # transform data to ranks
    for (j in 1:p)
        x[,j] <- rank(x[,j])

    if (is.null(simulated.distribution))
    {
        N <- 2000
        Bn <- .C("simulate_empirical_copula",
                 n = n,
                 Bn = double(N),
                 N = as.integer(N),
                 PACKAGE="LLAhclust")$Bn
    }
    else
    {
        Bn <- simulated.distribution$dist.independence
        N <- simulated.distribution$number.repetitions
    }
    
    s <- .C("similarity_empirical_copula",
            x = as.integer(x),
            n = n,
            p = p,
            s = double(p*(p - 1)/2),
            Bn = Bn,
            N = as.integer(N),
            PACKAGE="LLAhclust")$s

    attr(s, "Probabilistic") <- TRUE
    return(s)
}

##############################################################################

LLAsimvar <- function(x, method="LLAnumerical", upper=FALSE, simulated.distribution=NULL)
{
    METHODS <- c("LLAnumerical","LLAcategorical","LLAordinal","LLAboolean",
                 "chi.square","pearson.abs","spearman.abs","kendall.abs",
                 "empirical.copula")
    method <- pmatch(method, METHODS)
    if(is.na(method))
	stop("invalid similarity method")
    if(method == -1)
	stop("ambiguous similarity method")
    if (!is.numeric(x <- as.matrix(x)))
        stop("objects should be described by numerical values")
    if (sum(is.na(x)) > 0)
        warning("objects descriptions contain missing values")
    if (sum(sd(x,na.rm=TRUE) < 1e-16) > 0)
        stop("some variables seem to always take the same value")
    if (!is.null(simulated.distribution))
    {
        if (class(simulated.distribution) != "empcopula.simulation")
            stop("simulated.distribution should be obtained by means of the function empcop.simulate")
        if (simulated.distribution$sample.size != nrow(x))
            warning("simulated.distribution corresponds to a different sample size")
    }

    s <- switch(method,
                similarity.numerical(x),
                similarity.categorical(x),
                similarity.ordinal(x),
                similarity.boolean(x),
                similarity.chi.square(x),
                similarity.cor.abs(x,method="pearson"),
                similarity.cor.abs(x,method="spearman"),
                similarity.cor.abs(x,method="kendall"),
                similarity.empirical.copula(x,simulated.distribution))

    attr(s, "Size") <- ncol(x)
    attr(s, "Labels") <- dimnames(x)[[2]]
    attr(s, "Upper") <- upper
    attr(s, "method") <- paste("LLAsimvar:",METHODS[method])
    attr(s, "call") <- match.call()
    class(s) <- "LLAsim"
           
    return(s)
}

##############################################################################
