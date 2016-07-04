####################################################
# These functions are only for comparison purposes #
####################################################

# Incomplete Beta function
incompleteBeta_test <- function(x, alpha, beta) {
    stats::pbeta(x, shape1 = alpha, shape2 = beta) * exp(lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta))
}

# CDF computation
pRoy <- function(x, s, m, n) {
    # Initialize matrix
    A <- matrix(0, s, s)
    b <- rep_len(0, s)

    if (s != 1) {
        for (i in 1:(s-1)) {
            b[i] <- 0.5*incompleteBeta(x, m + i, n + 1)^2
            for (j in i:(s-1)) {
                b[j+1] <- ((m+j)*b[j] - incompleteBeta(x, 2*m+i+j, 2*n+2))/(m+j+n+1)
                A[i,j+1] <- incompleteBeta(x, m+i, n+1) * incompleteBeta(x, m+j+1, n+1) -
                    2*b[j+1]
            }
        }
    }

    if (s %% 2) {
        A <- rbind(A, 0)
        extraCol <- rep_len(0, s+1)
        for (i in 1:s) {
            extraCol[i] <- incompleteBeta(x, m+i, n+1)
        }
        A <- cbind(A, extraCol)
    }
    A <- A - t(A)

    # Compute scaling constant
    c1 <- 0.5*(1:s+2*m+2*n+s+2)
    c2 <- 0.5*(1:s)
    c3 <- 0.5*(1:s+2*m+1)
    c4 <- 0.5*(1:s+2*n+1)

    C <- pi^(0.5*s) * exp(sum(lgamma(c1) - lgamma(c2) -
                                  lgamma(c3) - lgamma(c4)))

    result <- C * sqrt(abs(det(A)))
    return(result)
}

# Multivariate Gamma function
Gm_2 <- function(x, m, log = FALSE) {
    if(log) {
        val <- 0.25*m*(m-1)*log(pi) + sum(lgamma(x - 0.5*0:(m-1)))
    } else {
        val <- pi^(0.25*m*(m-1))*prod(gamma(x - 0.5*0:(m-1)))
    }
    return(val)
}

Gm_ <- function(m,a, log=FALSE){
    if (!log){
        return(pi^(m*(m-1)/4)*Reduce('*',  zipfR::Cgamma(a - ((1:m)-1)/2 )))} else {
            return(log(pi^(m*(m-1)/4)) + Reduce('+',  zipfR::Cgamma(a -
                                                                 ((1:m)-1)/2, log=TRUE)) )
        }
}

#function to calculate CDF of the largest root of Wishart matrix
lrgWishart3<- function(nmin,nmax, x){
    if ((nmin/2)==round(nmin/2)){
        A=matrix(0, nmin, nmin)
        nmat<-nmin
    } else {
        A=matrix(0, nmin+1, nmin+1)
        nmat<-nmin+1
    }
    alph<-(nmax-nmin-1)/2
    b <- rep(0, nmin)

    Klog = Reduce('+', c(((nmin^2)/2)*log(pi) - (nmin*nmax/2)*log(2),
                         -Gm_(nmin, nmax/2, log=T) , -Gm_(nmin, nmin/2, log=T)))
    Kplog = Klog + (alph*nmat+nmat*(nmat+1)/2)*log(2) + Reduce('+',
                                                               zipfR::Cgamma(alph+(1:nmat), log=T))
    Llog<-Kplog*(2/nmat)
    L = exp(Kplog*(2/nmat))

    p<-zipfR::Rgamma(alph+(1:nmin), x/2);
    glog<-zipfR::Cgamma(alph+(1:nmin), log=T)
    qlog<-(log(2)*(-(2*alph+(1:(2*nmin-1))))+zipfR::Cgamma(2*alph+(1:(2*nmin-1)),
                                                    log=T)+zipfR::Rgamma(2*alph+(1:(2*nmin-1)), x, log=T))

    for (i in (1:(nmin-1))){
        b[i] <- (p[i]^2)/2
        for (j in i:(nmin-1)){
            b[j+1] <- b[j] - exp(qlog[i+j]-(glog[i]+glog[j+1]))
            A[i, j+1] = p[i]*p[j+1] - 2*b[j+1]
        }
    }

    if ((nmin/2)!=round(nmin/2)){
        alphi <- alph+(1:nmin)
        logLastcol<-log(2)*(-alph-nmin-1)+zipfR::Rgamma(alphi, x/2,
                                                 log=T)-zipfR::Cgamma(alph+nmin+1, log=T)
        minlogLastcol<-min(logLastcol)
        B<-A
        B[1:nmin,nmin+1]<-exp(logLastcol-minlogLastcol)

        B <- B - t(B);
        logdet <- nmat*Llog+2*minlogLastcol+log(det(B))
        return(list('F'= exp(logdet/2), 'A'=B))
    }

    A <- A - t(A)
    return(list('F'= exp(Kplog + log(det(A))/2), 'A'=A))
}

