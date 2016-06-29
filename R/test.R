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
