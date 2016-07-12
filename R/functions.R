#' Distribution of the largest root
#'
#' Computes the cumulative distribution function of the largest root in the
#' single and double Wishart setting.
#'
#' If \eqn{S} follows a Wishart(\eqn{p,n}) distribution, e.g. if we can write
#' \deqn{S = X^TX,} where \eqn{X} is an \eqn{n x p} matrix with i.i.d rows
#' coming from a \eqn{p}-variate standard normal, then \code{singleWishart}
#' gives the distribution of the largest root of \eqn{S}.
#'
#' As its name indicates, the double Wishart setting involves two Wishart
#' variables: let \eqn{A} and \eqn{B} be Wishart(\eqn{p,m}) and
#' Wishart(\eqn{p,n}), respectively. If \eqn{A+B} is invertible, then
#' \code{doubleWishart} gives the distribution of the largest root of
#' \deqn{(A+B)^-1B.} Alternatively, it gives the distribution of the largest
#' root of the determinental equation \deqn{det(B - \theta(A+B)).}
#'
#' @param x Vector of numeric values at which to compute the CDF.
#' @param p,n,m Parameters of the single and double Wishart settings. See
#'   details.
#' @param mprec Logical. Should we perform high precision computations?
#' @return Returns the value of the CDF at \code{x}.
#' @examples
#' x1 <- seq(0, 30, length.out = 50)
#' y1 <- singleWishart(x1, 5, 10, mprec = FALSE)
#' plot(x1, y1, type='l')
#'
#' x2 <- seq(0, 1, length.out = 50)
#' y2 <- doubleWishart(x2, 10, 10, 200, mprec = FALSE)
#' plot(x2, y2, type='l')
#' @export
#' @aliases doubleWishart singleWishart
#' @rdname largestRoot
singleWishart <- function(x, p, n, mprec = TRUE) {
    # Check input
    stopifnot(all(x >= 0), p > 0, n > 0, isWhole(p), isWhole(n))

    n_min <- min(p, n)
    n_max <- max(p, n)

    if (mprec) {
        result <- singleWishart_raw(x, n_min, n_max, mprec)
    } else {
        # There are some special cases
        specialCases <- list(c(2,2),
                             c(2,5),
                             c(3,3),
                             c(4,4))
        if(any(vapply(specialCases,
                      function(pair) all(pair == c(n_min, n_max)),
                      logical(1)))) {
            if (all(c(n_min, n_max) == c(2,2))) result <- F22(x)
            if (all(c(n_min, n_max) == c(2,5))) result <- F25(x)
            if (all(c(n_min, n_max) == c(3,3))) result <- F33(x)
            if (all(c(n_min, n_max) == c(4,4))) result <- F44(x)
        } else {
            result <- singleWishart_raw(x, n_min, n_max, FALSE)
        }
    }

    return(result)
}

#' @export
#' @rdname largestRoot
doubleWishart <- function(x, p, n, m, mprec = TRUE) {
    # Check input
    stopifnot(all(x >= 0), all(x <= 1),
              p > 0, isWhole(p),
              n > 0, isWhole(n),
              m > 0, isWhole(m))

    # Convert to Chiani's notation
    sC <- p
    mC <- 0.5*(abs(n - p) - 1)
    nC <- 0.5*(abs(m - p) - 1)

    doubleWishart_raw(x, sC, mC, nC, mprec)
}

#' @useDynLib rootWishart
#' @importFrom Rcpp sourceCpp
NULL
