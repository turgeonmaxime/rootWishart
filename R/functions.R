#' Distribution of the largest root
#'
#' Computes the cumulative distribution function of the largest root in the single and double Wishart setting.
#'
#' @param x Value at which to compute the CDF.
#' @param s,m,n Parameters of the double Wishart setting. See Details.
#' @param mprec Logical. Should we use gmp for high precision computations?
#' @return Returns the value of the CDF at \code{x}.
#' @export
#' @aliases doubleWishart singleWishart
#' @rdname largestRoot
doubleWishart <- function(x, s, m, n, mprec = TRUE) {
    doubleWishart_raw(x, s, m, n, mprec)
}

#' @param n_min,n_max Parameters of the single Wishart setting. See details.
#' @export
#' @rdname largestRoot
singleWishart <- function(x, n_min, n_max, mprec = TRUE) {
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
            result <- singleWishart_raw(x, n_min, n_max)
        }
    }

    return(result)
}

#' @useDynLib rootWishart
#' @importFrom Rcpp sourceCpp
NULL
