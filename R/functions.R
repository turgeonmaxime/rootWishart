#' Incomplete beta function
#'
#' Computes the incomplete beta function.
#'
#' This is defined as the integral from 0 to \code{x} of the function
#' \code{t^(a-1)(1-t)^(b-1)}.
#'
#' @param x Value between 0 and 1 at which to evaluate the incomplete beta
#'   function.
#' @param a non-negative numeric
#' @param b non-negative numeric
#' @return Returns the value of the incomplete beta function at x.
#' @export
incompleteBeta <- function(x, a, b){
    if(any(c(length(x), length(a), length(b)) != 1)) {
        stop("None of the parameters can be vectors of length > 1", call. = FALSE)
    }
    stopifnot(x >= 0, x <= 1, a >= 0, b >= 0)

    # Special cases
    if(a == 0 || b == 0) {
        val <- Inf
    } else {
        if(x == 0) val <- 0
        if(x == 1) val <- beta(a, b)
        if(x > 0 && x < 1) val <- incompleteBeta_C(x, alpha = a, beta = b)
    }

    return(val)
}

#' Distribution of the largest root
#'
#' Computes the cumulative distribution function of the largest root in the single and double Wishart setting.
#'
#' @param x Value at which to compute the CDF.
#' @param s,m,n Parameters of the double Wishart setting. See Details.
#' @return Returns the value of the CDF at \code{x}.
#' @export
#' @aliases doubleWishart singleWishart
#' @rdname largestRoot
doubleWishart <- function(x, s, m, n) {
    doubleWishart_C(x, s, m, n)
}

#' @param n_min,n_max Parameters of the single Wishart setting. See details.
#' @export
#' @rdname largestRoot
singleWishart <- function(x, n_min, n_max) {
    singleWishart_C(x, n_min, n_max)
}
