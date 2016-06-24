#' @useDynLib rootWishart
#' @export
incompleteBeta <- function(x, alpha, beta) {
    out <- .C("incompleteBeta_R",
              x, alpha, beta, numeric(1))
    res <- out[[4]]
    return(res)
}
