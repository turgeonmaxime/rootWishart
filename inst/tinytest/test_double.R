library(rootWishart)
x <- seq(0, 1, length.out = 50)

# no error in double wishart
y1 <- try(doubleWishart(x, p = 10, n = 10, m = 200),
          silent = TRUE)
y2 <- try(doubleWishart(x, p = 10, n = 10, m = 200,
                        type = "fixed"),
          silent = TRUE)
y3 <- try(doubleWishart(x, p = 10, n = 10, m = 200,
                        type = "arbitrary"),
          silent = TRUE)

expect_false(inherits(y1, "try-error"))
expect_false(inherits(y2, "try-error"))
expect_false(inherits(y3, "try-error"))
