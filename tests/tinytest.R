if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.0.0") {
    tinytest::test_package("rootWishart",
                           ncpu = getOption("Ncpus", 1))
}
