if (
  requireNamespace("tinytest", quietly = TRUE) &&
    identical(Sys.getenv("RUN_TESTS"), "true")
) {
  set.seed(42)
  tinytest::test_package("BVAR", at_home = TRUE, pattern = "^.*\\.[rR]$")
}
