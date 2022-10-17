library("testthat")
library("usethis")


usethis::use_testthat()
# use tset_file / test_dir to run tests

test_file("tests/testthat/test-metadata_functions.R")
