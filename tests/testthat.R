# This file is part of the standard testthat setup.
# When R CMD check runs, it sources this file, which then discovers and runs
# every test-*.R script inside tests/testthat/.
#
# See <https://testthat.r-lib.org/reference/test_package.html>.

library(testthat)
library(dosefitr)

test_check("dosefitr")
