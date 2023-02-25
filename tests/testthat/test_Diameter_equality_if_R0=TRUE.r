context("Diameter estimation")
library(TapeR)

test_that("Diameter are interpolated", {
  xm <- 1.3 / 27 # relative measurement height
  ym <- runif(1, 10, 90) # measured diameter
  xp <- xm # measurement height / tree height
  data("SK.par.lme", package = "TapeR")
  EBLUB_R0TRUE <- TapeR:::SK_EBLUP_LME.f(xm, ym, xp, SK.par.lme, 
                                         Rfn = list(fn="zero"))
  expect_equal(EBLUB_R0TRUE$yp[,,drop=TRUE], ym[1])
  expect_output(str(EBLUB_R0TRUE), "List of 9")
})
