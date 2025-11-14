test_that("singlearm function works", {
  expect_no_error(singlearmDesign(nmin = 50,
                                  nmax = 70,
                                  C = 10,
                                  p0 = 0.05,
                                  p1 = 0.15,
                                  power = 0.8,
                                  alpha = 0.05,
                                  progressBar = FALSE))
})
