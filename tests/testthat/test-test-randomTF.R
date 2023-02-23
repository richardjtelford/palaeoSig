test_that("randomTF works", {
  library(rioja)
  data("RLGH")
  data("SWAP")
  expect_no_condition(randomTF(spp = SWAP$spec, fos = RLGH$spec, env = SWAP$pH, fun = WA, col = "WA.inv", n = 10))
  expect_no_condition(randomTF(spp = SWAP$spec, fos = RLGH$spec, env = SWAP$pH, fun = MAT, col = "MAT.wm", n = 10, k = 5))
})
