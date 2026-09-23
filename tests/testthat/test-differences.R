test_that("common difference infrastructure implements q = 1 and q = 2", {
  D1 <- TRSC:::.trsc_difference_matrix(5, 1)
  D2 <- TRSC:::.trsc_difference_matrix(5, 2)

  expect_equal(D1 %*% (1:5), matrix(rep(1, 4), ncol = 1))
  expect_equal(D2 %*% (1:5), matrix(rep(0, 3), ncol = 1))

  quadratic <- (1:5)^2
  expect_equal(drop(D2 %*% quadratic), rep(2, 3))
})

test_that("movement inside a duplicate-donor null representation has zero energy", {
  X <- cbind(
    A = c(1, 2, 3, 4),
    B = c(0, 1, 0, 1),
    A_duplicate = c(1, 2, 3, 4)
  )
  g <- trsc_geometry(X)
  W <- rbind(
    c(1, 0.75, 0.50, 0.25, 0),
    c(0, 0, 0, 0, 0),
    c(0, 0.25, 0.50, 0.75, 1)
  )
  expect_equal(TRSC:::.trsc_temporal_energy(W, g, 1), 0, tolerance = 1e-10)
  expect_equal(TRSC:::.trsc_temporal_energy(W, g, 2), 0, tolerance = 1e-10)
})
