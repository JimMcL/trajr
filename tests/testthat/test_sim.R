library(trajr)

context("Simulation tests")

test_that("Split and merge", {
  trj <- TrajGenerate(n = 100)
  l <- TrajSplit(trj, c(12, 13, 30))
  expect_equal(length(l), 4)
  nt <- TrajMerge(l)
  expect_true(all(trj == nt))

  l <- TrajSplit(trj, c(30, 13, 12))
  expect_equal(length(l), 4)
  nt <- TrajMerge(l)
  expect_true(all(trj == nt))

  l <- TrajSplit(trj, c(-1, 0, 1))
  expect_equal(length(l), 1)
  expect_true(all(trj == TrajMerge(l)))

  l <- TrajSplit(trj, c(-1, 0, 1, nrow(trj) + 1))
  expect_equal(length(l), 1)
  expect_true(all(trj == TrajMerge(l)))

  l <- TrajSplit(trj, c(-1, 0, 1, nrow(trj) + 1, 50, 50))
  expect_equal(length(l), 2)
  expect_true(all(trj == TrajMerge(l)))
})

test_that("Boundary handling", {
  set.seed(1)
  trj <- TrajGenerate(n = 8)
  # Square arena
  boundary <- data.frame(x = c(-10, 10, 10, -10), y = c(-10, -10, 10, 10))

  inside <- TrajInPolygon(trj, boundary)
  expect_equal(inside, c(1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L))

  l <- TrajSplitAtFirstCrossing(trj, boundary)
  expect_equal(length(l), 2)
  expect_equal(TrajInPolygon(l[[1]], boundary), rep(1, 6))
  expect_equal(TrajInPolygon(l[[2]], boundary), rep(0, 3))
  expect_true(all(trj == TrajMerge(l)))
})
