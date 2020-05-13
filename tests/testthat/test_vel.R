library(trajr)

context("Velocity and acceleration tests")

test_that("Velocity", {
  set.seed(1)
  n <- 10
  # Argument to Trajgenerate is number of steps, whereas we want n to be number of points
  trj <- TrajGenerate(n = n - 1)
  vel <- TrajVelocity(trj)
  expect_equal(length(vel), n)
  expect_true(inherits(vel, "complex"))
  expect_equal(attr(vel, "trj"), trj)
  # Expect 1st and last to be NA, no others
  expect_true(is.na(vel[1]))
  expect_true(is.na(vel[n]))
  expect_equal(sum(is.na(vel)), 2)
})

test_that("Acceleration", {
  set.seed(1)
  n <- 10
  # Argument to Trajgenerate is number of steps, whereas we want n to be number of points
  trj <- TrajGenerate(n = n - 1)
  acc <- TrajAcceleration(trj)
  expect_equal(length(acc), n)
  expect_true(inherits(acc, "complex"))
  expect_equal(attr(acc, "trj"), trj)
  # Expect 1st and last to be NA, no others
  expect_true(is.na(acc[1]))
  expect_true(is.na(acc[n]))
  expect_equal(sum(is.na(acc)), 2)
})
