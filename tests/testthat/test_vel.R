library(trajr)

context("Velocity and acceleration tests")

test_that("Basic velocity", {
  # Crude test that it doesn't crash and returns the expected number of values
  set.seed(1)
  n <- 10
  # Argument to TrajGenerate is number of steps, whereas we want n to be number
  # of points, i.e. 1 more than number of steps
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

test_that("velocity", {
  # Test returned values are as expected.
  # Simple trajectory, constant velocity in x, exponentially increasing in y
  x <- 1:10
  y <- x ^ 2
  trj <- TrajFromCoords(data.frame(x, y), fps = 1)
  vel <- TrajVelocity(trj)
  # Discard first and last elements - we know they are NA (tested above)
  vel <- vel[2:(length(vel) - 1)]
  # X velocities should all be 1
  expect_true(all(Re(vel) == 1))
  # Y velocities should be increasing linearly
  expect_true(all(Im(vel) == seq(4, length.out = length(vel), by = 2)))
})

test_that("Basic acceleration", {
  set.seed(1)
  n <- 10
  # Argument to TrajGenerate is number of steps, whereas we want n to be number
  # of points, i.e. 1 more than number of steps
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

test_that("Acceleration", {
  # Test returned values are as expected.
  # Simple trajectory, constant velocity in x, exponentially increasing in y
  x <- 1:10
  y <- x ^ 2
  trj <- TrajFromCoords(data.frame(x, y), fps = 1)
  acc <- TrajAcceleration(trj)
  # Discard first and last elements - we know they are NA (tested above)
  acc <- acc[2:(length(acc) - 1)]
  # X acceleration should all be 0 (constant velocity)
  expect_true(all(Re(acc) == 0))
  # Y acceleration should all be 2 (constant but non-zero acceleration)
  expect_true(all(Im(acc) == 2))
})

test_that("Variable step times", {
  # Simple trajectory, constant step length of 5, speed oscillates between 1 & 5
  x <- 1:10 * 3
  y <- 1:10 * 4
  times <- cumsum(rep(c(1, 5), times = 5))
  trj <- TrajFromCoords(data.frame(x, y, times), timeCol = 3)
  # Calculate velocity using forward differences
  fvel <- TrajVelocity(trj, diff = "forward")
  expect_true(is.na(fvel[length(fvel)]))
  # Discard last element
  vel <- fvel[1:(length(fvel) - 1)]
  # Speeds should match those from TrajDerivates since the trajectory is 1-dimensional and we've used forward differences
  expect_true(all(Mod(vel) == TrajDerivatives(trj)$speed))

  # Backward should be the same as forward except NA is first rather than last element
  bvel <- TrajVelocity(trj, diff = "backward")
  expect_true(is.na(bvel[1]))
  expect_true(all(bvel[2:length(bvel)] == vel))
})

test_that("acc sanity", {
  # Perform a sanity check: simple implementation of acceleration by central
  # differences from fixed step times should yield the same result as
  # TrajAcceleration applied to fixed step times

  set.seed(2)
  n <- 100
  trj <- TrajGenerate(n = n - 1)
  acc <- TrajAcceleration(trj)

  # Recalculate using central diffs
  dt <- mean(diff(trj$time))
  ax <- stats::filter(trj$x, c(1, -2, 1)) / dt ^ 2
  ay <- stats::filter(trj$y, c(1, -2, 1)) / dt ^ 2

  expect_equal(as.numeric(Re(acc)), as.numeric(ax))
  expect_equal(as.numeric(Im(acc)), as.numeric(ay))
})

test_that("Variable steps acc", {
  # Step lengths are equal, but time decreases so trajectory is accelerating
  x <- c(1, 1, 1, 1, 1)
  y <- 1:5 * 8
  times <- cumsum(c(0, 8, 4, 2, 1))
  trj <- TrajFromCoords(data.frame(x, y, times), timeCol = 3)

  # Integrating velocity should give us distance travelled at each step
  vel <- TrajVelocity(trj, diff = "forward")
  vel <- head(vel, -1) # Remove trailing NA
  dx <- Re(vel) * diff(times)
  expect_equal(dx, diff(x))
  dy <- Im(vel) * diff(times)
  expect_equal(dy, diff(y))

  # Compare acceleration to hand-calculated acceleration
  acc <- TrajAcceleration(trj)
  acc <- na.omit(acc)
  expect_equal(as.numeric(Re(acc)), rep(0, 3))
  # Speeds should be 1, 2, 4, 8
  expect_equal(Im(vel), c(1, 2, 4, 8))
  # Acceleration should be (diffs in speed) / (diffs in times of step centres)
  expect_equal(as.numeric(Im(acc)), c(1, 2, 4) / c(6, 3, 1.5))
})

test_that("velocity & speed", {
  set.seed(5)
  trj <- TrajGenerate()
  # Compare speed from TrajDerivatives & TrajVelocity
  drv <- TrajDerivatives(trj)
  vel <- TrajVelocity(trj, diff = "backward")
  vspeed <- data.frame(speed = Mod(vel), speedTimes = trj$time)
  vspeed <- na.omit(vspeed)
  expect_equal(drv$speed, vspeed$speed)
  expect_equal(drv$speedTimes, vspeed$speedTimes)

  # Try a different, log trajectory
  trj <- TrajGenerate(n = 10000, linearErrorDist = rcauchy)
  drv <- TrajDerivatives(trj)
  vel <- TrajVelocity(trj, diff = "backward")
  vspeed <- na.omit(data.frame(speed = Mod(vel), speedTimes = trj$time))
  expect_equal(drv$speed, vspeed$speed)
  expect_equal(drv$speedTimes, vspeed$speedTimes)

  trj <- TrajGenerate(50000, angularErrorDist = function(n) runif(n, -pi, pi))
  drv <- TrajDerivatives(trj)
  vel <- TrajVelocity(trj, diff = "backward")
  vspeed <- na.omit(data.frame(speed = Mod(vel), speedTimes = trj$time))
  expect_equal(drv$speed, vspeed$speed)
  expect_equal(drv$speedTimes, vspeed$speedTimes)

})
