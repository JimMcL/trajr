library(trajr)

context("trajectory creation")

test_that("Trajectory creation", {
  csvFile <- "../testdata/096xypts.csv"
  expect_true(file.exists(csvFile))
  coords <- utils::read.csv(csvFile, stringsAsFactors = FALSE)
  expect_false(is.null(coords))
  trj <- TrajFromCoords(coords, fps = 1000)

  expect_false(is.null(trj))
  expect_equal(2030, nrow(trj))
  xRange <- c(997.31, 1541.549436)
  expect_equal(range(trj$x), xRange)
  yRange <- c(669.883810, 956.924828)
  expect_equal(range(trj$y), yRange)

  # Scaling
  scale <- 1 / 2500
  scaled <- TrajScale(trj, scale, "m")
  #plot(scaled)
  expect_false(is.null(scaled))
  expect_equal(nrow(trj), nrow(scaled))
  expect_equal(range(scaled$x), xRange * scale)
  expect_equal(range(scaled$y), yRange * scale)

  # Smoothing
  smoothed <- TrajSmoothSG(scaled, 3, 101)
  #plot(smoothed)
  expect_true(TrajLength(smoothed) < TrajLength(scaled))
  expect_true(abs(TrajDistance(smoothed) - TrajDistance(scaled)) < TrajDistance(scaled) / 10)

  # Derivatives
  derivs <- TrajDerivatives(smoothed)
  #plot(derivs$speed, type = 'l', col = 'red')
  #plot(derivs$acceleration, type = 'l')

  # Rediscretization
  rd <- TrajRediscretize(smoothed, .001)
  #plot(rd)

  expect_true(TrajStraightness(smoothed) < 1)
  expect_true(TrajStraightness(smoothed) > 0)

  corr <- TrajDirectionAutocorrelations(rd)
  # plot(corr, type='l')
  mn <- TrajDAFindFirstMinimum(corr, 10)
  # points(mn["deltaS"], mn["C"], pch = 16, col = "red", lwd = 2)
  # points(mn["deltaS"], mn["C"], col = "black", lwd = 2)
  mx <- TrajDAFindFirstMaximum(corr, 5)
  # points(mx["deltaS"], mx["C"], pch = 16, col = "green", lwd = 2)
  # points(mx["deltaS"], mx["C"], col = "black", lwd = 2)

  sinuosity <- TrajSinuosity(rd)

  emax <- TrajEmax(smoothed)
})

test_that("Speed intervals", {
  plotIntervalsByTime <- function(smoothed, slowerThan, fasterThan, intervals) {
    derivs <- TrajDerivatives(smoothed)
    speed <- derivs$speed
    plot(x = derivs$speedTimes, y = speed, type = 'l', xlab = 'Time (sec)', ylab = "Speed")
    abline(h = slowerThan, col = "red")
    abline(h = fasterThan, col = "green")
    if (nrow(intervals) > 0)
      rect(intervals$startTime, min(speed), intervals$stopTime, max(speed), col = "#0000FF1E", border = NA)
  }

  plotIntervalsByFrame <- function(smoothed, slowerThan, fasterThan, intervals) {
    derivs <- TrajDerivatives(smoothed)
    speed <- derivs$speed
    plot(speed, type = 'l', xlab = 'Frame number', ylab = "Speed")
    abline(h = slowerThan, col = "red")
    abline(h = fasterThan, col = "green")
    if (nrow(intervals) > 0)
      rect(intervals$startFrame, min(speed), intervals$stopFrame, max(speed), col = "#0000FF1E", border = NA)
  }

  # 1 Interval with no start and 1 stop
  set.seed(1)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(smoothed, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 1)

  # 1 Interval with 1 start and no stop
  set.seed(2)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(smoothed, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 1)

  # 0 intervals
  set.seed(3)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 200
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(smoothed, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 0)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 150
  fasterThan = 90
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(smoothed, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 3)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 50
  fasterThan = NULL
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(smoothed, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 3)

  # 2 intervals
  set.seed(4)
  trj <- TrajGenerate(20, random = TRUE)
  slowerThan = 92
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(trj, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 2)

  # Interval wholly contained within a segment
  set.seed(4)
  trj <- TrajGenerate(10, random = TRUE)
  slowerThan = 110
  fasterThan = 107
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plotIntervalsByTime(trj, slowerThan, fasterThan, intervals)
  expect_true(nrow(intervals) == 0)
})

test_that("Emax", {
  set.seed(1)
  trj1 <- TrajGenerate(1000, angularErrorSd = .5)
  trj2 <- TrajGenerate(1000, angularErrorSd = .2)
  # trj2 (black) should be straighter than trj1 (red)
  plot(trj2, asp = NULL, xlim = range(c(trj1$x, trj2$x)), ylim = range(c(trj1$y, trj2$y)))
  plot(trj1, col = "red", add = TRUE)

  TrajEmax(trj1)
  TrajEmax(trj2)
})
