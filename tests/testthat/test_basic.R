library(trajr)

context("trajectory creation")

test_that("Trajectory creation", {
  csvFile <- "../testdata/096xypts.csv"
  expect_true(file.exists(csvFile))
  coords <- read.csv(csvFile, stringsAsFactors = FALSE)
  expect_false(is.null(coords))
  trj <- TrajFromCoords(coords, fps = 1000)
  expect_false(is.null(trj))
  expect_equal(2030, nrow(trj))
  xRange <- c(997.31, 1541.549436)
  expect_equal(range(trj$x), xRange)
  yRange <- c(669.883810, 956.924828)
  expect_equal(range(trj$y), yRange)

  # Scaling
  scale <- .1
  scaled <- TrajScale(trj, scale, "m")
  #TrajPlot(scaled)
  expect_false(is.null(scaled))
  expect_equal(nrow(trj), nrow(scaled))
  expect_equal(range(scaled$x), xRange * scale)
  expect_equal(range(scaled$y), yRange * scale)

  # Smoothing
  smoothed <- TrajSmoothSG(scaled, 3, 101)
  #TrajPlot(smoothed)
  expect_true(TrajLength(smoothed) < TrajLength(scaled))
  expect_true(abs(TrajDistance(smoothed) - TrajDistance(scaled)) < TrajDistance(scaled) / 10)

  # Derivatives
  derivs <- TrajDerivatives(smoothed)
  plot(derivs$speed, type = 'l', col = 'red')
  plot(derivs$acceleration, type = 'l')

  # Rediscretization
  rd <- TrajRediscretize(smoothed, .05)
  #TrajPlot(rd)

  expect_true(TrajStraightness(smoothed) < 1)
  expect_true(TrajStraightness(smoothed) > 0)

  corr <- TrajDirectionAutocorrelations(TrajRediscretize(smoothed, .5))
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
