library(trajr)

context("Basic tests")

trjFromAnglesAndLengths <- function(angles, lengths) {
  coords <- c(0, cumsum(complex(modulus = lengths, argument = angles)))
  TrajFromCoords(data.frame(Re(coords), Im(coords)))
}

# Reads a set of points from a file. The points come from multiple tracks
# due to noise in the video conversion process.
# The longest track is the one we are interested in
#
# Value - data frame with values x & y, and an attribute "numFrames" which records the number of frames in the source video
.MreadPoints <- function(file, ...) {
  points <- read.csv(file, comment.char = '#')

  # Save the number of frames in the file in case the track doesn't extend until the end
  maxFrame <- max(points$Frame)

  # Save number of frames
  attr(points, 'numFrames') <- maxFrame

  points
}

test_that("Trajectory creation", {
  csvFile <- "../testdata/096xypts.csv"
  expect_true(file.exists(csvFile))
  coords <- utils::read.csv(csvFile, stringsAsFactors = FALSE)
  expect_false(is.null(coords))
  trj <- TrajFromCoords(coords, fps = 850)

  expect_false(is.null(trj))
  expect_equal(2030, nrow(trj))
  xRange <- c(997.31, 1541.549436)
  expect_equal(range(trj$x), xRange)
  yRange <- c(669.883810, 956.924828)
  expect_equal(range(trj$y), yRange)
  expect_equal(TrajGetFPS(trj), 850)
  expect_equal(TrajGetNCoords(trj), nrow(coords))

  # Scaling
  scale <- 1 / 2500
  scaled <- TrajScale(trj, scale, "m")
  #plot(scaled)
  expect_false(is.null(scaled))
  expect_equal(nrow(trj), nrow(scaled))
  expect_equal(range(scaled$x), xRange * scale)
  expect_equal(range(scaled$y), yRange * scale)

  # Duration
  expect_equal(TrajDuration(trj), (nrow(trj) - 1) / 850)
  # Velocity
  v <- TrajMeanVelocity(scaled)
  expect_equal(Mod(v), 0.0728938)
  expect_equal(Arg(v), 0.03861151)

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
  # Check it can be plotted without an error
  expect_error(plot(corr, type = 'l'), NA)
  mn <- TrajDAFindFirstMinimum(corr, 10)
  # points(mn["deltaS"], mn["C"], pch = 16, col = "red", lwd = 2)
  # points(mn["deltaS"], mn["C"], col = "black", lwd = 2)
  mx <- TrajDAFindFirstMaximum(corr, 5)
  # points(mx["deltaS"], mx["C"], pch = 16, col = "green", lwd = 2)
  # points(mx["deltaS"], mx["C"], col = "black", lwd = 2)

})

test_that("Creation ignores unimportant NAs", {

  # Silently ignore NAs in other columns
  trj <- TrajFromCoords(data.frame(0, 0, NA))
  expect_equal(nrow(trj), 1)

  # Silently ignore leading or trailing NAs
  trj <- TrajFromCoords(data.frame(c(0, 1, 2, 3, NA), c(NA, 1, 2, 3, 4), c(NA, "a", NA, "c", NA)))
  expect_equal(nrow(trj), 3)

  # Report error if NA is in the middle of the trajectory
  expect_error(TrajFromCoords(data.frame(c(0, 1, 2), c(0, NA, 2), c(NA, 1, NA))),
               "Trajectory contains missing coordinate or time values, first row with NA is 2")

  # Complain about NA in time
  expect_error(TrajFromCoords(data.frame(c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 4), c(NA, "a", NA, "c", NA)), timeCol = 3),
               "Trajectory contains missing coordinate or time values, first row with NA is 3")
  # Skip leading/trailing NAs in time
  trj <- TrajFromCoords(data.frame(c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 4), c(NA, 1, 2, 3, NA)), timeCol = 3)
  expect_equal(nrow(trj), 3)
  trj <- TrajFromCoords(data.frame(c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 4)), timeCol = 3)
  expect_equal(nrow(trj), 5)
})



test_that("Speed intervals", {

  # 1 Interval with no start and 1 stop
  set.seed(1)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_error(plot(intervals), NA)
  expect_equal(nrow(intervals), 1)

  # 1 Interval with 1 start and no stop
  set.seed(2)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_equal(nrow(intervals), 1)

  # 0 intervals
  set.seed(3)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 200
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_error(plot(intervals), NA)
  expect_equal(nrow(intervals), 0)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 150
  fasterThan = 90
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_equal(nrow(intervals), 3)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 50
  fasterThan = NULL
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_error(plot(intervals), NA)
  expect_equal(nrow(intervals), 3)

  # 2 intervals
  set.seed(4)
  trj <- TrajGenerate(20, random = TRUE)
  slowerThan = 92
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 2)

  # Interval wholly contained within a segment
  set.seed(4)
  trj <- TrajGenerate(10, random = TRUE)
  slowerThan = 110
  fasterThan = 107
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 0)

  set.seed(1)
  trj <- TrajGenerate(10, random = TRUE)
  slowerThan = NULL
  fasterThan = 110
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 2)

  # Central diffs gives lower speed
  intervals <- TrajSpeedIntervals(trj, diff = "central", slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 0)

  slowerThan = 107
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 3)

  # Central diffs gives lower speed
  intervals <- TrajSpeedIntervals(trj, diff = "central", slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 2)

  # Entire trajectory is a single interval
  slowerThan = NULL
  fasterThan = 50
  intervals <- TrajSpeedIntervals(trj, diff = "central", slowerThan = slowerThan, fasterThan = fasterThan)
  expect_equal(nrow(intervals), 1)

  # No intervals
  slowerThan = 50
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, diff = "central", slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_equal(nrow(intervals), 0)
})

test_that("Emax", {
  set.seed(1)
  trj1 <- TrajGenerate(1000, angularErrorSd = .5)
  trj2 <- TrajGenerate(1000, angularErrorSd = .2)
  # trj2 (black) should be straighter than trj1 (red), hence Emax(trj1) < Emax(trj2)
  #plot(trj2, asp = NULL, xlim = range(c(trj1$x, trj2$x)), ylim = range(c(trj1$y, trj2$y)))
  #plot(trj1, col = "red", add = TRUE)

  expect_true(TrajEmax(trj1) < TrajEmax(trj2))
})

test_that("Sinuosity", {
  set.seed(1)
  trj1 <- TrajGenerate(1000, angularErrorSd = .5)
  trj2 <- TrajGenerate(1000, angularErrorSd = .2)
  # trj2 (black) should be straighter than trj1 (red), hence Sinuosity(trj1) > Sinuosity(trj2)
  #plot(trj2, asp = NULL, xlim = range(c(trj1$x, trj2$x)), ylim = range(c(trj1$y, trj2$y)))
  #plot(trj1, col = "red", add = TRUE)

  expect_true(TrajSinuosity(trj1) > TrajSinuosity(trj2))
})

test_that("Directional change", {

  # Test that directional change as implemented gives the same results as the equation in the book
  L <- c(1, 1, 1, 1)
  A <- c(pi / 4, 0, -pi / 8, pi / 6)
  trj <- trjFromAnglesAndLengths(A, L)
  #plot(trj, turning.angles = "random")

  .bookCalc <- function(trj) {
    # Lengths between consecutive points
    lengths1 <- Mod(diff(trj$polar))
    # Lengths between points 2 apart
    lengths2 <- Mod(trj$polar[3:nrow(trj)] - trj$polar[1:(nrow(trj) - 2)])
    # Times between points 2 apart
    times2 <- trj$displacementTime[3:nrow(trj)] - trj$displacementTime[1:(nrow(trj) - 2)]

    sapply(1:(nrow(trj) - 2), function(i) {
      a <- lengths1[i]
      b <- lengths1[i+1]
      c <- lengths2[i]
      t <- times2[i]
      (180 - (180 / pi * acos((a ^ 2 + b ^ 2 - c ^ 2) / (2 * a * b)))) / t
    })
  }
  expect_equal(TrajDirectionalChange(trj), .bookCalc(trj))

  set.seed(1)
  trj <- TrajGenerate()
  expect_equal(TrajDirectionalChange(trj), .bookCalc(trj))

  # csvFile <- "../testdata/test-dc.tsv"
  # expect_true(file.exists(csvFile))
  # data <- read.table(csvFile)
  # names(data) <- c("time(s)", "x", "y", "immobile")
  # trj <- TrajFromCoords(data, xCol = "x", yCol = "y", timeCol = "time(s)", spatialUnits = "pixels")
  # expect_equal(TrajDirectionalChange(trj), .bookCalc(trj))

  #microbenchmark(TrajDirectionalChange(trj), .bookCalc(trj), times = 1000)
})

test_that("Reverse", {
  set.seed(1)
  trj <- TrajGenerate()
  rv <- TrajReverse(trj)
  expect_equal(nrow(rv), nrow(trj))
  expect_equal(rv$polar[1], trj$polar[nrow(trj)])
  expect_equal(rv$polar[nrow(rv)], trj$polar[1])
  expect_equal(TrajLength(rv), TrajLength(trj))
  expect_equal(TrajEmax(rv), TrajEmax(trj))
})

test_that("Translate", {
  set.seed(1)
  trj <- TrajGenerate()
  dx <- 10
  dy <- 15
  dt <- 2
  tt <- TrajTranslate(trj, dx, dy, dt)
  expect_equal(nrow(tt), nrow(trj))
  expect_equal(tt$x, trj$x + dx)
  expect_equal(tt$y, trj$y + dy)
  expect_equal(tt$time, trj$time + dt)
  expect_equal(tt$displacement, trj$displacement)
  expect_equal(TrajLength(tt), TrajLength(trj))
  expect_equal(TrajEmax(tt), TrajEmax(trj))

  tto <- TrajTranslate(tt, -tt$x[1], -tt$y[1], -tt$time[1])
  expect_equal(nrow(tto), nrow(trj))
  expect_equal(tto$polar, trj$polar)
  expect_equal(tto$x, trj$x)
  expect_equal(tto$y, trj$y)
  expect_equal(tto$time, trj$time)
  expect_equal(TrajLength(tto), TrajLength(trj))
  expect_equal(TrajEmax(tto), TrajEmax(trj))
})

test_that("Step lengths", {
  set.seed(1)
  nSteps <- 100
  nTrajs <- 4
  stepLength <- 1
  trjs <- lapply(1:nTrajs, TrajGenerate, n = nSteps, stepLength = stepLength)
  sl <- TrajsStepLengths(trjs)
  expect_equal(length(sl), nSteps * nTrajs)
  # Expect mean and median to be roughly equal to the specified step length
  expect_equal(mean(sl), stepLength, tolerance = 2e-2)
  expect_equal(median(sl), stepLength, tolerance = 2e-2)
})

test_that("Generate", {

  unifDist <- function(n) runif(n, -1, 1)

  set.seed(1)
  sd <- 0.5
  trj <- TrajGenerate(angularErrorSd = sd, linearErrorDist = unifDist)
  # Should NOT be able to reject the NULL hypothesis that turning angle errors are normally distributed
  expect_true(shapiro.test(TrajAngles(trj))$p.value > 0.05)
  expect_equal(sd(TrajAngles(trj)), sd, tolerance = 5e-2)
  # Should be able to reject the NULL hypothesis that linear errors are normally distributed
  expect_true(shapiro.test(TrajStepLengths(trj))$p.value <= 0.05)
  trj <- TrajGenerate(angularErrorDist = unifDist)
  # Should be able to reject the NULL hypothesis that turning angles are normally distributed
  expect_true(shapiro.test(TrajAngles(trj))$p.value <= 0.05)
  # Should NOT be able to reject the NULL hypothesis that linear errors are normally distributed
  expect_true(shapiro.test(TrajStepLengths(trj))$p.value > 0.05)
})

test_that("Smoothing", {
  set.seed(1)
  sd <- 0.5
  trj <- TrajGenerate(angularErrorSd = sd)
  smoothed <- TrajSmoothSG(trj, 3, 41)
  expect_true(TrajEmax(trj) < TrajEmax(smoothed))
  smoothed2 <- TrajSmoothSG(trj, 3, 101)
  expect_true(TrajEmax(smoothed) < TrajEmax(smoothed2))
})

test_that("Convenience", {
  tracks <- rbind(
    data.frame(file = "3527.csv", species = "Zodariid2 sp1", category = "spider"),
    data.frame(file = "3530.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3534.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3537.csv", species = "Myrmarachne erythrocephala", category = "mimic spider"),
    data.frame(file = NA, species = "", category = ""),
    data.frame(file = "3542.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3543.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3548.csv", species = "Crematogaster sp1", category = "ant"),
    data.frame(file = NA, species = "", category = ""),
    stringsAsFactors = FALSE
  )
  csvStruct <- list(x = "x", y = "y", time = "Time")

  # Expect to fail with a message when there are blank
  expect_error(TrajsBuild(tracks$file, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints),
               "Trajectory input file name is blank or NULL.*")
  tracks <- na.omit(tracks)
  trjs <- TrajsBuild(tracks$file, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints)

  expect_equal(length(trjs), nrow(tracks))
  expect_equal(TrajGetUnits(trjs[[2]]), "m")
  expect_equal(TrajGetTimeUnits(trjs[[2]]), "s")

  # Trajectories should start at origin
  expect_true(!any(sapply(trjs, function (t) c(t$x[1], t$y[1])) == 0))

  # Define a function which calculates some statistics
  # of interest for a single trajectory
  characteriseTrajectory <- function(trj) {
    # Measures of speed
    derivs <- TrajDerivatives(trj)
    mean_speed <- mean(derivs$speed)
    sd_speed <- sd(derivs$speed)

    # Measures of straightness
    sinuosity <- TrajSinuosity(trj)
    resampled <- TrajRediscretize(trj, .001)
    Emax <- TrajEmax(resampled)

    # Periodicity
    corr <- TrajDirectionAutocorrelations(resampled, round(nrow(resampled) / 4))
    first_min <- TrajDAFindFirstMinimum(corr)

    # Return a list with all of the statistics for this trajectory
    list(mean_speed = mean_speed,
         sd_speed = sd_speed,
         sinuosity = sinuosity,
         Emax = Emax,
         first_min_deltaS = first_min[1],
         first_min_C = first_min[2])
  }

  stats <- TrajsMergeStats(trjs, characteriseTrajectory)

  expect_true(any(is.na(stats)))
  stats <- TrajsStatsReplaceNAs(stats, "first_min_deltaS", flagColumn = "No_first_min")
  stats <- TrajsStatsReplaceNAs(stats, "first_min_C")
  expect_false(any(is.na(stats)))

  # Test translating to the origin
  trjs <- TrajsBuild(tracks$file, translateToOrigin = TRUE, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints)
  expect_true(all(sapply(trjs, function (t) c(t$x[1], t$y[1])) == 0))

})

test_that("Convenience-multi", {

  # Test building multiple trajectories from each "file"
  readFn <- function(filename, ...) {
    # Return 2 very different trajectories
    t1 <- TrajGenerate(50)
    t2 <- TrajGenerate(20, random = FALSE)
    list(t1 = t1[, c('x', 'y', 'time')], t2 = t2[, c('x', 'y', 'time')])
  }

  trjs <- TrajsBuild(c("one", "two"), csvReadFn = readFn, smoothN = 11)

  expect_equal(length(trjs), 4)
})

test_that("works with readr", {
  tracks <- rbind(
    data.frame(file = "3527.csv", species = "Zodariid2 sp1", category = "spider"),
    data.frame(file = "3530.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3534.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3537.csv", species = "Myrmarachne erythrocephala", category = "mimic spider"),
    data.frame(file = "3542.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3543.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3548.csv", species = "Crematogaster sp1", category = "ant"),
    stringsAsFactors = FALSE
  )
  csvStruct <- list(x = "x", y = "y", time = "Time")

  trjs <- TrajsBuild(tracks$file, csvStruct = csvStruct, smoothN = 11, rootDir = "..",
                     csvReadFn = readr::read_csv,
                     col_types = readr::cols(
                       Frame = readr::col_integer(),
                       Time = readr::col_double(),
                       TrackId = readr::col_integer(),
                       x = readr::col_double(),
                       y = readr::col_double(),
                       ValueChanged = readr::col_logical()
                     ))

  expect_equal(length(trjs), nrow(tracks))
})

test_that("Sinuosity", {
  set.seed(1)
  for(aa in seq(0, 2, by = .1)) {
    trj <- TrajGenerate(angularErrorSd = aa)
    # Don't expect equal, just close
    expect_equal(TrajSinuosity(trj), TrajSinuosity2(trj), tolerance = 0.2)
  }
})

test_that("fractal dimension", {
  set.seed(1)
  n <- 5
  angErr <- runif(n, 0, pi)
  trjs <- lapply(1:n, function(i) TrajGenerate(500, angularErrorSd = angErr[i]))
  range <- TrajLogSequence(1, 10, 10)
  fd <- sapply(trjs, function(trj) TrajFractalDimension(trj, range))

  # Test slope of regression
  l <- lm(fd ~ angErr)
  l$coefficients[2] %>% expect_gt(0.2) %>% expect_lt(2)
  })

test_that("Expected square displacement", {
  set.seed(1)
  n <- 200
  angErr <- runif(n, 0, pi)
  trjs <- lapply(1:n, function(i) TrajGenerate(500, angularErrorSd = angErr[i]))
  esd1 <- sapply(trjs, function(trj) TrajExpectedSquareDisplacement(trj, eqn1 = TRUE))
  esd2 <- sapply(trjs, function(trj) TrajExpectedSquareDisplacement(trj, eqn1 = FALSE))

  # plot(angErr, y = abs(esd1), log = 'xy', pch = 16, cex = .7)
  # points(angErr, y = esd2, pch = 16, cex = .6, col = "red")

  # Test slopes of regressions
  l <- lm(log(esd1) ~ log(angErr))
  l$coefficients[2] %>% expect_lt(-1.5) %>% expect_gt(-2)

  l <- lm(log(esd2) ~ log(angErr))
  l$coefficients[2] %>% expect_lt(-1.5) %>% expect_gt(-2)
})

test_that("straightness r", {
  set.seed(1)
  n <- 200
  angErr <- runif(n, 0, pi)
  trjs <- lapply(1:n, function(i) TrajRediscretize(TrajGenerate(500, angularErrorSd = angErr[i]), 2))
  sir <- sapply(trjs, function(trj) Mod(TrajMeanVectorOfTurningAngles(trj)))
  sid <- sapply(trjs, function(trj) TrajStraightness(trj))

  # plot(angErr, y = sid, pch = 16, cex = .7, ylim = range(c(sir, sid)))
  # points(angErr, y = sir, pch = 16, cex = .6, col = "red")

  # Test slopes of regressions
  l <- lm(sir ~ angErr)
  l$coefficients[2] %>% expect_lt(0) %>% expect_gt(-.5)

  l <- lm(sid ~ angErr)
  l$coefficients[2] %>% expect_lt(0) %>% expect_gt(-.5)

})

test_that("plots", {
  csvFile <- "../testdata/096xypts.csv"
  coords <- utils::read.csv(csvFile, stringsAsFactors = FALSE)
  trj <- TrajFromCoords(coords, fps = 850)

  # Scaling
  scale <- 1 / 2500
  scaled <- TrajScale(trj, scale, "m")

  # Smoothing
  smoothed <- TrajSmoothSG(scaled, 3, 101)


  # Expect no errors from plotting (weird syntax!)
  expect_error(plot(scaled), NA)
  expect_error(lines(smoothed, col = "red"), NA)
  expect_error(points(smoothed, pch = '.', col = 'green'), NA)

  # Plot a simple trajectory with turning angles
  set.seed(2)
  trj <- TrajGenerate(5)
  expect_error(plot(scaled, turning.angles = TRUE))
  expect_error(plot(trj, turning.angles = "random"), NA)
  expect_error(plot(trj, turning.angles = "directed"), NA)
})

test_that("rotation", {
  # Expect a value to be equal in the first and last points
  .expectSameCol <- function(trj, col = "x") {
    expect_equal(trj[1, col], tail(trj, 1)[, col])
  }

  set.seed(1)
  trj <- TrajGenerate(10)
  # All default parameters
  r <- TrajRotate(trj)
  .expectSameCol(r, "y")
  rotAngle <- pi / 2
  r <- TrajRotate(trj, rotAngle)
  .expectSameCol(r)
  # Test rotation of trajectory that doesn't start at the origin
  ttrj <- TrajTranslate(trj, 1, 1)
  tr <- TrajRotate(ttrj, rotAngle)
  .expectSameCol(tr)
  or <- TrajRotate(trj, rotAngle, origin = as.numeric(trj[1, c("x", "y")]))
  .expectSameCol(or)
  # Rotate about start point
  otr <- TrajRotate(ttrj, rotAngle, origin = as.numeric(ttrj[1, c("x", "y")]))
  .expectSameCol(otr)
  # Rotated start point should be same as unrotated start point
  expect_equal(otr[1, "x"], ttrj[1, "x"])
  expect_equal(otr[1, "y"], ttrj[1, "y"])

  vo <- TrajMeanVelocity(trj)
  vr <- TrajMeanVelocity(r)

  # Expect mean vector length and path length to be unchanged, but angle to be changed
  expect_equal(Mod(vr), Mod(vo))
  expect_equal(TrajLength(trj), TrajLength(r))
  expect_true(Arg(vr) != Arg(vo))

  # Test absolute rotation
  aotr <- TrajRotate(ttrj, rotAngle, origin = as.numeric(ttrj[1, c("x", "y")]), relative = FALSE)
  expect_equal(Arg(aotr$displacement[2]), Arg(ttrj$displacement[2]) + rotAngle)
  expect_true(aotr[1, "y"] != tail(aotr, 1)[, "y"])

  plot(trj, ylim = c(-5, 20))
  plot(ttrj, lty = 2, add = TRUE)
  plot(r, col = "blue", add = TRUE)
  plot(tr, col = "blue", lty = 2, add = TRUE)
  plot(or, col = "red", lty = 3, add = TRUE)
  plot(otr, col = "red", lty = 3, add = TRUE)
  plot(aotr, col = "green", lty = 3, lwd = 2, add = TRUE)
})

test_that("Convert times", {
  .checkTimes <- function(t, s) {
    seconds <- TrajConvertTime(t)
    expect_equal(length(t), length(seconds))
    expect_identical(seconds, s)
  }

  .checkTimes(c("0:00:00:000", "0:00:00:001", "0:00:01:000", "0:01:00:000", "1:00:00:000"),
              c(0, .001, 1, 60, 60 * 60))

  .checkTimes(c("0:00:00:001", "0:00:00:002", "0:00:00:003", "0:00:00:004"),
              c(.001, .002, .003, .004))

  .checkTimes(c("1:01:01:001", "2:02:02:002", "3:03:03:003", "4:04:04:004"),
              c(3661.001, 7322.002, 10983.003, 14644.004))
})

test_that("Resampling", {
  # Plot one trajectory over another
  plotTwoTrjs <- function(trj1, trj2) {
    plot(trj1, draw.start.pt = FALSE, lwd = 2)
    points(trj1, cex = .6, draw.start.pt = FALSE)
    lines(trj2, col = "red", lty = 2)
    points(trj2, col = "red", draw.start.pt = FALSE)
  }

  set.seed(1)
  # Give it a constant step length of 1 so that it has a constant speed to simplify some tests
  trj <- TrajGenerate(10, angularErrorSd = 1, stepLength = 1, linearErrorDist = function(n) rep(0, n), fps = 1)
  trjL <- TrajLength(trj)

  # These tests aren't strictly required to be true because resampled
  # trajectories are also smoothed, so may be shorter than expected
  ta <- TrajResampleTime(trj, .5)
  expect_true(trjL - TrajLength(ta) < .5)
  # plotTwoTrjs(trj, ta)
  tb <- TrajResampleTime(trj, .7)
  expect_true(trjL - TrajLength(tb) < .7)
  # plotTwoTrjs(trj, tb)
  tc <- TrajResampleTime(trj, 1)
  expect_true(trjL - TrajLength(tc) == 0)
  # plotTwoTrjs(trj, tc)
  td <- TrajResampleTime(trj, 2)
  expect_true(trjL - TrajLength(td) < 2)
  # plotTwoTrjs(trj, td)
  te <- TrajResampleTime(trj, 2.1)
  # This isn't TRUE, since resampled te is straighter than trj
  #expect_true(trjL - TrajLength(td) < 2)
  # plotTwoTrjs(trj, te)
})

test_that("Invalid parameter detection", {
  expect_error(TrajsBuild("short.csv", csvStruct = list(x = "x", y = "y", time="t"), rootDir = ".."),
               "Invalid smoothing parameter n \\(41): n must be less than the number of points in the trajectory \\(5)")
})

test_that("FPS calculation", {
  tracks <- rbind(
    data.frame(file = "3527.csv", species = "Zodariid2 sp1", category = "spider"),
    stringsAsFactors = FALSE
  )
  csvStruct <- list(x = "x", y = "y")

  # Expect an error if there's no time column or FPS specified
  expect_error(TrajsBuild(tracks$file, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints),
               ".*Cannot create a trajectory without times: one of fps or a time column must be specified")
  fps <- 50
  trjs <- TrajsBuild(tracks$file, fps = fps, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints)
  # Expect the time interval between every 50 frames to be 1 second (allow small tolerance when testing == 1)
  diffs <- diff(trjs[[1]]$Time, fps)
  expect_true(all(abs(diffs - 1) < .00001))

  # Now do the same thing with a frame rate of 240
  fps <- 240
  trjs <- TrajsBuild(tracks$file, fps = fps, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints)
  # Expect the time interval between every 50 frames to be 1 second (allow small tolerance when testing == 1)
  diffs <- diff(trjs[[1]]$Time, fps)
  expect_true(all(abs(diffs - 1) < .00001))
})

test_that("Empty trajectory", {
  trj <- TrajFromCoords(data.frame(numeric(), numeric()))
  expect_equal(nrow(trj), 0)
  trj <- TrajFromCoords(data.frame(0, 0))
  expect_equal(nrow(trj), 1)
})

test_that("Turning angles", {
  set.seed(1)
  nsteps <- 10000
  trj <- TrajGenerate(nsteps)
  expect_equal(nrow(trj), nsteps + 1)
  expect_equal(length(TrajAngles(trj)), nsteps - 1)
  expect_equal(length(TrajAngles(trj, compass.direction = 0)), nsteps)

  # # Now add in some 0-length segments
  idx <- sort(sample(seq_len(nsteps + 1), round(nsteps * 1.5), replace = T))
  trj <- TrajFromCoords(trj[idx, ])
  nsteps <- nrow(trj) - 1
  expect_equal(nrow(trj), nsteps + 1)
  expect_equal(length(TrajAngles(trj)), nsteps - 1)
  expect_equal(length(TrajAngles(trj, compass.direction = 0)), nsteps)

})
