context("3D tests")

# Hack together a random 3D path
gen3dPts <- function(n = 100) {
  pts <- TrajGenerate(n)
  dz <- rnorm(n + 1)
  pts$z[1] <- 10
  for (i in 2:(nrow(pts))) {
    pts$z[i] <- pts$z[i - 1] + dz[i]
  }
  pts
}

#library(rgl)
plot3dTrj <- function(t3, add = FALSE, col = 1, points = FALSE) {
  plot3d(t3$x, t3$y, t3$z, add = add, col = col, type = 'l', lwd = 2)
  if (points)
    plot3d(t3$x, t3$y, t3$z, add = T, col = col, type = 's', size = 0.6)
  else
    plot3d(t3$x[1], t3$y[1], t3$z[1], add = T, col = col, type = 's', size = 0.8)
}

test_that("3D basic", {
  track <- data.frame(x = c(0, 1, 1, 1),
                      y = c(0, 0, 1, 1),
                      z = c(0, 0, 0, 1))
  t3 <- Traj3DFromCoords(track)
  expect_equal(Traj3DLength(t3), 3)
  expect_equal(Traj3DDistance(t3), sqrt(3))
  expect_equal(Traj3DStraightness(t3), sqrt(3) / 3)
  expect_equal(Traj3DStepLengths(t3), rep(1, nrow(track) - 1))

  t32 <- Traj3DFromCoords(track + 0.6)
  expect_equal(Traj3DLength(t32), 3)
  expect_equal(Traj3DDistance(t32), sqrt(3))
  expect_equal(Traj3DStraightness(t32), sqrt(3) / 3)
  expect_equal(Traj3DStepLengths(t32), rep(1, nrow(track) - 1))

  set.seed(1)
  pts <- gen3dPts(100)
  tr <- Traj3DFromCoords(pts, zCol = "z", timeCol = "time")
  expect_gt(Traj3DDistance(tr), TrajDistance(tr))
  expect_gt(Traj3DLength(tr), TrajLength(tr))
  expect_gt(Traj3DLength(tr), Traj3DDistance(tr))

  sm <- Traj3DSmoothSG(tr)
  # A smoothed trajectory should be shorter
  expect_lt(Traj3DLength(sm), Traj3DLength(tr))
  # Expect start and end points to be close to each other, but not identical
  expect_equal(Traj3DDistance(sm), Traj3DDistance(tr), tolerance = .05)

  rs <- Traj3DResampleTime(sm, 0.05)
  # A resampled trajectory (with longer step times) should be even shorter
  expect_lt(Traj3DLength(rs), Traj3DLength(sm))
})

test_that("3d rediscretization", {
  track <- data.frame(x = c(0, 1, 1, 1),
                      y = c(0, 0, 1, 1),
                      z = c(0, 0, 0, 1))
  t3 <- Traj3DFromCoords(track)
  r <- Traj3DRediscretize(t3, 1)
  expect_true(inherits(r, .TRAJ_3D_CLASS))
  expect_equal(r[, c("x", "y", "z")], t3[, c("x", "y", "z")])

  r <- Traj3DRediscretize(t3, 0.5)
  expect_true(inherits(r, .TRAJ_3D_CLASS))
  expect_equal(nrow(r), 7)
  expect_equal(TrajGetUnits(r), TrajGetUnits(t3))
  expect_equal(Traj3DLength(r), Traj3DLength(t3))
  expect_equal(Traj3DDistance(r), Traj3DDistance(t3))
  expect_equal(Traj3DStepLengths(r), rep(0.5, 6))

  r <- Traj3DRediscretize(t3, 0.8)
  expect_true(inherits(r, .TRAJ_3D_CLASS))
  expect_equal(nrow(r), 4)
  expect_equal(TrajGetUnits(r), TrajGetUnits(t3))
  expect_lt(Traj3DLength(r), Traj3DLength(t3))
  expect_lt(Traj3DDistance(r), Traj3DDistance(t3))
  expect_equal(Traj3DStepLengths(r), rep(0.8, 3))

  r <- Traj3DRediscretize(t3, 1.1)
  expect_true(inherits(r, .TRAJ_3D_CLASS))
  expect_equal(nrow(r), 3)
  expect_equal(TrajGetUnits(r), TrajGetUnits(t3))
  expect_lt(Traj3DLength(r), Traj3DLength(t3))
  expect_lt(Traj3DDistance(r), Traj3DDistance(t3))
  expect_equal(Traj3DStepLengths(r), rep(1.1, 2))

  set.seed(1)
  tr <- gen3dPts(100)
  r <- Traj3DRediscretize(tr, 1.1)
  expect_equal(Traj3DStepLengths(r), rep(1.1, nrow(r) - 1))
  expect_lt(Traj3DLength(r), Traj3DLength(tr))
  # Can't say much about distance because it depends on which direction the two
  # ends are facing, but it shouldn't change by more than 1 step length
  expect_lt(abs(Traj3DDistance(r) - Traj3DDistance(tr)), 1.1)
  #plot3dTrj(tr); plot3dTrj(r, add = T, col = 2, points = T)
})

test_that("rediscretize with speed", {
  pts <- gen3dPts(100)
  tr <- Traj3DFromCoords(pts, zCol = "z", timeCol = "time")

  rd <- Traj3DRediscretize(tr, 2, simConstantSpeed = TRUE)
  #plot3dTrj(tr); plot3dTrj(rd, add = T, col = 2, points = T)
  # Start times should be equal
  expect_equal(rd$time[1], tr$time[1])
  # Average speed should be similar
  rdSp <- Traj3DLength(rd) / TrajDuration(rd)
  trjSp <- Traj3DLength(tr) / TrajDuration(tr)
  expect_lt(abs(log(rdSp / trjSp)), log(1.02))

  # Test that simulation without time throws exception
  rd2 <- Traj3DRediscretize(tr, 2, simConstantSpeed = FALSE)
  expect_error(Traj3DRediscretize(rd2, 4, simConstantSpeed = TRUE))

})

##### TODO
# test_that("Mean vector", {
#   track <- data.frame(x = c(0, 1, 1),
#                       y = c(0, 0, 1),
#                       z = c(0, 0, 0))
#   t3 <- Traj3DFromCoords(track)
#   mv <- Traj3DMeanVectorOfTurningAngles(t3)
#   mv2 <- TrajMeanVectorOfTurningAngles(t3)
#   expect_equal(norm(as.matrix(mv), type = "2"), Mod(mv2))
#
#   track <- data.frame(x = c(0, 1, 1, 1, 3),
#                       y = c(0, 0, 1, 2, 4),
#                       z = c(0, 0, 0, 0, 0))
#   t3 <- Traj3DFromCoords(track)
#   mv <- Traj3DMeanVectorOfTurningAngles(t3)
#   mv2 <- TrajMeanVectorOfTurningAngles(t3)
#   expect_equal(norm(as.matrix(mv), type = "2"), Mod(mv2))
#
#   track <- data.frame(x = c(0, 1, 1, 1),
#                       y = c(0, 0, 1, 1),
#                       z = c(0, 0, 0, 1))
#   t3 <- Traj3DFromCoords(track)
#   mv <- Traj3DMeanVectorOfTurningAngles(t3)
# })

