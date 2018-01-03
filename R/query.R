# Trajectory basic query functions

# ---- Trajectory query ----

#' Trajectory frames-per-second
#'
#' Returns the frames-per-second recorded for this trajectory.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetFPS <- function(trj) { attr(trj, .TRAJ_FPS) }

#' Trajectory number of coordinates
#'
#' Returns the number of coordinates recorded for this trajectory, i.e. 1 more than the number of steps.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetNCoords <- function(trj) { attr(trj, .TRAJ_NFRAMES) }

#' Trajectory spatial units
#'
#' Returns the spatial units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajScale}}, \code{\link{TrajGetTimeUnits}}.
#'
#' @export
TrajGetUnits <- function(trj) { attr(trj, .TRAJ_UNITS) }

#' Trajectory temporal units
#'
#' Returns the temporal units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajFromCoords}}, \code{\link{TrajGetUnits}}.
#'
#' @export
TrajGetTimeUnits <- function(trj) {
  attr(trj, .TRAJ_TIME_UNITS)
}

# ---- Trajectory analysis ----

#' Trajectory step lengths
#'
#' Returns the lengths of each step in a trajectory.
#'
#' @param trj Trajectory to query.
#'
#' @seealso \code{\link{TrajLength}}
#'
#' @export
TrajStepLengths <- function(trj) {
  Mod(utils::tail(trj$displacement, -1)) # First displacement is not a step (also usually 0)
}

#' Trajectory distance
#'
#' Calculates the distance between the start and end of a trajectory (or a
#' portion of a trajectory). Also called the diffusion distance, net distance,
#' displacement, or bee-line from start to finish.
#'
#' @param trj Trajectory whose distance is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#' @return Numeric distance from the start to the end of the trajectory.
#'
#' @export
TrajDistance <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  Mod(diff(trj$polar[c(startIndex, endIndex)]))
}

#' Trajectory length
#'
#' Calculates the cumulative length of a trajectory (or a portion of a
#' trajectory), which is the total distance travelled along the trajectory.
#'
#' @param trj Trajectory whose length is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#' @return Path Numeric length of the trajectory.
#'
#' @seealso \code{\link{TrajStepLengths}}
#'
#' @export
TrajLength <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  sum(Mod(diff(trj$polar[startIndex:endIndex])))
}

#' Turning angles of a Trajectory
#'
#' Calculates the step angles (in radians) of each segment, either relative to
#' the previous segment or relative to the specified compass direction.
#'
#' @param trj the trajectory whose whose angles are to be calculated.
#' @param lag Angles between every lag'th segment are calculated.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @return Step angles in radians, normalised so that \code{-pi < angle <= pi}.
#'
#' @seealso \code{\link{TrajStepLengths}}
#'
#' @export
TrajAngles <- function(trj, lag = 1, compass.direction = NULL) {
  if (is.null(compass.direction)) {
    angles <- diff(Arg(trj$displacement[2:nrow(trj)]), lag)
  } else {
    angles <- Arg(trj$displacement[2:nrow(trj)]) - compass.direction
  }
  # Normalise so that -pi < angle <= pi
  ii <- angles <= -pi
  angles[ii] <- angles[ii] + 2 * pi
  ii <- angles > pi
  angles[ii] <- angles[ii] - 2 * pi
  angles
}

#' Trajectory expected square displacement
#'
#' Calculates the expected square displacement for a trajectory assuming it is a
#' correlated random walk, using the formula in Kareiva & Shigesada, (1983).
#'
#' Note that Cheung, Zhang, Stricker, and Srinivasan (2007) define an
#' alternative formulation for expected maximum displacement, Emax (see
#' \code{\link{TrajEmax}}).
#'
#' @param trj A Trajectory.
#' @param n Number of steps to calculate.
#' @param eqn1 If \code{TRUE}, calculate using equation 1, otherwise using
#'   equation 2. Equation 2 applies when the mean of turning angles is 0,
#'   i.e.turns are unbiased.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @examples
#' set.seed(1)
#' # A random walk
#' trj <- TrajGenerate(200)
#' smoothed <- TrajSmoothSG(trj)
#'
#' # Calculate actual squared displacement at all points along the trajectory
#' sd2 <- sapply(2:nrow(smoothed), function(n) TrajDistance(smoothed, 1, n) ^ 2)
#' # Calculate expected squared displacement
#' ed2_1 <- sapply(2:nrow(smoothed), function(n) TrajExpectedSquareDisplacement(smoothed, n, TRUE))
#' ed2_2 <- sapply(2:nrow(smoothed), function(n) TrajExpectedSquareDisplacement(smoothed, n, FALSE))
#'
#' # Plot expected against actual. According to Kareiva & Shigesada, (1983), if actual
#' # (approximately) matches expected, the trajectory is probably a correlated random walk
#' par(mar = c(5, 5, 0.1, 0.1) + .1)
#' plot(2:nrow(smoothed), sd2, type = 'l', pch = 16, cex = .2, lwd = 2,
#'      xlab = 'Number of consecutive moves',
#'      ylab = expression('Squared displacement, ' * R[n]^2))
#' lines(2:nrow(smoothed), ed2_1, col = "grey", lwd = 2)
#' lines(2:nrow(smoothed), ed2_2, col = "pink", lwd = 2)
#'
#' legend("bottomright",
#'        c(expression("Actual displacement"^2),
#'          expression("Expected displacement"^2 * " (eqn 1)"),
#'          expression("Expected displacement"^2 * " (eqn 2)")),
#'        col = c('black', 'grey', 'pink'), lwd = 2,
#'        inset = c(0.01, 0.02))
#'
#' @seealso \code{\link{TrajEmax}}
#'
#' @references
#'
#' Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007). Animal
#' navigation: the difficulty of moving in a straight line. Biological
#' Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' Kareiva, P. M., & Shigesada, N. (1983). Analyzing insect movement as a
#' correlated random walk. Oecologia, 56(2), 234-238. doi:10.1007/bf00379695
#'
#' @export
TrajExpectedSquareDisplacement <- function(trj, n = nrow(trj), eqn1 = TRUE, compass.direction = NULL) {
  sl <- TrajStepLengths(trj)
  ta <- TrajAngles(trj, compass.direction = compass.direction)
  l <- mean(sl)
  l2 <- mean(sl ^ 2)
  c <- mean(cos(ta))
  s <- mean(sin(ta))
  s2 <- s^2

  if (eqn1) {
    # Eqn 1
    alpha <- atan2(s, c)
    gamma <- ((1 - c)^2 - s2) * cos((n + 1) * alpha) - 2 * s * (1 - c) * sin((n + 1) * alpha)
    n * l2 + 2 * l^2 * ((c - c^2 - s2) * n  - c) / ((1 - c)^2 + s2) +
      2 * l^2 * ((2 * s2 + (c + s2) ^ ((n + 1) / 2)) / ((1 - c)^2 + s2)^2) * gamma
  } else {
    # Eqn 2
    n * l2 + 2 * l^2 * c / (1 - c) * (n - (1 - c^n) / (1 - c))
  }
}
