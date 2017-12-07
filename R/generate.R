# Functions to generate trajectories

#' Generate a random trajectory
#'
#' Generates a trajectory - either a random search (corresponding to an animal
#' navigating without a compass) or a directed or oriented path (corresponding
#' to an animal navigating with a compass).
#'
#' @param n Number of steps in the trajectory.
#' @param random If TRUE, a random search trajectory is returned, otherwise a
#'   directory trajectory is returned.
#' @param stepLength Mean length of each step in the trajectory, in arbitrary
#'   length units.
#' @param angularErrorSd Standard deviation of angular errors in radians.
#' @param linearErrorSd Standard deviation of linear step length errors.
#' @param fps Simulated frames-per-second - used to generate times for each
#'   point in the trajectory.
#'
#' @return A new Trajectory.
#'
#' @references
#'
#' Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007). Animal
#' navigation: the difficulty of moving in a straight line. Biological
#' Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' @export
TrajGenerate <- function(n = 1000, random = FALSE, stepLength = 2, angularErrorSd = 0.5, linearErrorSd = 0.2, fps = 50) {
  angularErrors <- stats::rnorm(n, sd = angularErrorSd)
  linearErrors <- stats::rnorm(n, sd = linearErrorSd)
  steps <- complex(length.out = n, modulus = stepLength + linearErrors, argument = angularErrors)

  if (random) {
    # Angular errors accumulate
    coords <- complex(n + 1)
    angle <- 0
    for (i in 1:n) {
      angle <- angle + angularErrors[i]
      length <- stepLength + linearErrors[i]
      coords[i + 1] <- coords[i] + complex(modulus = length, argument = angle)
    }
  } else {
    # Angular error resets at every step
    coords <- c(complex(length.out = 1), cumsum(steps))
  }

  TrajFromCoords(data.frame(x = Re(coords), y = Im(coords)), fps = fps)
}
