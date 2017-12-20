# Functions to generate trajectories

#' Generate a random trajectory
#'
#' Generates a trajectory. If \code{random} is \code{TRUE}, the trajectory will
#' be a correllated random walk/idiothetic directed walk (Kareiva & Shigesada,
#' 1983), corresponding to an animal navigating without a compass (Cheung,
#' Zhang, Stricker, & Srinivasan, 2008). If \code{random} is \code{FALSE}, it
#' will be a directed walk/allothetic directed walk/oriented path, corresponding
#' to an animal navigating with a compass (Cheung, Zhang, Stricker, &
#' Srinivasan, 2007, 2008).
#'
#' By default, for both random and directed walks, errors are normally
#' distributed, unbiased, and independent of each other, so are \emph{simple
#' directed walks} in the terminology of Cheung, Zhang, Stricker, & Srinivasan,
#' (2008). This behaviour may be modified by specifying alternative values for
#' the \code{angularErrorDist} and/or \code{linearErrorDist} parameters.
#'
#' The initial angle (for a random walk) or the intended direction (for a
#' directed walk) is \code{0} radians. The starting position is \code{(0, 0)}.
#'
#' @param n Number of steps in the trajectory.
#' @param random If TRUE, a random search trajectory is returned, otherwise a
#'   directed trajectory (with direction = 0 radians) is returned.
#' @param stepLength Mean length of each step in the trajectory, in arbitrary
#'   length units.
#' @param angularErrorSd Standard deviation of angular errors in radians.
#' @param angularErrorDist Function which accepts a single argument - the number
#'   of values to return, and generates random deviates according to some
#'   distribution. The returned values are added to the previous step angle
#'   (when \code{random == TRUE}), or to \code{0} (is \code{random == FALSE}) to
#'   generate the step angle for each step in the trajectory. If the mean of the
#'   returned values is not zero, the walk will be biased.
#' @param linearErrorSd Standard deviation of linear step length errors.
#' @param linearErrorDist Function which accepts a single argument - the number
#'   of values to return, and generates random deviates according to some
#'   distribution. The returned values are added to \code{stepLength} to
#'   generate the lengths of each step.
#' @param fps Simulated frames-per-second - used to generate times for each
#'   point in the trajectory.
#'
#' @return A new Trajectory with \code{n} segments and \code{n + 1} coordinate
#'   pairs.
#'
#' @examples
#' # Generate a 1000 step correlated random walk
#' trj <- TrajGenerate()
#' plot(trj, main = "Correlated walk")
#'
#' # Generate a 1000 step levy flight - paths lengths follow a cauchy distribution
#' trj <- TrajGenerate(linearErrorDist = rcauchy)
#' plot(trj, main = "Levy flight")
#'
#' # Generate a short directed trajectory
#' trj <- TrajGenerate(n = 20, random = FALSE)
#' plot(trj, main = "Directed walk")
#'
#' # Generate an uncorrelated random walk
#' trj <- TrajGenerate(500, angularErrorDist = function(n) runif(n, -pi, pi))
#' plot(trj, main = "Uncorrelated walk")
#'
#' @references
#'
#' Kareiva, P. M., & Shigesada, N. (1983). Analyzing insect movement as a
#' correlated random walk. Oecologia, 56(2), 234-238. doi:10.1007/bf00379695
#'
#' Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007). Animal
#' navigation: the difficulty of moving in a straight line. Biological
#' Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2008). Animal
#' navigation: general properties of directed walks. Biological Cybernetics,
#' 99(3), 197-217. doi:10.1007/s00422-008-0251-z
#'
#' @export
TrajGenerate <- function(n = 1000, random = TRUE, stepLength = 2,
                         angularErrorSd = 0.5,
                         angularErrorDist = function(n) stats::rnorm(n, sd = angularErrorSd),
                         linearErrorSd = 0.2,
                         linearErrorDist = function(n) stats::rnorm(n, sd = linearErrorSd),
                         fps = 50) {
  angularErrors <- angularErrorDist(n)
  linearErrors <- linearErrorDist(n)
  stepLengths <- stepLength + linearErrors
  # Don't allow negative lengths
  stepLengths[stepLengths < 0] <- 0
  steps <- complex(length.out = n, modulus = stepLengths, argument = angularErrors)

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
    # Angular error resets at every step. First coordinate is at the origin
    coords <- c(complex(length.out = 1), cumsum(steps))
  }

  TrajFromCoords(data.frame(x = Re(coords), y = Im(coords)), fps = fps)
}
