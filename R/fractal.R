# Fractal index functions

#' Fractal dimension step sizes
#'
#' Calculates a vector of possible fractal dimension step sizes, suitable for
#' use when calculating the fractal dimension of a trajectory. Note that the
#' selected step sizes may determine the calculated fractal dimension, which may
#' render fractal dimension meaningless for animal trajectories.
#'
#' By default, returns a vector of step sizes ranging from 1/1000th of the path
#' length to 1/10th of the path length, with logarithmically increasing intervals.
#'
#' @param trj Trajectory to return step sizes for.
#' @param smallestStepFactor Multiplied by the path length of \code{trj} to
#'   determine the smallest step size in the returned range.
#' @param largestStepFactor Multiplied by the path length of \code{trj} to
#'   determine the largest step size in the returned range.
#' @param numSteps The number of step sizes in the returned vector.
#' @return A vector of step sizes which can be passed to
#'   \code{\link{TrajFractalDimension}}.
#'
#' @seealso \code{\link{TrajFractalDimension}}
#'
#' @export
TrajFractalStepSizes <- function(trj, smallestStepFactor = .0001, largestStepFactor = .1, numSteps = 100) {
  pathLen <- TrajLength(trj)
  exp(log(10) * seq(log10(smallestStepFactor * pathLen), log10(largestStepFactor * pathLen), length.out = numSteps))
}

#' Fractal dimension calculation
#'
#' Calculates path length (\eqn{L(\delta)}) for a range of step sizes (\eqn{\delta}).
#' @param trj Trajectory to calculate fractal dimension for.
#' @param stepSizes Vector of step sizes used to calculate path lengths (see \code{\link{TrajFractalStepSizes}}).
#' @return Data frame with columns \code{stepsize} (\eqn{\delta}) and \code{pathlength} ((\eqn{L(\delta)}).
#'
#' @seealso \code{\link{TrajFractalStepSizes}} for a possible calculation of \code{stepSizes}.
#' \code{\link{TrajFractalDimension}} for fractal dimension calculation.
#'
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate()
#' plot(TrajFractalDimensionValues(trj, TrajFractalStepSizes(trj)), log = "xy", pch = 16, cex = .5)
#'
#' @export
TrajFractalDimensionValues <- function(trj, stepSizes) {
  fi <- data.frame(t(sapply(stepSizes, function(ss) c(ss, TrajLength(TrajRediscretize(trj, ss))))))
  names(fi) <- c('stepsize', 'pathlength')
  fi
}

#' Fractal dimension of a trajectory
#'
#' Calculates the fractal dimension of a trajectory using the 'dividers' method
#' (Sugihara & May, 1990). This value may be meaningless for animal trajectories
#' as they may not be fractal curves - see Benhamou (2004) and Turchin (1996).
#'
#' @param trj Trajectory to calculate fractal dimension for.
#' @param stepSizes Vector of step sizes used to calculate path lengths (see \code{\link{TrajFractalStepSizes}}).
#' @return The fractal dimension of the trajectory for the given step sizes.
#' @seealso \code{\link{TrajFractalStepSizes}}, \code{\link{TrajFractalDimensionValues}}.
#'
#' @references Sugihara, G., & M. May, R. (1990). Applications of fractals in ecology. Trends in Ecology & Evolution, 5(3), 79-86. doi:10.1016/0169-5347(90)90235-6
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's path. Journal of Theoretical Biology, 229(2), 209-220. doi:10.1016/j.jtbi.2004.03.016
#'
#' Turchin, P. (1996). Fractal Analyses of Animal Movement: A Critique. Ecology, 77(7), 2086-2090. doi:10.2307/2265702
#'
#' @export
TrajFractalDimension <- function(trj, stepSizes) {
  # Calculate path length for a range of step sizes
  fi <- TrajFractalDimensionValues(trj, stepSizes)
  # Calculate slope of points on log/log axes
  l <- stats::lm(log(fi$pathlength) ~ log(fi$stepsize))
  slope <- l$coefficients[2]
  # Fractal dimension
  unname(1 - slope)
}
