# Fractal index functions


#' Logarithmically spaced sequence
#'
#' Convenience function to return a sequence of points which are regularly
#' spaced when plotted on a logarithmic axis.
#'
#' @param from Starting value of the sequence.
#' @param to End (maximal) value of the sequence.
#' @param length.out Desired length of the sequence (non-negative). Rounded up if fractional.
#'
#' @seealso \code{\link[base]{seq}}
#'
#' @export
TrajLogSequence <- function(from, to, length.out) {
  exp(log(10) * seq(log10(from), log10(to), length.out = length.out))
}

#' Fractal dimension calculation
#'
#' Calculates path length (\eqn{L(\delta)}) for a range of step sizes
#' (\eqn{\delta}). For a fractal (i.e. scale independent) curve,
#' \eqn{log(L(\delta))} grows linearly as \eqn{log(\delta)} grows smaller. In
#' other words, if the points returned by this function lie on a straight line
#' in a log-log plot, \code{trj} is a fractal curve.
#'
#' @param trj Trajectory to calculate fractal dimension for.
#' @param stepSizes Vector of step sizes used to calculate path lengths.
#' @param adjustD If TRUE, path length is adjusted to reduce truncation error
#'   (Nams, 2006).
#' @return Data frame with columns \code{stepsize} (\eqn{\delta}) and
#'   \code{pathlength} ((\eqn{L(\delta)}).
#'
#' @seealso \code{\link{TrajFractalDimension}} for fractal dimension
#'   calculation.
#'
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate()
#' muL <- mean(TrajStepLengths(trj))
#' # Use 20 step sizes from 1/2 mean step length to 5 * mean step length.
#' # For real use, biologically meaningful step sizes should be used.
#' stepSizes <- TrajLogSequence(0.5 * muL, 5 * muL, 20)
#' plot(TrajFractalDimensionValues(trj, stepSizes), log = "xy", pch = 16, cex = .5)
#'
#' @references
#'
#' Nams, V. O. (2006). Improving Accuracy and Precision in Estimating Fractal
#' Dimension of Animal movement paths. Acta Biotheoretica, 54(1), 1-11.
#' doi:10.1007/s10441-006-5954-8
#'
#' @export
TrajFractalDimensionValues <- function(trj, stepSizes, adjustD = TRUE) {

  .redLen <- function(ss) {
    red <- TrajRediscretize(trj, ss)
    len <- TrajLength(red)
    if (adjustD) {
      # Add the straight line distance from the end of the rediscretized
      # trajectory to the end of the original trajectory
      lastRed <- red$polar[length(red$polar)]
      lastTrj <- trj$polar[length(trj$polar)]
      len <- len + Mod(lastRed - lastTrj)
    }
    c(ss, len)
  }
  fi <- data.frame(t(sapply(stepSizes, .redLen)))
  names(fi) <- c('stepsize', 'pathlength')
  fi
}

#' Fractal dimension of a trajectory
#'
#' Calculates the fractal dimension (\code{D}) of a trajectory using the
#' 'dividers' method (Sugihara & May, 1990). By default, overestimation of
#' \code{D} is compensated for as recommended by Nams (2006), by walking the
#' dividers backwards and forwards, and by estimating the remaining path length
#' at the end of the last step.
#'
#' Fractal dimension may be meaningless for animal trajectories as they may not
#' be true fractal curves - see Benhamou (2004) and Turchin (1996), although it
#' may be useful for studies involving differences in behaviour at different
#' spatial scales (Nams, 2006).
#'
#' You can test whether a trajectory is a fractal curve for a range of step
#' sizes using the \code{\link{TrajFractalDimensionValues}} function. The
#' example code in its documentation demonstrates how to plot path length for a
#' range of step sizes. If the plotted points lie along straight line, then the
#' trajectory is a fractal curve for that range of step sizes. However, typical
#' trajectories result in a curve rather than a straight line.
#'
#' If you decide to use fractal dimension despite the warnings of Benhamou
#' (2004) and Turchin (1996), try to select a biologically meaningful range of
#' step sizes (and be prepared to justify your choice). If comparing fractal
#' dimensions across trajectories, be consistent in your choice of step sizes.
#'
#' @param trj Trajectory to calculate fractal dimension for.
#' @param stepSizes Vector of step sizes (aka divider sizes) used to calculate
#'   path lengths.
#' @param adjustD If \code{TRUE}, path length is adjusted for truncation error
#'   (Nams, 2006).
#' @param dMean If \code{TRUE}, the fractal dimension is calculated starting
#'   from the beginning of the trajectory, then re-calculated starting from the
#'   end and moving backwards. The value returned is the mean of the two fractal
#'   dimensions (Nams, 2006).
#'
#' @return The fractal dimension of the trajectory for the given step sizes.
#'
#' @seealso \code{\link{TrajLogSequence}} to create a logarithmically spaced
#'   sequence, \code{\link{TrajFractalDimensionValues}} for the function used
#'   internally to calculate a range of path lengths for different step sizes,
#'   \code{\link{TrajEmax}} and \code{\link{TrajSinuosity2}} for some alternate
#'   measures of trajectory tortuosity.
#'
#' @references
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's
#' path. Journal of Theoretical Biology, 229(2), 209-220.
#' doi:10.1016/j.jtbi.2004.03.016
#'
#' Nams, V. O. (2006). Improving Accuracy and Precision in Estimating Fractal
#' Dimension of Animal movement paths. Acta Biotheoretica, 54(1), 1-11.
#' doi:10.1007/s10441-006-5954-8
#'
#' Sugihara, G., & M. May, R. (1990). Applications of fractals in ecology.
#' Trends in Ecology & Evolution, 5(3), 79-86. doi:10.1016/0169-5347(90)90235-6
#'
#' Turchin, P. (1996). Fractal Analyses of Animal Movement: A Critique. Ecology,
#' 77(7), 2086-2090. doi:10.2307/2265702
#'
#' @export
TrajFractalDimension <- function(trj, stepSizes, adjustD = TRUE, dMean = TRUE) {

  .calcD <- function(trj) {
    # Calculate path length for a range of step sizes
    vals <- TrajFractalDimensionValues(trj, stepSizes, adjustD = adjustD)
    # Calculate slope of points on log/log axes
    l <- stats::lm(log(pathlength) ~ log(stepsize), data = vals)
    slope <- l$coefficients[2]
    # Fractal dimension
    unname(1 - slope)
  }

  if (dMean)
    # Take mean of fractal D in forward and reverse directions
    mean(c(.calcD(trj), .calcD(TrajReverse(trj))))
  else
    .calcD(trj)
}
