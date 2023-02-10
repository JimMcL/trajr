#' Resample a trajectory to a constant time interval.
#'
#' Constructs a new trajectory by resampling the input trajectory to a fixed
#' time interval. Points are linearly interpolated along the trajectory. Spatial
#' and time units are preserved.
#'
#' @param trj The trajectory to be resampled.
#' @param stepTime The resampled trajectory step time. Each step in the new
#'   trajectory will have this duration.
#' @param newFps Value to be stored as the FPS value in the new trajectory (see
#'   \code{\link{TrajGetFPS}}). It is not otherwise used by this function.
#' @return A new trajectory with a constant time interval for each step. Points
#'   in the new trajectory are calculated by linearly interpolating along
#'   \code{trj}.
#'
#' @examples
#' # Simulate a trajectory with steps every 5 hours
#' set.seed(46)
#' trj <- TrajGenerate(10, stepLength = 5, fps = 1/5, timeUnits = "hours", linearErrorSd = .8)
#'
#' # Resample to 1 hour steps
#' resampled <- TrajResampleTime(trj, 1)
#'
#' par(mar = c(5, 4, .5, .5))
#' plot(trj, lwd = 2)
#' points(trj, pch = 16)
#' points(resampled, col = "red", draw.start.pt = FALSE)
#'
#' @export
TrajResampleTime <- function(trj, stepTime, newFps = NULL) {
  if (stepTime <= 0) {
    stop("stepTime must be > 0")
  }
  # Determine times of new points
  times <- seq(from = min(trj$time), to = max(trj$time), by = stepTime)

  # Interpolate x and y separately. This works by treating x (or y) as a
  # function of time (so x axis is time, y axis is trajectory x or y), then
  # interpolating to required times. approx interpolates y for given values of x
  x <- stats::approx(trj$time, trj$x, times)$y
  y <- stats::approx(trj$time, trj$y, times)$y

  # Create the new trajectory, preserving metadata
  TrajFromCoords(data.frame(x, y, times), timeCol = 3,
                 fps = newFps,
                 spatialUnits = TrajGetUnits(trj),
                 timeUnits = TrajGetTimeUnits(trj))
}
