#' Resample a trajectory to a constant time interval.
#'
#' Constructs a new trajectory by resampling the input trajectory to a fixed
#' time interval. Points are linearly interpolated along the trajectory. Spatial
#' and time units are preserved.
#'
#' @param trj The trajectory to be resampled.
#' @param stepTime New sampled step time. Each step in the new trajectory will
#'   have this duration.
#' @param newFps Value to be stored as the FPS value in the new trajectory (see
#'   \code{\link{TrajGetFPS}}). It is not otherwise used by this function.
#' @return A new trajectory with a constant time interval for each step. Points
#'   in the new trajectory are calculated by linearly interpolating along
#'   \code{trj}.
#'
#' @export
TrajResampleTime <- function(trj, stepTime, newFps = NULL) {
  # Determine times of new points
  times <- seq(from = min(trj$time), to = max(trj$time), by = stepTime)
  # Interpolate x and y separately
  x <- stats::approx(trj$time, trj$x, times)$y
  y <- stats::approx(trj$time, trj$y, times)$y

  TrajFromCoords(data.frame(x, y, times), timeCol = 3,
                 fps = newFps,
                 spatialUnits = TrajGetUnits(trj),
                 timeUnits = TrajGetTimeUnits(trj))
}
