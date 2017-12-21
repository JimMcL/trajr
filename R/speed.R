# Functions relating to trajectory speed and acceleration

# Internal function which just throws an error if trj doesn't contain time information
.checkTrajHasTime <- function(trj) {
  if (!("displacementTime" %in% names(trj)))
    stop("Missing time information in trajectory. Perhaps this is a rediscretized trajectory?")
}

#' Calculates trajectory speed and acceleration
#'
#' Calculates speed and linear acceleration along a trajectory over time. Noisy
#' trajectories should be smoothed before being passed to this function, as
#' noise is effectively amplifed when calculating speed and acceleration.
#'
#' @param trj Trajectory whose speed and acceleration is to be calculated.
#'
#' @return A list with components: \item{speed}{numeric vector, speed between
#'   each pair of trajectory points.} \item{speedTimes}{numeric vector, times
#'   corresponding to values in \code{speed}.} \item{acceleration}{numeric
#'   vector.} \item{accelerationTimes}{numeric vector.}
#'
#' @seealso \code{\link{TrajSpeedIntervals}} for analysing intervals within the
#'   trajectory of low or high speed. \code{\link{TrajSmoothSG}} for smoothing a
#'   trajectory.
#'
#' @export
TrajDerivatives <- function(trj) {
  .checkTrajHasTime(trj)

  # Note that displacements are the (polar) displacements from 1 point to the next
  d <- Mod(trj$displacement)
  t <- trj$displacementTime

  # Calculate speed
  v <- d[2:length(d)] / diff(t)
  vt <- t[2:length(t)]
  # Calculate linear acceleration
  a <- diff(v) / diff(vt)
  at <- vt[2:length(vt)]

  list(speed = v, speedTimes = vt, acceleration = a, accelerationTimes = at)
}

# Linear interpolation of interval times for TrajSpeedIntervals
.linearInterpTimes <- function(slowerThan, fasterThan, speed, times, startFrames, startTimes, stopFrames, stopTimes) {

  if(is.null(slowerThan))
    slowerThan <- NA
  if(is.null(fasterThan))
    fasterThan <- NA

  .interp <- function(f) {
    if (f < length(speed)) {
      proportion <- (slowerThan - speed[f]) / (speed[f + 1] - speed[f])
      if (is.na(proportion) || proportion < 0 || proportion > 1)
        proportion <- (fasterThan - speed[f]) / (speed[f + 1] - speed[f])
      if (!is.na(proportion) && proportion >= 0 && proportion <= 1)
        return(times[f] + proportion * (times[f + 1] - times[f]))
    }
    times[f]
  }

  .frameIsInInterval <- function(f) {
    (!is.na(fasterThan) && speed[f] > fasterThan) ||
      (!is.na(slowerThan) && speed[f] < slowerThan)
  }

  for (i in 1:length(startFrames)) {
    # Special case if starting point was added as the start of an interval, leave its time unchanged
    if (!(startFrames[i] == 1 && .frameIsInInterval(1)))
      startTimes[i] <- .interp(startFrames[i])
  }
  for (i in 1:length(stopFrames)) {
    # Special case if stopping point was added as the end of an interval, leave its time unchanged
    numFrames <- length(speed)
    if (!(stopFrames[i] == numFrames && .frameIsInInterval(numFrames)))
      stopTimes[i] <- .interp(stopFrames[i])
  }
  list(startTimes, stopTimes)
}

#' Calculate speed time intervals
#'
#' Calculates and returns a list of time intervals during which speed is slower
#' and/or faster than specified values.
#'
#' @param trj Trajectory to be analysed.
#' @param fasterThan,slowerThan If not \code{NULL}, intervals will cover time
#'   periods where speed exceeds/is lower than this value.
#' @param interpolateTimes If \code{TRUE}, times will be linearly interpolated
#'   between frames.
#'
#' @return A data frame of class "TrajSpeedIntervals", each row is an interval,
#'   columns are: \item{startFrame}{Indices of frames at the start of each
#'   interval.} \item{stopFrame}{Indices of frames at the end of each interval.}
#'   \item{startTime}{Time at the start of each interval.} \item{stopTime}{Time
#'   at the end of each interval} \item{duration}{Duration of each interval.}
#'
#'   The data frame will also have non-standard attributes:
#'   \item{trajectory}{Value of the \code{trj} argument.}
#'   \item{slowerThan}{Value of the \code{slowerThan} argument.}
#'   \item{fasterThan}{Value of the \code{fasterThan} argument.}
#'   \item{derivs}{Value returned by calling \code{TrajDerivatives(trj)}.}
#'
#' @seealso \code{\link{TrajDerivatives}} for calculating trajectory speed and
#'   acceleration, \code{\link{plot.TrajSpeedIntervals}} for plotting speed over
#'   time with intervals highlighted.
#'
#' @examples
#' # Plot speed, highlighting intervals where speed drops below 50 units/sec
#' set.seed(4)
#' trj <- TrajGenerate(200, random = TRUE)
#' smoothed <- TrajSmoothSG(trj, 3, 101)
#' intervals <- TrajSpeedIntervals(smoothed, slowerThan = 50, fasterThan = NULL)
#' plot(intervals)
#'
#' # Report the duration of the maximum period of low speed
#' cat(sprintf("Duration of the longest low-speed interval was %g secs\n", max(intervals$duration)))
#'
#' @export
TrajSpeedIntervals <- function(trj, fasterThan = NULL, slowerThan = NULL, interpolateTimes = TRUE) {
  if (is.null(fasterThan) && is.null(slowerThan)) {
    stop("Parameters fasterThan and slowerThan are both NULL, one must be specified")
  }

  # Calculate trajectory speeds
  derivs <- TrajDerivatives(trj)
  speed <- derivs$speed
  times <- derivs$speedTimes

  # Calculate for each point whether it is within an interval
  flags <- rep(TRUE, length(speed))
  if(!is.null(fasterThan)) {
    flags <- flags & (speed > fasterThan)
  }
  if(!is.null(slowerThan)) {
    flags <- flags & (speed < slowerThan)
  }

  changes <- diff(flags)
  stopFrames <- which(changes == -1)
  startFrames <- which(changes == 1)

  # Handle situation where interval begins or ends outside of trajectory
  if (length(startFrames) > 0 || length(stopFrames) > 0) {
    # Assume interval started at beginning of trajectory, since we don't know what happened before that
    if (length(stopFrames) > 0 && (length(startFrames) == 0 || stopFrames[1] < startFrames[1]))
      startFrames <- c(1, startFrames)
    # Similarly, assume that interval can't extend past end of trajectory
    if (length(stopFrames) == 0 || startFrames[length(startFrames)] > stopFrames[length(stopFrames)])
      stopFrames <- c(stopFrames, length(speed))
  }

  stopTimes <- times[stopFrames]
  startTimes <- times[startFrames]

  # Maybe linearly interpolate times
  if (interpolateTimes && length(startFrames) > 0) {
    r <- .linearInterpTimes(slowerThan, fasterThan, speed, times, startFrames, startTimes, stopFrames, stopTimes)
    startTimes <- r[[1]]
    stopTimes <- r[[2]]
  }

  durations <- stopTimes - startTimes

  result <- data.frame(startFrame = startFrames, startTime = startTimes, stopFrame = stopFrames, stopTime = stopTimes, duration = durations)

  # Record some attributes
  attr(result, "slowerThan") <- slowerThan
  attr(result, "fasterThan") <- fasterThan
  attr(result, "derivs") <- derivs
  attr(result, "trajectory") <- trj

  # Give it a class
  class(result) <- c("TrajSpeedIntervals", class(result))

  result
}

#' Plot method for trajectory speed intervals
#'
#' Plots speed over time, with intervals of fast and/or slow speed highlighted.
#'
#' @param x An object of class "SpeedIntervals", as created by
#'   \code{\link{TrajSpeedIntervals}}.
#' @param slowerThanColour,fasterThanColour The colour of the horizontal line
#'   plotted at the "slower than" or "faster than" speed. Specify \code{NULL} to
#'   prevent the line from being plotted.
#' @param highlightColor Colour of the highlight rectangles.
#' @param xlab,ylab plotting parameters with useful defaults.
#' @param ... Additional arguments are passed to \code{\link[graphics]{plot}}.
#'
#' @seealso \code{\link{TrajSpeedIntervals}}
#'
#' @export
plot.TrajSpeedIntervals <- function(x,
                                    slowerThanColour = "red",
                                    fasterThanColour = "green",
                                    highlightColor = "#0000FF1E",
                                    xlab = sprintf("Time (%s)", TrajGetTimeUnits(attr(x, "trajectory"))),
                                    ylab = sprintf("Speed (%s/%s)", TrajGetUnits(attr(x, "trajectory")), TrajGetTimeUnits(attr(x, "trajectory"))),
                                    ...) {
  derivs <- attr(x, "derivs")
  speed <- derivs$speed
  graphics::plot(x = derivs$speedTimes, y = speed, type = 'l', xlab = xlab, ylab = ylab)
  graphics::abline(h = attr(x, "slowerThan"), col = slowerThanColour)
  graphics::abline(h = attr(x, "fasterThan"), col = fasterThanColour)
  if (nrow(x) > 0)
    graphics::rect(x$startTime, min(speed), x$stopTime, max(speed), col = highlightColor, border = NA)
}
