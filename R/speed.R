# Functions relating to trajectory speed and acceleration

#### Velocity, acceleration & speed

# Internal function which just throws an error if trj doesn't contain time information
.checkTrajHasTime <- function(trj) {
  if (!("displacementTime" %in% names(trj)))
    stop("Missing time information in trajectory. Perhaps this is a rediscretized trajectory?")
}

# Returns sums of adjacent pairs of elements. This is similar to diff(v), except
# that elements are added rather than subtracted
.sumPairs <- function(v) {
  i <- utils::head(seq_along(v), -1)
  v[i] + v[i + 1]
}

#' Calculates trajectory speed and change of speed
#'
#' Calculates speed and change of speed along a trajectory over time. These are
#' the first and second order derivatives of distance travelled over time. Noisy
#' trajectories should be smoothed before being passed to this function, as
#' noise is effectively amplifed when taking derivatives.
#'
#' The value returned as \code{acceleration} is \emph{not} technically
#' acceleration. In mechanics, acceleration is a vector. This value is a scalar
#' quantity: change of speed, which is sometimes known informally as
#' acceleration. This value corresponds to the acceleration in a 1-dimensional
#' trajectory, with the sign indicating the direction of acceleration relative
#' to the current direction of velocity. See \code{\link{TrajAcceleration}} for
#' an approximation of (vector) acceleration, and \code{\link{TrajVelocity}} for
#' an approximation of velocity.
#'
#' @param trj Trajectory whose speed and change in speed is to be calculated.
#'
#' @return A list with components: \item{speed}{numeric vector, speed between
#'   each pair of trajectory points, i.e. the speed of each step.}
#'   \item{speedTimes}{numeric vector, times corresponding to values in
#'   \code{speed}, i.e. the time from the start of the trajectory to the end of
#'   each step.} \item{acceleration}{numeric vector, change in speed between
#'   steps. Despite the name, this is not acceleration as defined by mechanics.}
#'   \item{accelerationTimes}{numeric vector, time from the start of the
#'   trajectory to the end of the second step in each pair.}
#'
#' @seealso \code{\link{TrajSpeedIntervals}} for analysing intervals of low or
#'   high speed within the trajectory. \code{\link{TrajSmoothSG}} for smoothing
#'   a trajectory. \code{\link{TrajAcceleration}} for calculating acceleration,
#'   and \code{\link{TrajVelocity}} for calculating velocity.
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

#' Approximates the acceleration of a trajectory
#'
#' Returns an approximation of the acceleration of a trajectory at each point
#' using the second-order central finite differences.
#'
#' `trajr` trajectories, which consist of straight line displacements between
#' sampled locations, do not contain enough information to correctly derive
#' velocity or acceleration. Since we have to assume a constant velocity at each
#' step, the first derivative is discontinuous. Acceleration, therefore, is zero
#' during each step and infinite at each change of velocity. The approximation
#' implemented by this function assumes that acceleration occurs over a period
#' of time: half the duration of the previous step plus half the duration of the
#' next step.
#'
#' @param trj Trajectory whose acceleration is to be calculated.
#'
#' @return Vector of complex numbers. The modulus (\code{Mod(a)}) is the
#'   magnitude of the acceleration at each point, and the argument
#'   (\code{Arg(a)}) is the direction of the acceleration. The vector has an
#'   attribute, \code{trj}, with the trajectory as its value. The first and last
#'   values will always be \code{NA}, since acceleration cannot be estimated for
#'   those points.
#'
#' @seealso \code{\link{TrajVelocity}} for calculating velocity,
#'   \code{\link{TrajResampleTime}} and \code{\link{TrajRediscretize}} to
#'   resample a trajectory to fixed time or length steps.
#'
#' @examples
#' # A function to plot acceleration as arrows (scaled in length)
#' AccArrows <- function(acc, scale = .001, trj = attr(acc, "trj"), ...) {
#'   graphics::arrows(trj$x, trj$y, trj$x + Re(acc) * scale, trj$y + Im(acc) * scale, ...)
#' }
#'
#' # Generate and plot a random trajectory
#' set.seed(101)
#' trj <- TrajGenerate(30)
#' plot(trj)
#'
#' # Calculate acceleration
#' acc <- TrajAcceleration(trj)
#' # Plot acceleration as red arrows at each point. They need to be scaled down to
#' # fit in the plot, and the arrowhead lengths need to be shortened to look good
#' AccArrows(acc, scale = .001, col = "red", length = .1)
#'
#' @export
TrajAcceleration <- function(trj) {
  .checkTrajHasTime(trj)

  # Note that there's no point in calculating backward or forward differences as
  # they provide no different information

  # If we were guaranteed a constant step time, h, we could use the more elegant
  # ax <- stats::filter(trj$x, c(1, -2, 1)) / h^2

  x <- trj$x
  y <- trj$y
  h <- diff(trj$time)

  # Calculate velocities using forward/backward diffs
  vx <- diff(x) / h
  # calculate acceleration from velocity and time
  ax <- diff(vx) / .sumPairs(h / 2)
  vy <- diff(y) / h
  ay <- diff(vy) / .sumPairs(h / 2)

  acc <- complex(real = c(NA, ax, NA), imaginary = c(NA, ay, NA))
  attr(acc, "trj") <- trj
  acc
}

#' Velocity of a trajectory
#'
#' The velocity is approximated at each point of the trajectory using
#' first-order finite differences. Central, forward or backward differences can
#' be used. Central differences yield a more accurate approximation if the
#' velocity is smooth. As a practical guide, if velocity doesn't change much
#' between steps, use central differences. If it changes substantially (and not
#' just as an artifact of recording noise), then use either forward or backward
#' differences.
#'
#' Intuitively, think of the central difference velocity at a point as the mean
#' of the velocities of the two adjacent steps. Forward difference velocity is
#' the velocity of the step starting at the point. Backward difference is the
#' velocity of the step ending at the point.
#'
#' @param trj Trajectory whose velocity is to be calculated.
#' @param diff Type of difference to be calculated, one of "central" (the
#'   default), "forward" or "backward".
#'
#' @return A vector of complex numbers representing the velocity at each point
#'   along the trajectory. The modulus (\code{Mod(v)}) is the magnitude of the
#'   velocity, i.e. the speed; the argument (\code{Arg(v)}) is the direction of
#'   the velocity; the real part (\code{Re(v)}) is velocity in the X direction;
#'   and the imaginary part (\code{Im(v)}) is velocity in the Y direction. The
#'   vector has an attribute, \code{trj}, with the trajectory as its value. If
#'   \code{diff} is \code{"central"}, the first and last velocity values are
#'   \code{NA} since velocity cannot be calculated for them. If \code{diff} is
#'   \code{"forward"}, the last value will be NA, and if \code{diff} is
#'   \code{"backward"}, the first value will be NA.
#'
#' @seealso \code{\link{TrajAcceleration}} for calculating acceleration;
#' \code{\link{TrajResampleTime}} and \code{\link{TrajRediscretize}} to resample
#' a trajectory to fixed time or length steps; \code{\link{TrajSpeedIntervals}}
#' for calculating when speed crosses some threshold; Finite differences on
#' \href{https://en.wikipedia.org/wiki/Finite_difference}{Wikipedia}.
#'
#' @examples
#' set.seed(11)
#' trj <- TrajGenerate(100)
#' # calculate velocity
#' vel <- TrajVelocity(trj)
#'
#' # Obtain speed over time, with NAs removed
#' speed <- na.omit(data.frame(speed = Mod(vel), time = trj$time))
#'
#' plot(speed ~ time, speed, type = 'l')
#'
#' @export
TrajVelocity <- function(trj, diff = c("central", "forward", "backward")) {
  .checkTrajHasTime(trj)
  diff <- match.arg(diff)

  x <- trj$x
  y <- trj$y
  h <- diff(trj$time)

  if (diff == "central") {

    # Central diffs (this is the "double-interval" central difference)
    # Sum time (h) for each adjacent pair of steps
    dt <- c(NA, .sumPairs(h), NA)
    vx <- stats::filter(x, c(1, 0, -1)) / dt
    vy <- stats::filter(y, c(1, 0, -1)) / dt
  } else {

    # Forward or backward diffs, variable step times
    vx <- diff(x) / h
    vy <- diff(y) / h
    if (diff == "forward") {
      # Forward diffs - speed at last point is unknown
      vx <- c(vx, NA)
      vy <- c(vy, NA)
    } else {
      # Backward diffs - speed at first point is unknown
      vx <- c(NA, vx)
      vy <- c(NA, vy)
    }
  }

  vel <- complex(real = vx, imaginary = vy)
  attr(vel, "trj") <- trj
  vel
}


##########################################################################
#### Speed intervals

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
#' and/or faster than specified values. Speed is calculated by taking the
#' modulus of velocity (\code{\link{TrajVelocity}}).
#'
#' @param trj Trajectory to be analysed.
#' @param fasterThan,slowerThan If not \code{NULL}, intervals will cover time
#'   periods where speed exceeds/is lower than this value.
#' @param interpolateTimes If \code{TRUE}, times will be linearly interpolated
#'   between frames.
#' @param diff Method used to calculate speed, see \code{\link{TrajVelocity}}
#'   for details. The default is \code{"backward"} to maintain backwards
#'   compatibility; in general, \code{"central"} provides a more accurate
#'   estimate of velocity.
#'
#' @return A data frame of class "TrajSpeedIntervals", each row is an interval,
#'   columns are: \item{startFrame}{Indices of frames at the start of each
#'   interval.} \item{stopFrame}{Indices of frames at the end of each interval.}
#'   \item{startTime}{Time since start of trajectory at the start of each
#'   interval.} \item{stopTime}{Time since start of trajectory at the end of
#'   each interval} \item{duration}{Duration of each interval.}
#'
#'   The data frame will also have non-standard attributes:
#'   \item{trajectory}{Value of the \code{trj} argument.}
#'   \item{slowerThan}{Value of the \code{slowerThan} argument.}
#'   \item{fasterThan}{Value of the \code{fasterThan} argument.}
#'   \item{speed}{Data frame with columns \code{speed} and \code{time}.}
#'   \item{derivs}{Value returned by calling \code{TrajDerivatives(trj)}.
#'   Provided for backwards-compatibility; use of \code{speed} is now preferred.}
#'
#' @seealso \code{\link{TrajVelocity}} for calculating trajectory velocity,
#'   \code{\link{plot.TrajSpeedIntervals}} for plotting speed over time with
#'   intervals highlighted.
#'
#' @examples
#' # Plot speed, highlighting intervals where speed drops below 50 units/sec
#' set.seed(4)
#' trj <- TrajGenerate(200, random = TRUE)
#' smoothed <- TrajSmoothSG(trj, 3, 101)
#' intervals <- TrajSpeedIntervals(smoothed, diff = "central", slowerThan = 50, fasterThan = NULL)
#' plot(intervals)
#'
#' # Report the duration of the longest period of low speed
#' cat(sprintf("Duration of the longest low-speed interval was %g secs\n", max(intervals$duration)))
#'
#' @export
TrajSpeedIntervals <- function(trj, fasterThan = NULL, slowerThan = NULL, interpolateTimes = TRUE, diff = c("backward", "central", "forward")) {
  if (is.null(fasterThan) && is.null(slowerThan)) {
    stop("Parameters fasterThan and slowerThan are both NULL, one must be specified")
  }
  diff <- match.arg(diff)

  # Calculate trajectory speeds
  vel <- TrajVelocity(trj, diff = diff)
  df <- stats::na.omit(data.frame(speed = Mod(vel), time = trj$time - trj$time[1]))
  speed <- df$speed
  times <- df$time

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

  # Handle special cases: situation where interval begins or ends outside of trajectory
  if (length(startFrames) > 0 || length(stopFrames) > 0) {
    # Assume interval started at beginning of trajectory, since we don't know what happened before that
    if (length(stopFrames) > 0 && (length(startFrames) == 0 || stopFrames[1] < startFrames[1]))
      startFrames <- c(1, startFrames)
    # Similarly, assume that interval can't extend past end of trajectory
    if (length(stopFrames) == 0 || startFrames[length(startFrames)] > stopFrames[length(stopFrames)])
      stopFrames <- c(stopFrames, length(speed))

  } else if (all(flags) && length(flags) > 0) {
    # The entire trajectory is an interval
    startFrames <- 1
    stopFrames <- length(speed)
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
  attr(result, "speed") <- df
  attr(result, "derivs") <- TrajDerivatives(trj) # This is only here for backwards compatibility
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
#' @param xlab,ylab,type Plotting parameters with useful defaults.
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
                                    type = "l",
                                    ...) {
  trj <- attr(x, "trajectory")
  df <- attr(x, "speed")
  graphics::plot(x = df$time + trj$time[1], y = df$speed, type = type, xlab = xlab, ylab = ylab, ...)
  graphics::abline(h = attr(x, "slowerThan"), col = slowerThanColour)
  graphics::abline(h = attr(x, "fasterThan"), col = fasterThanColour)
  if (nrow(x) > 0) {
    plotExtents <- graphics::par("usr")
    graphics::rect(x$startTime + trj$time[1], plotExtents[3], x$stopTime + trj$time[1], plotExtents[4], col = highlightColor, border = NA)
  }
}
