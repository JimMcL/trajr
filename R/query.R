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

#' Trajectory number of frames
#'
#' Returns the number of frames recorded for this trajectory
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetNFrames <- function(trj) { attr(trj, .TRAJ_NFRAMES) }

#' Trajectory spatial units
#'
#' Returns the spatial units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajScale}}
#'
#' @export
TrajGetUnits <- function(trj) { attr(trj, .TRAJ_UNITS) }

# ---- Trajectory analysis ----

#' Trajectory step lengths
#'
#' Returns the lengths of each step in a trajectory.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajStepLengths <- function(trj) Mod(tail(trj$displacement, -1)) # First displacement is not a step (also usually 0)

#' Turning angles of a Trajectory
#'
#' Calculates the step angles (in radians) of each segment, either relative to
#' the previous segment or relative to the specified compass direction.
#'
#' @param trj the trajectory whose whose angles are to be calculated.
#' @param lag Angles between every lag'th segment is calculated.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @return Step angles in radians, normalised so that \code{-pi < angle <= pi}.
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

#' Calculates trajectory speed and acceleration
#'
#' Calculates speed and linear acceleration along a trajectory over time.
#'
#' @param trj Trajectory whose speed and acceleration is to be calculated.
#'
#' @return A list with components: \item{speed}{numeric vector, speed between
#'   each pair of trajectory points.} \item{speedTimes}{numeric vector, times
#'   corresponding to values in \code{speed}.} \item{acceleration}{numeric
#'   vector.} \item{accelerationTimes}{numeric vector.}
#'
#' @seealso \code{\link{TrajSpeedIntervals}} for analysing intervals within the
#'   trajectory of low or high speed.
#'
#' @export
TrajDerivatives <- function(trj) {
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
    # Special case if starting point was added as the start of an interval, leave it unchanged
    if (!(startFrames[i] == 1 && .frameIsInInterval(1)))
      startTimes[i] <- .interp(startFrames[i])
  }
  for (i in 1:length(stopFrames)) {
    # Special case if stopping point was added as the end of an interval, leave it unchanged
    lastStop <- length(stopFrames)
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
#' @param fasterThan If not \code{NULL}, intervals will cover time periods where
#'   speed exceeds this value.
#' @param slowerThan If not \code{NULL}, intervals will cover time periods where
#'   speed is lower than this value.
#' @param interpolateTimes If \code{TRUE}, times will be linearly interpolated
#'   between frames.
#'
#' @return data.frame, each row is an interval, with columns:
#'   \item{startFrame}{Indices of frames at the start of each interval.}
#'   \item{stopFrame}{Indices of frames at the end of each interval.}
#'   \item{startTime}{Time at the start of each interval.} \item{stopTime}{Time
#'   at the end of each interval} \item{duration}{Duraction of each interval.}.
#'
#' @seealso \code{\link{TrajDerivatives}} for calculating trajectory speed and
#'   acceleration.
#'
#' @examples
#' plotIntervals <- function(smoothed, slowerThan = NULL, fasterThan = NULL) {
#'   derivs <- TrajDerivatives(smoothed)
#'   speed <- derivs$speed
#'   plot(x = derivs$speedTimes, y = speed, type = 'l',
#'        xlab = 'Time (sec)', ylab = "Speed")
#'   abline(h = slowerThan, col = "red")
#'   abline(h = fasterThan, col = "green")
#'   intervals <- TrajSpeedIntervals(smoothed, fasterThan = fasterThan, slowerThan = slowerThan)
#'   rect(intervals$startTime, min(speed), intervals$stopTime, max(speed),
#'        col = "#0000FF1E", border = NA)
#'   intervals
#' }
#'
#' # Plot speed, highlighting intervals where speed drops below 50 units/sec
#' set.seed(4)
#' trj <- TrajGenerate(200, random = TRUE)
#' smoothed <- TrajSmoothSG(trj, 3, 101)
#' intervals <- plotIntervals(smoothed, 50, NULL)
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

  data.frame(startFrame = startFrames, startTime = startTimes, stopFrame = stopFrames, stopTime = stopTimes, duration = durations)
}


#' Path length
#'
#' Calculates the cumulative length of a track, which is the total distance
#' travelled along the trajectory.
#'
#' @param trj Trajectory whose length is to be calculated.
#' @return Path Numeric length of the trajectory.
#'
#' @export
TrajLength <- function(trj) {
  sum(Mod(diff(trj$polar)))
}

#' Calculates the distance between the start and end of a trajectory.
#' Also called the diffusion distance, net distance, or bee-line from start to finish.
#'
#' @param trj Trajectory whose distance is to be calculated.
#' @return Numeric distance from the start to the end of the trajectory.
#'
#' @export
TrajDistance <- function(trj) {
  Mod(diff(trj$polar[c(1,length(trj$polar))]))
}
