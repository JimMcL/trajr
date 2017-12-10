# Trajectory plotting and basic query functions


# ---- Trajectory plotting ----

#' Plot method for trajectories
#'
#' The plot method for Trajectory objects.
#'
#' @param x An object of class "Trajectory", the trajectory to be plotted.
#' @param draw.start.pt if TRUE, draws a dot at the start point of the
#'   trajectory.
#' @param add If TRUE, the trajectory is added to the current plot.
#' @param turning.angles If \code{random} or \code{directed}, draws step turning
#'   angles. \code{directed} assumes errors are relative to the first recorded
#'   step angle. \code{random} assumes errors are relative to the previous step.
#' @param type,xlim,ylim,xlab,ylab,asp plotting parameters with useful defaults.
#' @param ... Additional arguments are passed to \code{\link{plot}}.
#'
#' @seealso \code{\link{TrajFromCoords}}
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate(angularErrorSd = 1.3)
#' plot(trj)
#'
#' @export
plot.Trajectory <- function(x, draw.start.pt = TRUE, add = FALSE, turning.angles = NULL,
                            type = 'l',
                            xlim = extendrange(x$x), ylim = extendrange(x$y),
                            xlab = "x", ylab = "y",
                            asp = 1, ...) {
  if (!add) {
    graphics::plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, asp = asp, ...)
  }
  graphics::lines(y ~ x, data = x, type = type, ...)
  if (draw.start.pt)
    graphics::points(x$x[1], x$y[1], pch = 16, cex = .8)

  if (!is.null(turning.angles)) {
    # There are n steps, but n+1 coordinates
    n <- nrow(x) - 1
    steps <- x[1:n,]
    angles <- x[2:(n+1),]
    meanStepLength <- TrajMeanStepLength(x)
    textDisplacement <- 0.3 * meanStepLength
    labels <- parse(text= paste("Delta[", 1:n, "]", sep=""))

    if (tolower(turning.angles) == "directed") {
      # Plot angles which represent angular errors, which reset at each step
      segments(steps$x, steps$y, steps$x + .8 * meanStepLength, steps$y, col = "darkgrey", lty = 2)

      textAngle <- Arg(angles$displacement) +
        ifelse(Arg(x$displacement[1]) < Arg(angles$displacement), pi / 4, -pi / 6)
      text(steps$x + textDisplacement * cos(textAngle), steps$y + textDisplacement * sin(textAngle),
           labels = labels)

      plotrix::draw.arc(steps$x, steps$y,
                        angle1 = Arg(x$displacement[1]), angle2 = Arg(angles$displacement),
                        radius = 0.4 * meanStepLength)

    } else if (tolower(turning.angles) == "random") {
      # Plot angles which represent angular errors, which accumulate
      segments(steps$x, steps$y, steps$x + .8 * meanStepLength * cos(Arg(steps$displacement)),
               steps$y + 1.7 * sin(Arg(steps$displacement)), col = "darkgrey", lty = 2)

      textAngle <- Arg(angles$displacement) +
        ifelse(Arg(steps$displacement) < Arg(angles$displacement), pi / 4, -pi / 6)
      text(steps$x + textDisplacement * cos(textAngle), steps$y + textDisplacement * sin(textAngle),
           labels = labels)

      plotrix::draw.arc(steps$x, steps$y,
                        angle1 = Arg(steps$displacement), angle2 = Arg(angles$displacement),
                        radius = 0.4 * meanStepLength)

    } else {
      stop(sprintf("Invalid turning.angles (%s), must be one of 'random' or 'directed'", turning.angles))
    }
  }
}

# ---- Trajectory query ----

#' Trajectory frmes-per-second
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

# ---- Trajectory analysis ----

#' Mean trajectory step length
#'
#' Returns the mean segment length of a trajectory
TrajMeanStepLength <- function(trj) mean(Mod(trj$displacement))

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
    angles <- diff(Arg(trj$displacement), lag)
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

#' Calculate speed time intervals
#'
#' Calculates and returns a list of time intervals during which speed is slower
#' and/or faster than specified values.
#'
#' @param trj Trajectory to be analysed.
#' @param fasterThan If not NULL, intervals will cover time periods where speed
#'   exceeds this value.
#' @param slowerThan If not NULL, intervals will cover time periods where speed
#'   is lower than this value.
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
#' plotIntervals <- function(smoothed, intervals, slowerThan = NULL, fasterThan = NULL) {
#'   derivs <- TrajDerivatives(smoothed)
#'   speed <- derivs$speed
#'   plot(x = derivs$speedTimes, y = speed, type = 'l',
#'        xlab = 'Time (sec)', ylab = "Speed")
#'   abline(h = slowerThan, col = "red")
#'   abline(h = fasterThan, col = "green")
#'   rect(intervals$startTime, min(speed), intervals$stopTime, max(speed),
#'        col = "#0000FF1E", border = NA)
#' }
#'
#' # Plot speed, highlighting intervals where speed drops below 50 units/sec
#' set.seed(4)
#' trj <- TrajGenerate(200, random = TRUE)
#' smoothed <- TrajSmoothSG(trj, 3, 101)
#' intervals <- TrajSpeedIntervals(smoothed, slowerThan = 50)
#' plotIntervals(smoothed, intervals, 50, NULL)
#'
#' # Report the duration of the maximum period of low speed
#' cat(sprintf("Duration of the longest low-speed interval was %g secs\n", max(intervals$duration)))
#'
#' @export
TrajSpeedIntervals <- function(trj, fasterThan = NULL, slowerThan = NULL) {
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
  startTimes <- derivs$speedTimes[startFrames]
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
