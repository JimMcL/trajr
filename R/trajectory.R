# Sesiidae trajectory analysis


.MTA_FPS <- 'fps'
.MTA_NFRAMES <- 'numFrames'

#' Create a Trajectory Object
#'
#' \code{TrajFromCoords} creates a new trajectory object from a set of
#' 2-dimensional cartesian coordinates and some metadata.
#'
#' If \code{track} contains a \code{time} (or \code{Time}) column, it is assumed
#' to contain the time (in seconds) of each data point. Otherwise times are
#' calculated for each point as \code{(frame - 1) / fps} where \code{frame} is
#' the index of the point.
#'
#' @param track data frame containing x & y columns.
#' @param xCol Name or index of the \code{x} column in \code{track} (default 1).
#' @param yCol Name or index of the \code{y} column in \code{track} (default 2).
#' @param fps Frames per second - used to calculate relative frame times if
#'   \code{track} does not contain a \code{time} column.
#'
#' @return A new trajectory object.
TrajFromCoords <- function(track, xCol = 1, yCol = 2, fps = 50) {
  # Get polar coordinates
  track$polar <- complex(real = track$x, imaginary = track$y)

  # Calculate displacements from each point to the next
  track$displacement <- c(0, diff(track$polar))

  # Rename time column to standard name
  if (!('time' %in% names(track)) && 'Time' %in% names(track)) {
      names(track)[names(track) == 'Time'] <- "time"
  }

  # Allocate times if they aren't already known
  if (!('time' %in% names(track))) {
    # Assign times to each frame, starting at 0
    track$time <- 0:(nrow(track) - 1) / fps
  }

  # Get times associated with displacements, with the first point at time 0,
  # i.e. time at each point in displacement, not time between points
  track$displacementTime <- track$time[1:nrow(track)] - track$time[1]

  # Save number of frames
  attr(track, .MTA_NFRAMES) <- nrow(track)

  # Save frame rate
  attr(track, .MTA_FPS) <- fps

  track
}

TrajGetFPS <- function(track) { attr(track, .MTA_FPS) }
TrajGetNFrames <- function(track) { attr(track, .MTA_NFRAMES) }

TrajPlot <- function(track, type = 'l', add = FALSE, ...) {
  if (!add) {
    plot(NULL, xlim = range(track$x), ylim = range(track$y), asp = 1, ...)
  }
  lines(y ~ x, data = track, type = type, ...)
  points(track$x[1], track$y[1], pch = 16, cex = .8)
}

# Returns angles (in radians) of each segment relative to the previous segment
TrajAngles <- function(track, lag = 1) {
  angles <- diff(Arg(track$displacement), lag)
  # Normalise so that -pi < angle <= pi
  ii <- angles <= -pi
  angles[ii] <- angles[ii] + 2 * pi
  ii <- angles > pi
  angles[ii] <- angles[ii] - 2 * pi
  angles
}

# Calculates speed and linear acceleration along a track over time.
#
# value - list with values speed, speedTimes, acceleration, accelerationTimes
TrajDerivatives <- function(track) {
  # Note that displacements are the (polar) displacements from 1 point to the next
  d <- Mod(track$displacement)
  t <- track$displacementTime

  # Calculate speed
  v <- d[2:length(d)] / diff(t)
  vt <- t[2:length(t)]
  # Calculate linear acceleration
  a <- diff(v) / diff(vt)
  at <- vt[2:length(vt)]

  list(speed = v, speedTimes = vt, acceleration = a, accelerationTimes = at)
}

# Calculates the cumulative length of a track, which is the distance travelled
TrajLength <- function(track) {
  sum(Mod(diff(track$polar)))
}

# Calculates the distance between the start and end of a set of points in polar form.
# Also called the diffusion distance, net distance, or bee-line from start to finish.
TrajDistance <- function(track) {
  Mod(diff(track$polar[c(1,length(track$polar))]))
}

TrajSinuosity <- function(track) {
  TrajDistance(track) / TrajLength(track)
}

# Rotates a track so that angle(finish - start) == angle
TrajRotate <- function(track, angle = 0) {
  # Calculate current orientation
  orient <- Arg(track$polar[length(track$polar)] - track$polar[1])
  # Calculate required rotation
  alpha <- angle - orient
  # Rotation matrix
  rm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)

  # New track is old track rotated
  nt <- as.data.frame(t(rm %*% (t(track[,c('x', 'y')]))))
  colnames(nt) <- c('x', 'y')
  nt <- TrajFromCoords(nt)

  nt
}

# Directional change (DC) for an entire track
# track - Track to calculate DC for
# nFrames - frame delta to process - if 1, every frame is process, if 2, every 2nd frame is processed...
#
# Kitamura, T., & Imafuku, M. (2015). Behavioural mimicry in flight path of Batesian intraspecific polymorphic butterfly Papilio polytes. Proceedings of the Royal Society B: Biological Sciences, 282(1809). doi:10.1098/rspb.2015.0483
TrajDirectionalChange <- function(track, nFrames = 1) {
  TrajAngles(track, nFrames) / diff(track$displacementTime, nFrames)
}

# Smooths a track using a Savitzky-Golay smoothing filter
TrajSmooth <- function(track, p = 3, n = p + 3 - p%%2) {
  sm <- data.frame(x = sgolayfilt(track$x, p, n), y = sgolayfilt(track$y, p, n))
  TrajFromCoords(sm, TrajGetFPS(track))
}



