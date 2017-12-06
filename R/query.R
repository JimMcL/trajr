# Trajectory plotting and basic query functions


# ---- Trajectory plotting ----

#' Plot method for trajectories
#'
#' The plot method for Trajectory objects.
#'
#' @param trj The trajectory to be plotted.
#' @param draw.start.pt if TRUE, draws a dot at the start point of the trajectory.
#' @param add If TRUE, the trajectory is added to the current plot.
#' @param type,xlim,ylim,asp plotting parameters with useful defaults.
#'
#' @export
plot.Trajectory <- function(trj, draw.start.pt = TRUE, add = FALSE,
                            type = 'l',
                            xlim = range(trj$x), ylim = range(trj$y),
                            asp = 1, ...) {
  if (!add) {
    NextMethod(NULL, xlim = xlim, ylim = ylim, asp = 1, ...)
  }
  lines(y ~ x, data = trj, type = type, ...)
  if (draw.start.pt)
    points(trj$x[1], trj$y[1], pch = 16, cex = .8)
}

# ---- Trajectory query ----

#' Returns the frames-per-second recorded for this trajectory
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetFPS <- function(trj) { attr(trj, .MTA_FPS) }

#' Returns the number of frames recorded for this trajectory
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetNFrames <- function(trj) { attr(trj, .MTA_NFRAMES) }

# ---- Trajectory analysis ----

#' Turning angles of a Trajectory
#'
#' Calculates the angles (in radians) of each segment relative to the previous segment.
#'
#' @param trj the trajectory whose whose angles are to be calculated.
#' @param lag Angles between every lag'th segment is calculated.
#'
#' @export
TrajAngles <- function(trj, lag = 1) {
  angles <- diff(Arg(trj$displacement), lag)
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
#' @return A list with components: \item{speed}{numeric vector, speed between each pair of trajectory points.}
#'   \item{speedTimes}{numeric vector, times corresponding to values in \code{speed}.}
#'   \item{acceleration}{numeric vector.}
#'   \item{accelerationTimes}{numeric vector.}
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
