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

#' Trajectory number of coordinates
#'
#' Returns the number of coordinates recorded for this trajectory, i.e. 1 more than the number of steps.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetNCoords <- function(trj) { attr(trj, .TRAJ_NFRAMES) }

#' Trajectory spatial units
#'
#' Returns the spatial units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajScale}}, \code{\link{TrajGetTimeUnits}}.
#'
#' @export
TrajGetUnits <- function(trj) { attr(trj, .TRAJ_UNITS) }

#' Trajectory temporal units
#'
#' Returns the temporal units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajFromCoords}}, \code{\link{TrajGetUnits}}.
#'
#' @export
TrajGetTimeUnits <- function(trj) {
  attr(trj, .TRAJ_TIME_UNITS)
}

# ---- Trajectory analysis ----

#' Trajectory step lengths
#'
#' Returns the lengths of each step in a trajectory.
#'
#' @param trj Trajectory to query.
#'
#' @seealso \code{\link{TrajLength}}
#'
#' @export
TrajStepLengths <- function(trj) {
  Mod(utils::tail(trj$displacement, -1)) # First displacement is not a step (also usually 0)
}

#' Trajectory distance
#'
#' Calculates the distance between the start and end of a trajectory.
#' Also called the diffusion distance, net distance, or bee-line from start to finish.
#'
#' @param trj Trajectory whose distance is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#' @return Numeric distance from the start to the end of the trajectory.
#'
#' @export
TrajDistance <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  Mod(diff(trj$polar[c(startIndex, endIndex)]))
}

#' Trajectory length
#'
#' Calculates the cumulative length of a track (or a portion of a track), which
#' is the total distance travelled along the trajectory.
#'
#' @param trj Trajectory whose length is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#' @return Path Numeric length of the trajectory.
#'
#' @seealso \code{\link{TrajStepLengths}}
#'
#' @export
TrajLength <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  sum(Mod(diff(trj$polar[startIndex:endIndex])))
}

#' Turning angles of a Trajectory
#'
#' Calculates the step angles (in radians) of each segment, either relative to
#' the previous segment or relative to the specified compass direction.
#'
#' @param trj the trajectory whose whose angles are to be calculated.
#' @param lag Angles between every lag'th segment are calculated.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @return Step angles in radians, normalised so that \code{-pi < angle <= pi}.
#'
#' @seealso \code{\link{TrajStepLengths}}
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
