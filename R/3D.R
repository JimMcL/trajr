# 3-dimensional trajectories


# ---- Trajectory creation and modification ----

#' Create a 3D Trajectory Object
#'
#' \code{Traj3DFromCoords} creates a new 3-dimensional trajectory object from a
#' set of 3-dimensional cartesian coordinates, times and some metadata. A 3D
#' trajectory is a 2D trajectory (see \code{\link{TrajFromCoords}}) with the
#' addition of a \code{z} coordinate. \code{trajr} functions that expect a 2D
#' trajectory will work on a 3D trajectory by simply ignoring the \code{z}
#' dimension, so should only be used with care. A small number of functions
#' operate on 3D trajectories, and are prefixed by \code{Traj3D}. Existing
#' \code{trajr} functions that that do not deal with spatial data (e.g.
#' \code{\link{TrajDuration}}, \code{\link{TrajGetNCoords}} etc.) can safely be
#' used.
#'
#' The coordinates are sometimes referred to as "relocations". Rows containing
#' \code{NA} coordinate or time values at the start or end of the trajectory are
#' discarded. \code{NA} coordinates or times in the middle of the trajectory
#' generate an error.
#'
#' If \code{timeCol} is specified, \code{track[,timeCol]} is expected to contain
#' the time (in some numeric units) of each coordinate. Otherwise, times are
#' calculated for each point as \code{(coord - 1) / fps} where \code{coord} is
#' the index of the point; in other words, sampling at constant time intervals
#' is assumed. Time values require conversion if they are not numeric. It may be
#' possible to use `strptime` for this purpose, or \code{\link{TrajConvertTime}}
#' can be used to convert multiple field time values.
#'
#' \code{x} and \code{y} must be square units. Longitude and latitude are not
#' suitable for use as \code{x} and \code{y} values, since in general, \code{1°
#' lat != 1° lon}. To create a trajectory from positions in latitude and
#' longitude, it is first necessary to transform the positions to a suitable
#' spatial projection such as UTM (possibly by using \code{spTransform} from the
#' \code{rgdal} package).
#'
#' Leading and trailing rows with \code{NA} coordinate values are discarded.
#' \code{NA} coordinate values within a trajectory generate an error.
#'
#' @param track data frame containing cartesian coordinates and optionally times
#'   for the points in the trajectory.
#' @param xCol Name or index of the \code{x} column in \code{track} (default 1).
#' @param yCol Name or index of the \code{y} column in \code{track} (default 2).
#' @param zCol Name or index of the \code{z} column in \code{track} (default 3).
#' @param timeCol optional name or index of the column which contains coordinate
#'   times.
#' @param fps Frames per second - used to calculate relative coordinate times if
#'   \code{track} does not contain a \code{time} column. Time intervals between
#'   coordinate are assumed to be constant throught the entire track.
#' @param spatialUnits Abbreviation for the x, y and z units.
#' @param timeUnits Abbreviation for the units that time is recorded in.
#'
#' @return An object with class "\code{Trajectory3D}", which is a data.frame
#'   with the following components: \item{x}{X coordinates of trajectory
#'   points.} \item{y}{Y coordinates of trajectory points.} \item{time}{Time (in
#'   \code{timeUnits}) for each point. if \code{timeCol} is specified, values
#'   are \code{track[,timeCol]}, otherwise values are calculated from
#'   \code{fps}.} \item{displacementTime}{Relative frame/observation times, with
#'   frame/observation 1 at time \code{0}.} \item{polar}{X and y coordinates
#'   represented as complex numbers, to simplify working with 2D segment angles.
#'   Note that the z dimension is not represented.}
#'   \item{displacement}{2-dimensional displacement vectors (represented as
#'   complex numbers) between each pair of consecutive points. Note that the z
#'   dimension is not represented.}
#'
#' @seealso \code{\link{Traj3DLength}}, \code{\link{Traj3DStepLengths}},
#'   \code{\link{Traj3DDistance}}, \code{\link{Traj3DStraightness}},
#'   \code{\link{Traj3DSmoothSG}}, \code{\link{Traj3DResampleTime}},
#'   \code{\link{Traj3DRediscretize}}, \code{\link{TrajFromCoords}} for creating
#'   2-dimensional trajectories
#'
#' @export
Traj3DFromCoords <- function(track, xCol = 1, yCol = 2, zCol = 3,
                             timeCol = NULL, fps = 50,
                             spatialUnits = "m", timeUnits = "s") {

  trj3d <- track

  # Check Z coordinates are valid and remove leading/trailing NAs
  trj3d <- .checkCoords(trj3d, cols = zCol)

  # Create 2D trajectory
  trj3d <- TrajFromCoords(trj3d, xCol = xCol, yCol = yCol, timeCol = timeCol, fps = fps, spatialUnits = spatialUnits, timeUnits = timeUnits)

  # Rename Z column
  trj3d <- .renameCol(zCol, "z", trj3d)

  # Give it a special class
  class(trj3d) <- c(.TRAJ_3D_CLASS, class(trj3d))

  trj3d
}

#' 3D trajectory length
#'
#' Calculates the cumulative length of a 3-dimensional trajectory, or a portion
#' of a trajectory. The length is the total distance travelled along the trajectory.
#'
#' @param trj3d 3D trajectory whose length is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric length of the trajectory.
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link{Traj3DStepLengths}}, \code{\link{TrajLength}}
#'
#' @export
Traj3DLength <- function(trj3d, startIndex = 1, endIndex = nrow(trj3d)) {
  # Restrict to the rows and columns we are interested in
  pts <- trj3d[seq(startIndex, endIndex, 1), c("x", "y", "z")]
  # Calculate Euclidean distance between each pair of consecutive points, then add them all up
  sum(sapply(seq_len(nrow(pts) - 1), function(i) stats::dist(pts[i:(i + 1), ])))
}

#' 3D trajectory step lengths
#'
#' Returns the lengths of each step in a 3-dimensional trajectory (or part of a trajectory).
#'
#' @param trj3d Trajectory to query.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Vector of step lengths. The vector will have length \code{1 - TrajGetNCoords(trj3d)}.
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link{Traj3DLength}}, \code{\link{TrajStepLengths}}
#'
#' @export
Traj3DStepLengths <- function(trj3d, startIndex = 1, endIndex = nrow(trj3d)) {
  # Restrict to the rows and columns we are interested in
  pts <- trj3d[seq(startIndex, endIndex, 1), c("x", "y", "z")]
  # Calculate Euclidean distance between each pair of consecutive points, then add them all up
  sapply(seq_len(nrow(pts) - 1), function(i) stats::dist(pts[i:(i + 1), ]))
}

#' 3D Trajectory distance
#'
#' Calculates the distance between the start and end of a 3-dimensional trajectory (or a
#' portion of a trajectory). Also called the diffusion distance, net distance,
#' displacement, or bee-line from start to finish.
#'
#' @param trj3d 3-dimensional trajectory whose distance is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric distance from the start to the end of the trajectory.
#'
#' @seealso \code{\link{Traj3DFromCoords}}
#'
#' @export
Traj3DDistance <- function(trj3d, startIndex = 1, endIndex = nrow(trj3d)) {
  # Restrict to the rows and columns we are interested in
  pts <- trj3d[c(startIndex, endIndex), c("x", "y", "z")]
  as.numeric(stats::dist(pts))
}

#' Straightness of a 3D Trajectory
#'
#' Calculates the straightness index of a 3-dimensional trajectory, \eqn{D / L},
#' where \code{D} is the beeline distance between the first and last points in
#' the trajectory,and \code{L} is the path length travelled (Batschelet, 1981).
#' Benhamou (2004) considers the straightness index to be a reliable measure of
#' the efficiency of a directed walk, but inapplicable to random trajectories.
#' The straightness index of a random walk tends towards zero as the number of
#' steps increases, hence should only be used to compare the tortuosity of
#' random walks consisting of a similar number of steps.
#'
#' The straightness index is also known as the net-to-gross displacement ratio.
#' According to Batschelet (1981), this value (termed \emph{d}) is an
#' approximation of \emph{r}, which is the length of the mean vector of turning
#' angles of a constant step-length trajectory (see
#' \code{\link{TrajMeanVectorOfTurningAngles}} and
#' \code{\link{TrajRediscretize}} for creating a constant step-length
#' trajectory).
#'
#' @param trj3d 3-dimensional trajectory to calculate straightness of.
#' @return The straightness index of \code{trj}, which is a value between 0
#'   (infinitely tortuous) to 1 (a straight line).
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link{Traj3DDistance}} for
#'   trajectory distance (or displacement), and \code{\link{Traj3DLength}} for
#'   trajectory path length, \code{\link{Traj3DStraightness}} for the
#'   straightness of a 2D trajectory.
#'
#' @references
#'
#' Batschelet, E. (1981). Circular statistics in biology. ACADEMIC PRESS, 111
#' FIFTH AVE., NEW YORK, NY 10003, 1981, 388.
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's
#' path. Journal of Theoretical Biology, 229(2), 209-220.
#' doi:10.1016/j.jtbi.2004.03.016
#'
#' @export
Traj3DStraightness <- function(trj3d) {
  Traj3DDistance(trj3d) / Traj3DLength(trj3d)
}

#' Smooth a 3D trajectory using a Savitzky-Golay filter
#'
#' Smooths a 3-dimensional trajectory using a Savitzky-Golay smoothing filter.
#'
#' Consider carefully the effects of smoothing a trajectory with temporal gaps
#' in the data. If the smoothed trajectory is used to derive speed and/or
#' acceleration, it may be advisable to fill in the gaps before smoothing,
#' possibly by calling \code{Traj3DResampleTime}.
#'
#' @param trj3d The 3=dimensional trajectory to be smoothed.
#' @param p polynomial order (passed to \code{\link[signal]{sgolayfilt}}).
#' @param n Filter length (or window size), must be an odd number.  Passed to
#'   \code{\link[signal]{sgolayfilt}}.
#' @param ... Additional arguments are passed to
#'   \code{\link[signal]{sgolayfilt}}.
#' @return A new trajectory which is a smoothed version of the input trajectory.
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link[signal]{sgolayfilt}}, \code{\link{TrajSmoothSG}}
#'
#' @export
Traj3DSmoothSG <- function(trj3d, p = 3, n = p + 3 - p %% 2, ...) {
  if (n %% 2 != 1)
    stop(sprintf("Invalid smoothing parameter n (%d): n must be odd", n))
  if (n > nrow(trj3d))
    stop(sprintf("Invalid smoothing parameter n (%d): n must be less than the number of points in the trajectory (%d)", n, nrow(trj3d)))
  trj3d$x <- signal::sgolayfilt(trj3d$x, p, n, ...)
  trj3d$y <- signal::sgolayfilt(trj3d$y, p, n, ...)
  trj3d$z <- signal::sgolayfilt(trj3d$z, p, n, ...)
  .fillInTraj(trj3d, "3D")
}

#' Resample a 3D trajectory to a constant time interval
#'
#' Constructs a new 3-dimensional trajectory by resampling the input trajectory
#' to a fixed time interval. Points are linearly interpolated along the
#' trajectory. Spatial and time units are preserved.
#'
#' @param trj3d The 3-dimensional trajectory to be resampled.
#' @param stepTime The resampled trajectory step time. Each step in the new
#'   trajectory will have this duration.
#' @param newFps Value to be stored as the FPS value in the new trajectory (see
#'   \code{\link{TrajGetFPS}}). It is not otherwise used by this function.
#'
#' @return A new 3-dimensional trajectory with a constant time interval for each
#'   step. Points in the new trajectory are calculated by linearly interpolating
#'   along \code{trj3d}.
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link{TrajResampleTime}}
#'
#' @export
Traj3DResampleTime <- function(trj3d, stepTime, newFps = NULL) {
  # Determine times of new points
  times <- seq(from = min(trj3d$time), to = max(trj3d$time), by = stepTime)

  # Interpolate x, y and z separately. This works by treating x (y or z) as a
  # function of time (so x axis is time, y axis is trajectory x or y), then
  # interpolating to required times. approx interpolates y for given values of x
  x <- stats::approx(trj3d$time, trj3d$x, times)$y
  y <- stats::approx(trj3d$time, trj3d$y, times)$y
  z <- stats::approx(trj3d$time, trj3d$z, times)$y

  # Create the new trajectory, preserving metadata
  Traj3DFromCoords(data.frame(x, y, z, times), timeCol = 4,
                   fps = newFps,
                   spatialUnits = TrajGetUnits(trj3d),
                   timeUnits = TrajGetTimeUnits(trj3d))
}

##############################################################################

# This function is the guts of Traj3DRediscretize
#
# @param points list of points, with x, y & z values
# @param R rediscretization step length
#
# @return data frame with x, y, & z columns, which are the rediscretized path.
.Traj3DRediscretizePoints <- function(points, R) {

  # Define dot product for vectors
  dp <- function(x, y) sum(x * y)

  # Throw out other columns
  p <- points[, c("x", "y", "z")]

  # result will contain the points in discretized path points (X, Y, Z)
  np <- nrow(points)
  result <- data.frame(x = NA, y = NA, z = NA)
  # Start with the first point
  p0 <- result[1, ] <- p[1, ]

  # We will traverse each step, looking for the next point ion the rediscretized path
  step <- 1

  # For each point...
  while (step < np) {

    # Find the first point k for which |p[k] - p0| >= R
    k <- NA
    for (j in (step + 1):np) {
      # Calculate distance from p[j] to p0
      d <- stats::dist(rbind(p0, p[j,]))
      # cat(sprintf("step = %d, j = %d, d = %g\n", step, j, d))
      if (d >= R) {
        k <- j
        break;
      }
    }

    # cat(sprintf("Got k = %d\n", k))
    if (is.na(k)) {
      # We have reached the end of the path
      break
    }

    # The next point may lie on the same segment
    step <- k - 1

    # The point lies on the segment p[k-1], p[k], and also on the sphere
    # with origin p0 and radius R. Find the intersection of the line and the
    # sphere (https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection).
    # Use variables names from the Wikipedia article
    u <- p[k,] - p[k - 1,]    # Direction of line
    o <- p[k - 1, ]           # Origin of line
    c <- p0                   # Centre of sphere

    # Get roots of quadratic (-b +- sqrt(b^2 - 4ac)) / 2a. We always want the
    # intersection in the positive direction (i.e. going from p[k-1] to p[k]).
    aq <- dp(u, u)
    bq <- 2 * dp(u, o - c)
    cq <- dp(o - c, o - c) - R^2

    d <- (-bq + sqrt(bq^2 - 4 * aq * cq)) / (2 * aq)

    p0 <- o + d * u
    # print(p0)
    # Save the point
    result <- rbind(result, p0)
    # print(dist(result[(nrow(result) - 1): nrow(result), ]))
  }

  result
}

#' Resample a 3D trajectory to a constant step length
#'
#' Constructs a new 3-dimensional trajectory by resampling the input
#' 3-dimensional trajectory to a fixed step (or segment) length. By
#' default, timing of frames is lost, so speed and acceleration cannot
#' be calculated on a rediscretized trajectory. However, a constant
#' speed may be applied to the rediscretized trajectory
#' (\code{simConstantSpeed = TRUE}), in which case the returned
#' trajectory will have (almost) constant speed, with average speed
#' approximately equal to the average speed of \code{trj3d}.
#'
#' Unfortunately this operation is slow for large trajectories.
#'
#' Based on the appendix in Bovet and Benhamou, (1988), extended to 3
#' dimensions.
#'
#' @param trj3d The 3-dimensional trajectory to be resampled.
#' @param R rediscretization step length, in the spatial units of \code{trj}.
#' @param simConstantSpeed If TRUE, speeds are interpolated along the new
#'   trajectory so that average speed is approximately the same as that of
#'   \code{trj3d}.
#'
#' @return A new 3-dimensional trajectory with a constant segment length which
#'   follows \code{trj3d}.
#'
#' @seealso \code{\link{Traj3DFromCoords}}, \code{\link{TrajRediscretize}}
#'
#' @references Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals'
#'   movements using a correlated random walk model. Journal of Theoretical
#'   Biology, 131(4), 419-433. doi:10.1016/S0022-5193(88)80038-9
#'
#' @export
Traj3DRediscretize <- function(trj3d, R, simConstantSpeed = FALSE) {
  rt <- .Traj3DRediscretizePoints(trj3d, R)

  # Sanity check
  if (length(rt) < 2) {
    stop(sprintf("Step length %g is too large for path (path length %g)", R, TrajLength(trj3d)))
  }

  # Fill in other track stuff
  rownames(rt) <- NULL
  #.fillInTraj(rt, "3D")
  rt <- Traj3DFromCoords(rt, zCol = "z", spatialUnits = TrajGetUnits(trj3d))

  # Remove time columns
  rt$time <- NULL
  rt$displacementTime <- NULL

  # Optionally simulate a fixed speed trajectory with the same average speed as the original trajectory
  if (simConstantSpeed) {
    if (!"time" %in% names(trj3d))
      stop("Unable to simulate constant speed: missing time column in original trajectory")

    # Fill in displacementTime and time so that average speed is close to that of trj3d
    avgSpeed <- Traj3DLength(trj3d) / TrajDuration(trj3d)
    newDuration <- Traj3DLength(rt) / avgSpeed

    rt$displacementTime <- seq(0, newDuration, length.out = nrow(rt))
    rt$time <- rt$displacementTime + trj3d$time[1]

    # Infer a frame rate
    attr(rt, .TRAJ_FPS) <- mean(diff(rt$time))

    # Copy original time units
    attr(rt, .TRAJ_TIME_UNITS) <- TrajGetTimeUnits(trj3d)
  }

  rt
}

############################################################################################

#' Approximates the acceleration of a 3-dimensional trajectory
#'
#' Returns an approximation of the acceleration of a trajectory at each point
#' using the second-order central \href{https://en.wikipedia.org/wiki/Finite_difference}{finite
#' differences}.
#'
#' \code{trajr} trajectories, which consist of straight line displacements
#' between sampled locations, do not contain enough information to correctly
#' derive velocity or acceleration. Since we have to assume a constant velocity
#' at each step, the first derivative is discontinuous. Acceleration, therefore,
#' is zero during each step and infinite at each change of velocity. The
#' approximation implemented by this function assumes that acceleration occurs
#' over a period of time: half the duration of the previous step plus half the
#' duration of the next step.
#'
#' The function \code{\link{Traj3DSpeed}}, despite its name, can be used to
#' calculate the magnitude of the acceleration vectors returned by this function.
#'
#' @param trj3d 3-dimensional trajectory whose acceleration is to be calculated.
#'
#' @return Numeric matrix of 3D acceleration vectors. Each row represents the
#'   acceleration at a point in the trajectory. Columns are named "`x`", "`y`"
#'   and "`z`". The vector has an attribute, \code{trj}, with the trajectory as
#'   its value. The first and last values will always be \code{NA}, since
#'   acceleration cannot be estimated for those points.
#'
#' @seealso \code{\link{Traj3DVelocity}} for calculating velocity,
#'   \code{\link{Traj3DResampleTime}} and \code{\link{Traj3DRediscretize}} to
#'   resample a trajectory to fixed time or length steps.
#'
#' @examples
#' \dontrun{
#' library(rgl)
#'
#' # Plot a trajectory and its acceleration in 3D, using the rgl package
#'
#' # Function to add acceleration vectors as arrows to a 3D trajectory plot
#' Acc3DArrows <- function(acc, scale = 0.0001, trj3d = attr(acc, "trj3d"), ...) {
#'     cols <- c("x", "y", "z")
#'     sapply(seq_len(nrow(t3) - 2) + 1, function(r) {
#'            arrow3d(t3[r, cols], t3[r, cols] + acc[r, ] * scale, type = "extrusion", ...)
#'     })
#' }
#' plot3d(trj3d$x, trj3d$y, trj3d$z, type = 'l')
#' Acc3DArrows(Traj3DAcceleration(trj3d), col = 2)
#' }
#'
#' @export
Traj3DAcceleration <- function(trj3d) {
  .checkTrajHasTime(trj3d)

  # Note that there's no point in calculating backward or forward differences as
  # they provide no different information

  # If we were guaranteed a constant step time, h, we could use the more elegant
  # ax <- stats::filter(trj$x, c(1, -2, 1)) / h^2

  x <- trj3d$x
  y <- trj3d$y
  z <- trj3d$z
  h <- diff(trj3d$time)

  # Calculate velocities using forward/backward diffs
  vx <- diff(x) / h
  vy <- diff(y) / h
  vz <- diff(z) / h
  # calculate acceleration from velocity and time
  h_2 <- .sumPairs(h / 2)
  ax <- diff(vx) / h_2
  ay <- diff(vy) / h_2
  az <- diff(vz) / h_2

  acc <- rbind(c(NA, NA, NA),
               cbind(ax, ay, az),
               c(NA, NA, NA))
  attr(acc, "trj3d") <- trj3d
  acc
}


#' Velocity of a trajectory
#'
#' The velocity, as a 3-dimensional vector, is approximated at each point of the
#' trajectory using first-order finite differences. Central, forward or backward
#' differences can be used. Central differences yield a more accurate
#' approximation if the velocity is smooth. As a practical guide, if velocity
#' doesn't change much between steps, use central differences. If it changes
#' substantially (and not just as an artifact of recording noise), then use
#' either forward or backward differences.
#'
#' Intuitively, think of the central difference velocity at a point as the mean
#' of the velocities of the two adjacent steps. Forward difference velocity is
#' the velocity of the step starting at the point. Backward difference is the
#' velocity of the step ending at the point.
#'
#' Speed (i.e. the magnitude of the velocity) can be derived from velocity by
#' calling \code{\link{Traj3DSpeed}}.
#'
#' @param trj3d 3-dimensional trajectory whose velocity is to be calculated.
#' @param diff Type of difference to be calculated, one of "central" (the
#'   default), "forward" or "backward".
#'
#' @return A numeric matrix. Each row represents the 3D velocity vector at each
#'   point along the trajectory, with `x`, `y` and `z` columns. If \code{diff}
#'   is \code{"central"}, the first and last velocity values will be \code{NA}
#'   since velocity cannot be calculated for them. If \code{diff} is
#'   \code{"forward"}, the last value will be NA, and if \code{diff} is
#'   \code{"backward"}, the first value will be NA.
#'
#' @seealso \code{\link{Traj3DSpeed}} for calculating scalar speed, which is the
#'   magnitude of the velocity vector at each step;
#'   \code{\link{Traj3DResampleTime}} and \code{\link{Traj3DRediscretize}} to
#'   resample a trajectory to fixed time or length steps; Finite differences on
#'   \href{https://en.wikipedia.org/wiki/Finite_difference}{Wikipedia}.
#'
#' @export
Traj3DVelocity <- function(trj3d, diff = c("central", "forward", "backward")) {
  .checkTrajHasTime(trj3d)
  diff <- match.arg(diff)

  x <- trj3d$x
  y <- trj3d$y
  z <- trj3d$z
  h <- diff(trj3d$time)

  if (diff == "central") {

    # Central diffs (this is the "double-interval" central difference)
    # Sum time (h) for each adjacent pair of steps
    dt <- c(NA, .sumPairs(h), NA)
    vx <- stats::filter(x, c(1, 0, -1)) / dt
    vy <- stats::filter(y, c(1, 0, -1)) / dt
    vz <- stats::filter(z, c(1, 0, -1)) / dt
  } else {

    # Forward or backward diffs, variable step times
    vx <- diff(x) / h
    vy <- diff(y) / h
    vz <- diff(z) / h
    if (diff == "forward") {
      # Forward diffs - speed at last point is unknown
      vx <- c(vx, NA)
      vy <- c(vy, NA)
      vz <- c(vz, NA)
    } else {
      # Backward diffs - speed at first point is unknown
      vx <- c(NA, vx)
      vy <- c(NA, vy)
      vz <- c(NA, vz)
    }
  }

  vel <- cbind(vx = vx, vy = vy, vz = vz)
  attr(vel, "trj3d") <- trj3d
  vel
}

#' Speed along a 3-dimensional trajectory
#'
#' Speed is calculated as the magnitude of velocity. The returned speed will
#' contain leading and/or trailing `NA` values, depending on the type of
#' differences used to calculate velocity.
#'
#' If the trajectory has constant time steps, then the average speed of the
#' trajectory is `mean(Traj3DSpeed(Traj3DVelocity(trj)), na.rm = TRUE)`. A
#' trajectory can be resampled by \code{\link{Traj3DResampleTime}} so that it
#' has constant time steps.
#'
#' This function is implemented to simply interpret each row of a matrix as a 3D
#' vector, and return their lengths. Accordingly, this function can also be used
#' to calculate the magnitude of acceleration as returned by
#' \code{\link{Traj3DAcceleration}}.
#'
#' @param vel Velocity of a 3-dimensional trajectory, as returned by
#'   \code{\link{Traj3DVelocity}}.
#'
#' @returns Numeric vector with speed at each point along the trajectory.
#'
#' @seealso \code{\link{Traj3DVelocity}} to calculate the 3D velocity,
#'   \code{\link{Traj3DAcceleration}} to calculate the 3D acceleration,
#'   \code{\link{Traj3DLength}}, \code{\link{TrajDuration}},
#'   \code{\link{Traj3DResampleTime}}
#'
#' @examples
#' \dontrun{
#' # Get the velocity along the trajectory
#' speed <- Traj3DSpeed(Traj3DVelocity(trj3d))
#'
#' # Alternative method to calculate mean speed of a portion of a trajectory
#' si <- 10
#' ei <- 400
#' # Mean speed is displacement divided by time
#' speedMean <- Traj3DLength(trj3d, startIndex = si, endIndex = ei) /
#'                 TrajDuration(trj3d, startIndex = si, endIndex = ei)
#' }
#'
#' @export
Traj3DSpeed <- function(vel) {
  # This is faster than using dist or norm.
  # See https://stackoverflow.com/a/63763823 to explain scaling (it prevents overflow)
  apply(X = vel, MARGIN = 1, FUN = function(x) {
    mx <- max(abs(x))
    if (!is.na(mx) && mx == 0)
      0
    else
      mx * sqrt(sum((x / mx)^2))
  })
}
