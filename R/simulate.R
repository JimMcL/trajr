# ---- Trajectory simulation helper functions ----

# Functions for handling bounded regions

#' Test whether each of the points in a trajectory lie inside a polygon
#'
#' Simply a wrapper around \code{\link[sp]{point.in.polygon}}. The \code{sp}
#' package must be installed for this function to be called. \code{sp} is not
#' automatically installed as a dependency of trajr.
#'
#' @param trj Trajectory to test
#' @param boundary A polygon defining the region to be tested against. Can be
#'   any structure that \code{\link{xy.coords}} can handle, such as a data frame
#'   with \code{x} and \code{y} columns.
#'
#' @return Integer array with a value for each point in the trajectory. Values
#'   are: 0: point is strictly exterior to boundary; 1: point is strictly
#'   interior to boundary; 2: point lies on the relative interior of an edge of
#'   boundary; 3: point is a vertex of boundary
#'
#' @seealso \code{\link[sp]{point.in.polygon}}, \code{\link{xy.coords}}
#'
#' @examples
#' # Square arena
#' boundary <- data.frame(x = c(-10, 10, 10, -10), y = c(-10, -10, 10, 10))
#'
#' # Generate a random trajectory
#' set.seed(1)
#' trj <- TrajGenerate(n = 10, stepLength = 2, angularErrorSd = .15)
#'
#' # Test which points lie inside the boundary
#' print(TrajInPolygon(trj, boundary))
#' ## [1] 1 1 1 1 1 1 0 0 0 0 0
#'
#' @export
TrajInPolygon <- function(trj, boundary) {
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package \"sp\" is needed for the TrajInPolygon function to work. Please install it.",
         call. = FALSE)
  }

  poly <- grDevices::xy.coords(boundary)
  sp::point.in.polygon(trj$x, trj$y, poly$x, poly$y)
}

#' Split a trajectory into multiple sections
#'
#' Every point in \code{trj} will belong to exactly one of the returned
#' sections. Note that this function will happily create single point
#' trajectories.
#'
#' @param trj The trajectory to be split
#' @param idx Indices of splits. Each new section starts at one of these
#'   indices.
#'
#' @return A list containing one or more trajectories. The first trajectory in
#'   the list contains the first points from \code{trj}. Remaining trajectories
#'   contain the points starting from each of the \code{idx} values, in
#'   ascending order.
#'
#' @seealso \code{\link{TrajMerge}}, \code{\link{TrajSplitAtFirstCrossing}}
#'
#' @export
TrajSplit <- function(trj, idx) {
  # Ignore indices below 2 or > length
  idx <- idx[idx > 1 & idx <= nrow(trj)]
  if (length(idx) < 1)
    # No breaks, just return the original trajectory
    return(list(trj))

  # Calculate section bounds
  idx <- sort(unique(idx))
  starts <- c(1, idx)
  ends <- c(idx - 1, nrow(trj))

  # Create a new trajectory for each section
  lapply(seq_along(starts), function(i) {
    .fillInTraj(trj[starts[i]:ends[i], ])
  })
}

#' Combine multiple trajectories into a single whole trajectory
#'
#' This is the inverse of \code{\link{TrajSplit}}.
#'
#' @param parts A list containing one or more trajectories. The trajectories are
#'   concatenated together in order.
#'
#' @return A single trajectory.
#'
#' @seealso \code{\link{TrajSplit}}
#'
#' @examples
#' trj <- TrajGenerate(n = 20)
#' ntrj <- TrajMerge(TrajSplit(trj, c(3, 9, 20)))
#' print(all(trj == ntrj))
#' ## [1] TRUE
#'
#' @export
TrajMerge <- function(parts) {
  .fillInTraj(do.call(rbind, parts))
}

#' Split a trajectory into two parts, separated at the first boundary crossing
#'
#' This is basically a wrapper around \code{\link{TrajInPolygon}} and
#' \code{\link{TrajSplit}}.
#'
#' @param trj The trajectory to split.
#' @param boundary A polygon defining the boundary.  Can be any structure that
#'   \code{\link{xy.coords}} can handle, such as a data frame with \code{x} and
#'   \code{y} columns.
#'
#' @return A list with 1 or 2 elements. If \code{trj} lies entirely inside or
#'   outside \code{boundary}, then the list simply contains \code{trj}. If
#'   \code{trj} crosses the boundary, then the list contains 2 trajectories. The
#'   first is the longest part of \code{trj} that lies entirely inside or
#'   outside \code{boundary}, and the second is the remainder of \code{trj}.
#'
#' @seealso \code{\link{TrajInPolygon}}, \code{\link{TrajSplit}}
#'
#' @examples
#' # Square arena
#' boundary <- data.frame(x = c(-10, 10, 10, -10), y = c(-10, -10, 10, 10))
#'
#' # Generate a random trajectory
#' set.seed(1)
#' trj <- TrajGenerate(n = 8, stepLength = 3, angularErrorSd = .4)
#' # Split the trajectory where it crosses the boundary
#' l <- TrajSplitAtFirstCrossing(trj, boundary)
#'
#' # Plot the boundary and the two trajectories
#' plot(NULL, xlim = range(c(boundary$x, trj$x)), ylim = range(c(boundary$y, trj$y)), asp = 1)
#' polygon(boundary, border = "brown", lwd = 2)
#' lines(l[[1]], col = "#2040ff80", lwd = 3)
#' lines(l[[2]], col = "#ff204080", lwd = 3)
#'
#' @export
TrajSplitAtFirstCrossing <- function(trj, boundary) {
  b <- TrajInPolygon(trj, boundary)
  rl <- rle(b)
  # First trajectory is all inside (or all outside) the boundary. The first
  # point of the 2nd trajectory is outside
  TrajSplit(trj, rl$lengths[1] + 1)
}

