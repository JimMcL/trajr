# ---- Trajectory plotting ----

# ---- Private functions ----

.drawTurningAngles <- function(x, turning.angles, text.displacement.angles = c(pi / 6, -pi / 8)) {

  .textAngle <- function(positive) ifelse(positive, text.displacement.angles[1], text.displacement.angles[2])

  # There are n steps, but n+1 coordinates
  n <- nrow(x) - 1
  steps <- x[1:n,]
  angles <- x[2:(n+1),]
  meanStepLength <- mean(TrajStepLengths(x))
  segLen <- 0.8 * meanStepLength
  textDisplacement <- 0.3 * meanStepLength
  labels <- parse(text= paste("Delta[", 1:n, "]", sep=""))

  if (tolower(turning.angles) == "directed") {
    # Plot angles which represent angular errors, which reset at each step
    graphics::segments(steps$x, steps$y, steps$x + .8 * meanStepLength, steps$y, col = "darkgrey", lty = 2)

    textAngle <- Arg(angles$displacement) +
      .textAngle(Arg(x$displacement[1]) < Arg(angles$displacement))
    graphics::text(steps$x + textDisplacement * cos(textAngle), steps$y + textDisplacement * sin(textAngle),
                   labels = labels)

    plotrix::draw.arc(steps$x, steps$y,
                      angle1 = Arg(x$displacement[1]), angle2 = Arg(angles$displacement),
                      radius = 0.4 * meanStepLength)

  } else if (tolower(turning.angles) == "random") {
    # Plot angles which represent angular errors, which accumulate
    graphics::segments(steps$x, steps$y,
                       steps$x + segLen * cos(Arg(steps$displacement)),
                       steps$y + segLen * sin(Arg(steps$displacement)),
                       col = "darkgrey", lty = 2)

    textAngle <- Arg(angles$displacement) +
      .textAngle(Arg(steps$displacement) < Arg(angles$displacement))
    graphics::text(steps$x + textDisplacement * cos(textAngle), steps$y + textDisplacement * sin(textAngle),
                   labels = labels)

    plotrix::draw.arc(steps$x, steps$y,
                      angle1 = Arg(steps$displacement), angle2 = Arg(angles$displacement),
                      radius = 0.4 * meanStepLength)

  } else {
    stop(sprintf("Invalid turning.angles (%s), must be one of 'random' or 'directed'", turning.angles))
  }
}

# Private
.drawTrajExtras <- function(x, draw.start.pt = TRUE, start.pt.cex = .8, turning.angles = NULL, ...) {
  if (draw.start.pt) {
    graphics::points(x$x[1], x$y[1], pch = 16, cex = start.pt.cex)
  }

  if (!is.null(turning.angles)) {
    .drawTurningAngles(x, turning.angles)
  }
}


# ---- Public functions ----

#' Plot method for trajectories
#'
#' The \code{plot} method for Trajectory objects.
#'
#' @param x An object of class "Trajectory", the trajectory to be plotted.
#' @param add If TRUE, the trajectory is added to the current plot.
#' @param draw.start.pt If TRUE, draws a dot at the start point of the
#'   trajectory.
#' @param start.pt.cex Scale to apply when drawing the start point dot.
#' @param turning.angles If \code{random} or \code{directed}, draws step turning
#'   angles. \code{directed} assumes errors are relative to the first recorded
#'   step angle. \code{random} assumes errors are relative to the previous step.
#' @param xlim,ylim,xlab,ylab,asp plotting parameters with useful defaults.
#' @param ... Additional arguments are passed to \code{\link[graphics]{plot}}.
#'
#' @seealso \code{\link{TrajFromCoords}}
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate(angularErrorSd = 1.3)
#' plot(trj)
#'
#' @importFrom grDevices extendrange
#'
#' @export
plot.Trajectory <- function(x, add = FALSE,
                            draw.start.pt = TRUE, start.pt.cex = TRUE, turning.angles = NULL,
                            xlim = grDevices::extendrange(x$x), ylim = grDevices::extendrange(x$y),
                            xlab = ifelse(is.null(TrajGetUnits(x)), "x", sprintf("x (%s)", TrajGetUnits(x))),
                            ylab = ifelse(is.null(TrajGetUnits(x)), "y", sprintf("y (%s)", TrajGetUnits(x))),
                            asp = 1, ...) {
  if (!add) {
    graphics::plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, asp = asp, ...)
  }
  graphics::lines(x, draw.start.pt = draw.start.pt, start.pt.cex = start.pt.cex, turning.angles = turning.angles, ...)
}

#' Add Trajectory lines to a plot
#'
#' The \code{lines} method for Trajectory objects.
#'
#' @param x An object of class "Trajectory", the trajectory to be plotted.
#' @param draw.start.pt If TRUE, draws a dot at the start point of the
#'   trajectory.
#' @param start.pt.cex Scale to apply when drawing the start point dot.
#' @param turning.angles If \code{random} or \code{directed}, draws step turning
#'   angles. \code{directed} assumes errors are relative to the first recorded
#'   step angle. \code{random} assumes errors are relative to the previous step.
#' @param ... Additional arguments are passed to \code{\link[graphics]{lines}}.
#'
#' @export
lines.Trajectory <- function(x, draw.start.pt = TRUE, start.pt.cex = 0.8, turning.angles = NULL, ...) {
  graphics::lines(y ~ x, data = x, ...)
  .drawTrajExtras(x, draw.start.pt, start.pt.cex, turning.angles, ...)
}

#' Add Trajectory points to a plot
#'
#' The \code{points} method for Trajectory objects.
#'
#' @param x An object of class "Trajectory", the trajectory to be plotted.
#' @param draw.start.pt If TRUE, draws a dot at the start point of the
#'   trajectory.
#' @param turning.angles If \code{random} or \code{directed}, draws step turning
#'   angles. \code{directed} assumes errors are relative to the first recorded
#'   step angle. \code{random} assumes errors are relative to the previous step.
#' @param ... Additional arguments are passed to \code{\link[graphics]{points}}.
#'
#' @export
points.Trajectory <- function(x, draw.start.pt = TRUE, turning.angles = NULL, ...) {
  graphics::points(y ~ x, data = x, ...)
  .drawTrajExtras(x, draw.start.pt, turning.angles, ...)
}
