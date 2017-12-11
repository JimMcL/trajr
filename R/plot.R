# ---- Trajectory plotting ----

.drawTurningAngles <- function(x, turning.angles) {
  # There are n steps, but n+1 coordinates
  n <- nrow(x) - 1
  steps <- x[1:n,]
  angles <- x[2:(n+1),]
  meanStepLength <- TrajMeanStepLength(x)
  segLen <- 0.8 * meanStepLength
  textDisplacement <- 0.3 * meanStepLength
  labels <- parse(text= paste("Delta[", 1:n, "]", sep=""))

  if (tolower(turning.angles) == "directed") {
    # Plot angles which represent angular errors, which reset at each step
    graphics::segments(steps$x, steps$y, steps$x + .8 * meanStepLength, steps$y, col = "darkgrey", lty = 2)

    textAngle <- Arg(angles$displacement) +
      ifelse(Arg(x$displacement[1]) < Arg(angles$displacement), pi / 4, -pi / 6)
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
      ifelse(Arg(steps$displacement) < Arg(angles$displacement), pi / 4, -pi / 6)
    graphics::text(steps$x + textDisplacement * cos(textAngle), steps$y + textDisplacement * sin(textAngle),
                   labels = labels)

    plotrix::draw.arc(steps$x, steps$y,
                      angle1 = Arg(steps$displacement), angle2 = Arg(angles$displacement),
                      radius = 0.4 * meanStepLength)

  } else {
    stop(sprintf("Invalid turning.angles (%s), must be one of 'random' or 'directed'", turning.angles))
  }
}

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
                            xlim = grDevices::extendrange(x$x), ylim = grDevices::extendrange(x$y),
                            xlab = "x", ylab = "y",
                            asp = 1, ...) {
  if (!add) {
    graphics::plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, asp = asp, ...)
  }
  graphics::lines(y ~ x, data = x, type = type, ...)
  if (draw.start.pt)
    graphics::points(x$x[1], x$y[1], pch = 16, cex = .8)

  if (!is.null(turning.angles)) {
    .drawTurningAngles(x, turning.angles)
  }
}
