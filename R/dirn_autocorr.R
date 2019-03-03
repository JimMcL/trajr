# Directional autocorrelation functions ##########################################################

# Autocorrelation based on (Shamble et al., 2017)
#
# I have tried to use variable names which are suggestive
# of the names in the article

# Private functions ====================================================

# Returns indices of local maxima.
#
# To obtain local minima, call .JLocalMaxima(-x).
#
# A maximimum is defined here as a point equal to the greatest value within a
# window. Hence, two or more equal contiguous values will count as maxima if
# there are no higher vaues within each window.
#
# @param v vector of values.
# @param window Number of points on each side which defines what counts as a
#   local maxima.
# @param startIndex Index of first point which can qualify as a maximum.
# @param endIndex Index of last point which can qualify as a maximum.
.JLocalMaxima <- function(v, window = 1, startIndex = 1, endIndex = length(v))
{
  getWindow <- function(i) {
    # Don't try to look past the ends of the data
    si <- max(1, i - window)
    ei <- min(length(v), i + window)
    v[si : ei]
  }

  maxima <- numeric(length(v) / 2)
  nm <- 0
  for (i in startIndex:endIndex) {
    # Is this point a maximum?
    if (v[i] == max(getWindow(i))) {
      nm <- nm + 1
      maxima[nm] <- i
    }
  }

  utils::head(maxima, nm)
}


# Public functions =====================================================

#' Direction autocorrelation
#'
#' Calculates the autocorrelation of the track for \eqn{\Delta}s ranging from 1
#' to \code{deltaSMax}, based on Shamble et al. (2017). \code{trj} must have a
#' constant step length (see \code{\link{TrajRediscretize}}) i.e. all segments
#' in the trajectory must be the same length. deltaS is specified in number of
#' segments. Call \code{\link{TrajDAFindFirstMinimum}} to locate the first local
#' minimum which may be used to characterise directional periodicity in a
#' trajectory (note that the first local minimum may not exist).
#'
#' @param trj The trajectory to calculate the directional autocorrelations for.
#' @param deltaSMax Maximum delta s to calculate, default is \eqn{1/4} the
#'   number of segments in the trajectory.
#' @return A data frame with class \code{TrajDirectionAutocorrelations} and 2
#'   columns, \code{deltaS} and \code{C}. Plotting this object displays a graph
#'   of the direction autocorrelation function, optionally with the location of
#'   the first local minimum marked
#'
#' @seealso \code{\link{TrajDAFindFirstMinimum}},
#'   \code{\link{plot.TrajDirectionAutocorrelations}}
#'
#' @references Shamble, P. S., Hoy, R. R., Cohen, I., & Beatus, T. (2017).
#'   Walking like an ant: a quantitative and experimental approach to
#'   understanding locomotor mimicry in the jumping spider Myrmarachne
#'   formicaria. Proceedings of the Royal Society B: Biological Sciences,
#'   284(1858). doi:10.1098/rspb.2017.0308
#'
#' @export
TrajDirectionAutocorrelations <- function(trj, deltaSMax = round(nrow(trj) / 4)) {

  deltaSs <- 1:round(deltaSMax)

  # The guts of the autocorrelation function
  # Calculates autocorrelation for a single delta s
  .deltaCorr <- function(deltaS, trj) {
    # Calculate difference in angle for every pair of segments which are deltaS apart,
    # take cos, then mean
    c <- sapply(deltaSs, function(offset) {
      t <- trj[offset:nrow(trj),]
      cos(TrajAngles(t, deltaS))
    })
    mean(unlist(c))
  }

  r <- data.frame(deltaS = deltaSs,
                  C = sapply(deltaSs, .deltaCorr, trj))
  class(r) <- c("TrajDirectionAutocorrelations", class(r))
  r
}

#' First direction autocorrelation minimum/maximum
#'
#' Determines the coordinates of the first local minimum/maximum of \code{C} in
#' the direction autocorrelation function of a trajectory as returned by
#' \code{\link{TrajDirectionAutocorrelations}}. The end point is excluded from
#' consideration as a minimum, similarly the start point will not be returned as
#' a maximum. if the trajectory does not osciallate in direction, there will not
#' be a local minimum/maximum, and \code{NULL} is returned.
#'
#' @param corr a \code{TrajDirectionAutocorrelations} object, i.e. the direction
#'   autocorrelation of a trajectory.
#' @param windowSize Size of window used to define what constitutes a local
#'   mimimum/maximum.
#' @return Numeric vector with 2 values, \code{deltaS} and \code{C}, or NULL if
#'   there is no local minimum/maximum.
#' @name TrajDAMinMax
#'
#' @seealso \code{\link{TrajDirectionAutocorrelations}}
#'
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate(600, angularErrorSd = 1)
#' smoothed <- TrajSmoothSG(trj, 3, 11)
#'
#' # Resample to fixed path length
#' resampled <- TrajRediscretize(smoothed, 1)
#' # Calculate direction autocorrelation for resampled trajectory
#' corr <- TrajDirectionAutocorrelations(resampled, 100)
#' # Extract first local minimum from autocorrelation
#' minPt <- TrajDAFindFirstMinimum(corr, 20)
#'
#' # Plot the autocorrelation function
#' plot(corr, type ='l')
#' # Plot a red dot with a black outline at the first minimum
#' points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
#' points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
#'
NULL

#' @rdname TrajDAMinMax
#'
#' @export
TrajDAFindFirstMinimum <- function(corr, windowSize = 10) {
  # Ignore local minimum if it's the end of the track
  windowSize <- min(length(corr$C) - 1, windowSize)
  minima <- .JLocalMaxima(-corr$C, windowSize, endIndex = length(corr$C) - windowSize)
  if (length(minima) > 0) {
    c(deltaS = corr$deltaS[minima][1], C = corr$C[minima][1])
  }
}

#' @rdname TrajDAMinMax
#'
#' @export
TrajDAFindFirstMaximum <- function(corr, windowSize = 10) {
  windowSize <- min(length(corr$C) - 1, windowSize)
  # Ignore local maxima if it's the start of the track
  maxima <- .JLocalMaxima(corr$C, windowSize, startIndex = 2)
  if (length(maxima) > 1) {
    c(deltaS = corr$deltaS[maxima][1], C = corr$C[maxima][1])
  }
}

#' Plot method for direction autocorrelation
#'
#' The \code{plot} method for \code{TrajDirectionAutocorrelations} objects.
#' Plots the direction autocorrelation function as returned by a call to
#' \code{link{TrajDirectionAutocorrelations}}, with a optional dot at the first
#' local minimum.
#'
#' @param x Trajectory to be plotted.
#' @param firstMinWindowSize If not NULL, specifies a window size used to
#'   calculate the first local minimum, which is then plotted as a point.
#' @param type,xlab,ylab Defaults for plotting.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @export
plot.TrajDirectionAutocorrelations <- function(x,
                                               firstMinWindowSize = 10,
                                               type = 'l',
                                               ylab = expression('C('*Delta*s*')'),
                                               xlab = expression(Delta*s),
                                               ...) {
  # Plot the autocorrelation function
  graphics::plot(x[,1], x[,2], type = type, xlab = xlab, ylab = ylab, ...)

  # Optionally plot first minimum
  if (!is.null(firstMinWindowSize)) {
    # Extract first local minimum from autocorrelation
    minPt <- TrajDAFindFirstMinimum(x, firstMinWindowSize)
    # Plot a red dot with a black outline at the first minimum
    graphics::points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
    graphics::points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
  }
}

