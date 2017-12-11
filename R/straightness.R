# Contains functions to calculate various measures of trajectory straightness or tortuosity

.rad2deg <- function(rad) { rad * 180 / pi}
.deg2rad <- function(deg) { deg * pi / 180}


#' Straightness of a Trajectory
#'
#' Calculates the straightness index of a trajectory, \eqn{D / L}, where
#' \code{D} is the beeline distance between the first and last points in the
#' trajectory,and \code{L} is the path length travelled (Batschelet, 1981).
#' Benhamou (2004) considers the straightness index to be a reliable measure of
#' the efficiency of a directed walk, but inapplicable to random trajectories.
#'
#' @param trj Trajectory to calculate straightness of.
#' @return The straightness index of \code{trj}, which is a value between 0
#'   (infinitely tortuous) to 1 (a straight line).
#'
#' @seealso \code{\link{TrajDistance}} for trajectory distance, and
#'   \code{\link{TrajLength}} for trajectory path length.
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
TrajStraightness <- function(trj) {
  TrajDistance(trj) / TrajLength(trj)
}

#' Directional change (DC)
#'
#' Calculates the time variation of directional change (DC) of a trajectory
#' \emph{sensu} Kitamura & Imafuku (2015). Directional change is defined as the
#' angular change (in degrees) between any two points in the trajectory, divided
#' by the time difference between the two points.
#'
#' This function returns the DC for each pair of consecutive points. Kitamura &
#' Imafuku (2015) used the mean and the standard deviation of DC for portions of
#' trajectories as index values of nonlinearity and irregularity respectively.
#'
#' @param trj Track to calculate DC for.
#' @param nFrames Frame delta to process: if 1, every frame is processed, if 2,
#'   every 2nd frame is processed, and so on. Default is 1.
#' @return The directional change (DC) between every pair of consecutive points
#'   in the trajectory, i.e. the returned vector will have length
#'   \code{(nrow(trj) - 1)}.
#'
#' @examples
#' set.seed(42)
#' trj <- TrajGenerate()
#' SD = mean(TrajDirectionalChange(trj))
#' SDDC = sd(TrajDirectionalChange(trj))
#'
#' @references Kitamura, T., & Imafuku, M. (2015). Behavioural mimicry in flight
#'   path of Batesian intraspecific polymorphic butterfly Papilio polytes.
#'   Proceedings of the Royal Society B: Biological Sciences, 282(1809).
#'   doi:10.1098/rspb.2015.0483
#'
#' @export
TrajDirectionalChange <- function(trj, nFrames = 1) {
  # Calculating this way is almost 1 order of magnitude faster than using the documented equation
  abs(.rad2deg(TrajAngles(trj, nFrames))) / diff(trj$displacementTime, 2 * nFrames)
}


# Directional autocorrelation functions ##########################################################

# Autocorrelation based on (Shamble et al., 2017)
#
# I have tried to use variable names which are suggestive
# of the names in the article

# Private functions ====================================================

# Returns indices of local maxima.
# To obtain local minima, call .JLocalMaxima(-x)
#
# @param v vector of values.
# @param window Number of points on each side which defines what counts as a local maxima.
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
#' segments.
#'
#' @param trj The trajectory to calculate the directional autocorrelations for.
#' @param deltaSMax Maximum delta s to calculate, default is \eqn{1/4} the
#'   number of segments in the trajectory.
#' @return a data frame with 2 columns, \code{deltaS} and \code{C}.
#'
#' @seealso TrajDAFindFirstMinimum, TrajPlotDirectionAutocorrelations
#'
#' @references Shamble, P. S., Hoy, R. R., Cohen, I., & Beatus, T. (2017).
#'   Walking like an ant: a quantitative and experimental approach to
#'   understanding locomotor mimicry in the jumping spider Myrmarachne
#'   formicaria. Proceedings of the Royal Society B: Biological Sciences,
#'   284(1858). doi:10.1098/rspb.2017.0308
#'
#' @export
TrajDirectionAutocorrelations <- function(trj, deltaSMax = round(nrow(trj) / 4)) {

  # The guts of the autocorrelation function
  # Calculates autocorrelation for a single delta s
  .deltaCorr <- function(deltaS, trj) {
    # Calculate difference in angle for every pair of segments which are deltaS apart,
    # take cos, then mean
    c <- sapply(seq(1, length.out = deltaS), function(offset) {
      t <- trj[offset:nrow(trj),]
      cos(TrajAngles(t, deltaS))
    })
    mean(unlist(c))
  }

  deltaSs <- 1:deltaSMax
  data.frame(deltaS = deltaSs,
             C = sapply(deltaSs, .deltaCorr, trj))
}

#' First direction autocorrelation minimum/maximum
#'
#' Determines the coordinates of the first local minimum/maximum of \code{C} in
#' the direction autocorrelation function of a trajectory as returned by
#' \code{\link{TrajDirectionAutocorrelations}}. The end point is excluded from
#' consideration as a minimum, similarly the start point will not be returned as a
#' maximum.
#'
#' @param corr Direction autocorrelation of a trajectory.
#' @param windowSize Size of window used to define what constitutes a local
#'   mimimum/maximum.
#' @return Numeric vector with 2 values, \code{deltaS} and \code{C}, or else
#'   NULL if there is no local minimum/maximum.
#' @name DAMinMax
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

#' @rdname DAMinMax
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

#' @rdname DAMinMax
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

#' Plot direction autocorrelation function
#'
#' Calculate the direction autocorrelation for a trajectory, then plot the
#' result, with a dot at the first local minimum. \code{trj} must have a
#' constant step length (see \code{\link{TrajDirectionAutocorrelations}} for
#' further details).
#'
#' @param trj Trajectory to be plotted.
#' @param deltaSMax Maximum delta s to be calculated, see
#'   \code{\link{TrajDirectionAutocorrelations}} for details.
#' @param firstMinWindowSize If not NULL, specifies a window size used to
#'   calculate the first local minimum, which is then plotted as a point.
#' @param type,xlab,ylab Defaults for plotting.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @export
TrajPlotDirectionAutocorrelations <- function(trj,
                                              deltaSMax = round(nrow(trj) / 4),
                                              firstMinWindowSize = 10,
                                              type = 'l',
                                              ylab = expression('C('*Delta*s*')'),
                                              xlab = expression(Delta*s),
                                              ...) {
  corr <- TrajDirectionAutocorrelations(trj, deltaSMax)

  # Plot the autocorrelation function
  graphics::plot(corr, type = type, xlab = xlab, ylab = ylab, ...)

  # Optionally plot first minimum
  if (!is.null(firstMinWindowSize)) {
    # Extract first local minimum from autocorrelation
    minPt <- TrajDAFindFirstMinimum(corr, firstMinWindowSize)
    # Plot a red dot with a black outline at the first minimum
    graphics::points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
    graphics::points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
  }
}

# Sinuosity ########################################################################################

#' Sinuosity of a trajectory
#'
#' Calculates the sinuosity of a trajectory as defined by Bovet & Benhamou
#' (1988), which is: \eqn{S = 1.18\sigma / \sqrtq} where \eqn{\sigma} is the
#' standard deviation of the step turning angles and \eqn{q} is the mean step
#' length.
#'
#' @param trj Trajectory to calculate sinuosity of.
#' @param compass.direction if not \code{NULL}, turning angles are calculated
#'   for a directed walk, assuming the specified compass direction (in radians).
#'   Otherwise, a random walk is assumed.
#' @return The sinuosity of \code{trj}.
#'
#' @seealso \code{\link{TrajAngles}} for the turning angles in a trajectory,
#'   \code{\link{TrajMeanStepLength}} for the mean step length.
#'
#' @references
#'
#' Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals' movements
#' using a correlated random walk model. Journal of Theoretical Biology, 131(4),
#' 419-433. doi:10.1016/S0022-5193(88)80038-9
#'
#' @export
TrajSinuosity <- function(trj, compass.direction = NULL) {
  segLen <- TrajMeanStepLength(trj)
  1.18 * stats::sd(TrajAngles(trj, compass.direction = compass.direction)) / sqrt(segLen)
}


# E-MAX ########################################################################################

#' Trajectory straightness index, E-max
#'
#' \eqn{E[max]} is a single-valued measure of straightness defined by (Cheung,
#' Zhang, Stricker, & Srinivasan, 2007).
#'
#' @param trj Trajectory to be analysed.
#' @param eMaxB If TRUE, calculates and returns E-max b, otherwise return E-max.
#' @param compass.direction if not \code{NULL}, turning angles are calculated
#'   for a directed walk, assuming the specified compass direction (in radians).
#'   Otherwise, a random walk is assumed.
#' @return E-max for \code{trj}.
#'
#' @references Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007).
#'   Animal navigation: the difficulty of moving in a straight line. Biological
#'   Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' @export
TrajEmax <- function(trj, eMaxB = FALSE, compass.direction = NULL) {

  # E(cos(angles)) = mean(cos(angles))
  b <- mean(cos(TrajAngles(trj, compass.direction = compass.direction)))

  # If it's E max b, multiply by mean step length
  f <- ifelse(eMaxB, TrajMeanStepLength(trj), 1)

  f * b / (1 - b)
}
