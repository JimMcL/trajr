# Contains functions to calculate various measures of trajectory straightness or tortuosity

.rad2deg <- function(rad) { rad * 180 / pi}
.deg2rad <- function(deg) { deg * pi / 180}


#' Straightness of a Trajectory
#'
#' Calculates the straightness index of a trajectory, \eqn{D / L}, where
#' \code{D} is the beeline distance between the first and last points in the
#' trajectory,and \code{L} is the path length travelled.
#'
#' @param trj Trajectory to calculate straightness of.
#' @return The straightness index of \code{trj}, which is a value between 0
#'   (infinitely tortuous) to 1 (a straight line).
#'
#' @seealso \code{\link{TrajDistance}} for trajectory distance, and
#'   \code{\link{TrajLength}} for trajectory path length.
#'
#' @references Batschelet, E. (1981). Circular statistics in biology. ACADEMIC
#' PRESS, 111 FIFTH AVE., NEW YORK, NY 10003, 1981, 388.
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's
#' path. Journal of Theoretical Biology, 229(2), 209-220.
#' doi:10.1016/j.jtbi.2004.03.016

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
#' SD = mean(TrajDirectionalChange(trj))
#' SDDC = sd(TrajDirectionalChange(trj))
#'
#' @references Kitamura, T., & Imafuku, M. (2015). Behavioural mimicry in flight
#'   path of Batesian intraspecific polymorphic butterfly Papilio polytes.
#'   Proceedings of the Royal Society B: Biological Sciences, 282(1809).
#'   doi:10.1098/rspb.2015.0483
TrajDirectionalChange <- function(trj, nFrames = 1) {
  .rad2deg(TrajAngles(trj, nFrames)) / diff(trj$displacementTime, nFrames)
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

  head(maxima, nm)
}

# Returns the mean segment length of a trajectory
.TrajSegLen <- function(trj) mean(Mod(trj))


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
#' # Resample to fixed path length
#' resampled <- TrajRediscretize(trj, .005)
#' # Calculate direction autocorrelation for resampled trajectory
#' corr <- TrajDirectionAutocorrelations(resampled)
#' # Extract first local minimum from autocorrelation
#' minPt <- TrajDAFindFirstMinimum(corr, 50)
#'
#' # Plot the autocorrelation function
#' plot(corr)
#' # Plot a red dot with a black outline at the first minimum
#' points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
#' points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
#'
NULL

#' @rdname DAMinMax
TrajDAFindFirstMinimum <- function(corr, windowSize = 10) {
  # Ignore local minimum if it's the end of the track
  windowSize <- min(length(corr$C) - 1, windowSize)
  minima <- .JLocalMaxima(-corr$C, windowSize, endIndex = length(corr$C) - windowSize)
  if (length(minima) > 0) {
    c(deltaS = corr$deltaS[minima][1], C = corr$C[minima][1])
  }
}

#' @rdname DAMinMax
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
#' result. \code{trj} must have a constant step length (see
#' \code{\link{TrajDirectionAutocorrelations}} for further details).
#'
#' @param trj Trajectory to be plotted.
#' @param deltaSMax Maximum delta s to calculated, see
#'   \code{\link{TrajDirectionAutocorrelations}} for details.
#' @param firstMinWindowSize If not NULL, specifies a window size used to
#'   calculate the first local minimum, which is then plotted as a point.
#' @param type,xlab,ylab Defaults for plotting.
TrajPlotDirectionAutocorrelations <- function(trj,
                                              deltaSMax = round(nrow(trj) / 4),
                                              firstMinWindowSize = 10,
                                              type = 'l',
                                              ylab = expression('C('*Delta*s*')'),
                                              xlab = expression(Delta*s),
                                              ...) {
  corr <- TrajDirectionAutocorrelations(trj, deltaSMax)

  # Plot the autocorrelation function
  plot(corr, type = type, xlab = xlab, ylab = ylab, ...)

  # Optionally plot first minimum
  if (!is.null(firstMinWindowSize)) {
    # Extract first local minimum from autocorrelation
    minPt <- TrajDAFindFirstMinimum(corr, firstMinWindowSize)
    # Plot a red dot with a black outline at the first minimum
    points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
    points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
  }
}

# Sinuosity ########################################################################################

#' Sinuosity of a trajectory
#'
#' Calculates the sinuosity of a trajectory as defined by Benhamou (2004), which
#' is: \deqn{S = 1.18\sigma_{q} / \sqrtq}.
#'
#' @param trj Trajectory to calculate sinuosity of.
#' @return The sinuosity of \code{trj}.
#'
TrajSinuosity <- function(trj) {
  # Discard initial 0-length segment
  trj <- tail(trj$displacement, -1)
  segLen <- .TrajSegLen(trj)
  1.18 * sd(diff(Arg(trj))) / sqrt(segLen)
}


# E-MAX ########################################################################################

#' E-max from Cheung et al., (2007)
#'
#' Not yet tested
TrajEmax <- function(trj, eMaxB = FALSE) {

  .beta <- function(points) {
    # Calculate difference in angle for every pair of segments which are deltaS apart,
    # take cos, then mean
    mean(cos(diff(Arg(points))))
  }

  b <- .beta(trj$polar)

  # If it's E max b, multiply by mean path length
  f <- ifelse(eMaxB, .TrajSegLen(trj), 1)

  f * b / (1 - b)
}
