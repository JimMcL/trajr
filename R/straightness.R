
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
#' @param deltaSMax Maximum delta s to calculate, default is 1/2 the number of
#'   segments in the trajectory.
#' @return a data frame with 2 columns, \code{deltaS} and \code{C}.
#'
#' @references Shamble, P. S., Hoy, R. R., Cohen, I., & Beatus, T. (2017). Walking like an ant: a quantitative and experimental approach to understanding locomotor mimicry in the jumping spider Myrmarachne formicaria. Proceedings of the Royal Society B: Biological Sciences, 284(1858). doi:10.1098/rspb.2017.0308
TrajDirectionAutocorrelations <- function(trj, deltaSMax = round(nrow(trj) / 2)) {

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

#' First direction autocorrelation minimum
#'
#' Determines the coordinates of the first local minimum for \code{C} in the
#' direction autocorrelation function of a trajectory as returned by
#' \code{\link{TrajDirectionAutocorrelations}}. The end point is excluded from
#' consideration as minimum.
#'
#' @param corr Direction autocorrelation of a trajectory.
#' @param windowSize Size of window used to define what constitutes a local mimimum.
#' @return Numeric vector with 2 values, \code{deltaS} and \code{C}, or else
#'   NULL if there is no local minimum (other than the end point).
#'
#' @seealso \code{\link{TrajDAFindFirstMaximum}}, \code{\link{TrajDirectionAutocorrelations}}
#'
#' @examples
#' # Calculate direction autocorrelation for trj
#' corr <- TrajDirectionAutocorrelations(trj)
#' # Extract first local minimum from autocorrelation
#' minPt <- TrajDAFindFirstMinimum(corr, 50)
#'
#' # Plot the autocorrelation function
#' plot(corr)
#' # Plot a red dot with a black outline at the first minimum
#' points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
#' points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
#'
TrajDAFindFirstMinimum <- function(corr, windowSize = 10) {
  # Ignore local minimum if it's the end of the track
  windowSize <- min(length(corr$C) - 1, windowSize)
  minima <- .JLocalMaxima(-corr$C, windowSize, endIndex = length(corr$C) - windowSize)
  if (length(minima) > 0) {
    c(deltaS = corr$deltaS[minima][1], C = corr$C[minima][1])
  }
}

#' First direction autocorrelation maximum
#'
#' Determines the coordinates of the first local maximum for \code{C} in the
#' direction autocorrelation function of a trajectory as returned by
#' \code{\link{TrajDirectionAutocorrelations}}. The start point is excluded from
#' consideration as maximum.
#'
#' @param corr Direction autocorrelation of a trajectory.
#' @param windowSize Size of window used to define what constitutes a local maximum.
#' @return Numeric vector with 2 values, \code{deltaS} and \code{C}, or else
#'   NULL if there is no local minimum (other than the start point).
#'
#' @seealso \code{\link{TrajDAFindFirstMinimum}}, \code{\link{TrajDirectionAutocorrelations}}
#'
#' @examples
#' # Calculate direction autocorrelation for trj
#' corr <- TrajDirectionAutocorrelations(trj)
#' # Extract first local maximum from autocorrelation
#' minPt <- TrajDAFindFirstMaximum(corr, 50)
#'
#' # Plot the autocorrelation function
#' plot(corr)
#' # Plot a red dot with a black outline at the first maximum
#' points(minPt["deltaS"], minPt["C"], pch = 16, col = "red", lwd = 2)
#' points(minPt["deltaS"], minPt["C"], col = "black", lwd = 2)
#'
TrajDAFindFirstMaximum <- function(corr, windowSize = 10) {
  windowSize <- min(length(corr$C) - 1, windowSize)
  # Ignore local maxima if it's the start of the track
  maxima <- .JLocalMaxima(corr$C, windowSize, startIndex = 2)
  if (length(maxima) > 1) {
    c(deltaS = corr$deltaS[maxima][1], C = corr$C[maxima][1])
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

  .Mbeta <- function(points) {
    # Calculate difference in angle for every pair of segments which are deltaS apart,
    # take cos, then mean
    mean(cos(diff(Arg(points))))
  }

  b <- .Mbeta(trj$polar)

  # If it's E max b, multiply by mean path length
  f <- ifelse(eMaxB, .TrajSegLen(trj), 1)

  f * b / (1 - b)
}
