# Contains functions to calculate various measures of trajectory straightness or tortuosity

.rad2deg <- function(rad) { rad * 180 / pi}
.deg2rad <- function(deg) { deg * pi / 180}


#' Mean vector of turning angles
#'
#' Returns the mean vector of the turning angles, as defined by Batschelet,
#' (1981). A unit vector is created for each turning angle in the trajectory,
#' and the centre-of-mass/mean vector is returned.
#'
#' According to Batschelet (1981), \code{r} may serve as a straightness index
#' ranging from 0 to 1, where \code{r} is the length of the mean vector of
#' turning angles of a trajectory with constant step length. Values of \code{r}
#' near 1 indicating straighter paths. Hence, \code{r =
#' Mod(TrajMeanVectorOfTurningAngles(trj))}, assuming that \code{trj} has a
#' constant step length (e.g. has been rediscretized).
#'
#' @param trj Trajectory object.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.

#' @return A complex number \code{r} which represents the mean vector,
#'   \code{Mod(r)} is the length of the mean vector which varies between 0 and
#'   1, \code{Arg(r)} is the angle.
#'
#' @seealso \code{\link{TrajStraightness}}, \code{\link{TrajAngles}},
#'   \code{\link{TrajRediscretize}} for resampling a trajectory to a constant
#'   step length, \code{\link{TrajResampleTime}} for resampling a trajectory to
#'   a constant step time.
#'
#' @references
#'
#' Batschelet, E. (1981). Circular statistics in biology. ACADEMIC PRESS, 111
#' FIFTH AVE., NEW YORK, NY 10003, 1981, 388.
#'
#' @export
TrajMeanVectorOfTurningAngles <- function(trj, compass.direction = NULL) {

  angles <- TrajAngles(trj, compass.direction = compass.direction)

  # The value as defined in Batschelet, value is identical but this calculation
  # is about twice as slow:
  # phi <- atan2(sum(sin(angles)), sum(cos(angles)))
  # r <- sqrt(sum(cos(angles)) ^ 2 + sum(sin(angles)) ^ 2) / length(angles)
  # complex(modulus = r, argument = phi)

  mean(complex(modulus = 1, argument = angles), na.rm = TRUE)
}

#' Straightness of a Trajectory
#'
#' Calculates the straightness index of a trajectory, \eqn{D / L}, where
#' \code{D} is the beeline distance between the first and last points in the
#' trajectory,and \code{L} is the path length travelled (Batschelet, 1981).
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
#' @param trj Trajectory to calculate straightness of.
#' @return The straightness index of \code{trj}, which is a value between 0
#'   (infinitely tortuous) to 1 (a straight line).
#'
#' @seealso \code{\link{TrajDistance}} for trajectory distance (or
#'   displacement), and \code{\link{TrajLength}} for trajectory path length.
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
#' angular change (in degrees) between two steps in the trajectory, divided by
#' the time difference between the two steps.
#'
#' This function returns the DC for each pair of consecutive steps. Kitamura &
#' Imafuku (2015) used the mean and the standard deviation of DC for portions of
#' trajectories as index values of nonlinearity and irregularity respectively.
#'
#' @param trj Track to calculate DC for.
#' @param nFrames Frame delta to process: if 1, every frame is processed, if 2,
#'   every 2nd frame is processed, and so on. Default is 1.
#' @return The directional change (DC) in degrees between every pair of
#'   consecutive segments in the trajectory, i.e. if \code{nFrames} is 1, the
#'   returned vector will have length \code{nrow(trj) - 2}.
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
  .checkTrajHasTime(trj)

  # Calculating this way is about 1 order of magnitude faster than using the
  # documented equation. There is a unit test (in test_basic.R) which checks
  # that this gives the same results as the documented equation
  abs(.rad2deg(TrajAngles(trj, nFrames))) / diff(trj$displacementTime, 2 * nFrames)
}


# Sinuosity ########################################################################################

#' Sinuosity of a trajectory
#'
#' Calculates the sinuosity of a (constant step length) trajectory as defined by
#' Bovet & Benhamou (1988), eqn 2, which is: \eqn{S = 1.18\sigma / \sqrt q} where
#' \eqn{\sigma} is the standard deviation of the step turning angles and \eqn{q}
#' is the mean step length. A corrected sinuosity index is available as the
#' function \code{\link{TrajSinuosity2}} which handles a wider range of
#' variations in step angles.
#'
#' If your trajectory does not have a constant step length, it should be
#' _rediscretized_ by calling \code{\link{TrajRediscretize}} before calling this
#' function.
#'
#' @param trj Trajectory to calculate sinuosity of.
#' @param compass.direction if not \code{NULL}, turning angles are calculated
#'   for a directed walk, assuming the specified compass direction (in radians).
#'   Otherwise, a random walk is assumed.
#' @return The sinuosity of \code{trj}.
#'
#' @seealso \code{\link{TrajSinuosity2}} for a corrected version of sinuosity,
#'   \code{\link{TrajAngles}} for the turning angles in a trajectory,
#'   \code{\link{TrajStepLengths}} for the step lengths, and
#'   \code{\link{TrajRediscretize}} for resampling to a constant step length.
#'
#' @references
#'
#' Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals' movements
#' using a correlated random walk model. Journal of Theoretical Biology, 131(4),
#' 419-433. doi:10.1016/S0022-5193(88)80038-9
#'
#' @export
TrajSinuosity <- function(trj, compass.direction = NULL) {
  segLen <- mean(TrajStepLengths(trj))
  1.18 * stats::sd(TrajAngles(trj, compass.direction = compass.direction), na.rm = TRUE) / sqrt(segLen)
}

#' Sinuosity of a trajectory
#'
#' Calculates the sinuosity of a trajectory as defined by Benhamou (2004),
#' equation 8. This is a corrected version of the sinuosity index defined in
#' Bovet & Benhamou (1988), which is suitable for a wider range of turning angle
#' distributions, and does not require a constant step length.
#'
#' This function implements equation 8 from Benhamou (2004):
#' \deqn{S = 2[p(\frac{1 + c}{1 - c} + b^2)]^{-0.5}}
#' where \eqn{p} is the mean step length, \eqn{c} is the mean
#' cosine of turning angles, and \eqn{b} is the coefficient of variation of
#' the step length.
#'
#' @param trj A Trajectory object.
#' @param compass.direction if not \code{NULL}, turning angles are calculated
#'   for a directed walk, assuming the specified compass direction (in radians).
#'   Otherwise, a random walk is assumed.
#'
#' @seealso \code{\link{TrajSinuosity}} for the uncorrected sinuosity index.
#'
#' @references
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's
#' path. Journal of Theoretical Biology, 229(2), 209-220.
#' doi:10.1016/j.jtbi.2004.03.016
#'
#' @export
TrajSinuosity2 <- function(trj, compass.direction = NULL) {
  stepLengths <- TrajStepLengths(trj)
  p <- mean(stepLengths)
  # Coefficient of variation of step length
  b <- stats::sd(stepLengths) / p
  # Mean cosine of turning angles = length of mean vector of angles
  c <- mean(cos(TrajAngles(trj, compass.direction = compass.direction)), na.rm = TRUE)
  # Same calculation but a bit slower
  #c <- Mod(TrajMeanVectorOfTurningAngles(trj, compass.direction = compass.direction))

  2 / sqrt(p * (((1 + c) / (1 - c)) + b^2))
}


# E-MAX ########################################################################################

#' Trajectory straightness index, E-max
#'
#' Emax, the maximum expected displacement, is a single-valued measure of
#' straightness defined by (Cheung, Zhang, Stricker, & Srinivasan, 2007). Emax-a
#' is a dimensionless, scale-independent measure of the maximum possible
#' expected displacement. Emax-b is \code{Emax-a * mean step length}, and gives
#' the maximum possible expected displacement in spatial units. Values closer to
#' 0 are more sinuous, while larger values (approaching infinity) are
#' straighter.
#'
#' @param trj Trajectory to be analysed.
#' @param eMaxB If TRUE, calculates and returns Emax-b, otherwise returns
#'   Emax-a.
#' @param compass.direction if not \code{NULL}, turning angles are calculated
#'   for a directed walk, assuming the specified compass direction (in radians).
#'   Otherwise, a random walk is assumed.
#' @return Emax (-a or -b) for \code{trj}.
#'
#' @references Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007).
#'   Animal navigation: the difficulty of moving in a straight line. Biological
#'   Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' @export
TrajEmax <- function(trj, eMaxB = FALSE, compass.direction = NULL) {

  # E(cos(angles)) = mean(cos(angles))
  b <- mean(cos(TrajAngles(trj, compass.direction = compass.direction)), na.rm = TRUE)

  # If it's E max b, multiply by mean step length
  f <- ifelse(eMaxB, mean(TrajStepLengths(trj)), 1)

  f * b / (1 - b)
}
