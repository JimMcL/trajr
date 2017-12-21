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
#' @return The directional change (DC) in degrees between every pair of
#'   consecutive points in the trajectory, i.e. the returned vector will have
#'   length \code{(nrow(trj) - 1)}.
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

    # Calculating this way is almost 1 order of magnitude faster than using the documented equation
  abs(.rad2deg(TrajAngles(trj, nFrames))) / diff(trj$displacementTime, 2 * nFrames)
}


# Sinuosity ########################################################################################

#' Sinuosity of a trajectory
#'
#' Calculates the sinuosity of a trajectory as defined by Bovet & Benhamou
#' (1988), which is: \eqn{S = 1.18\sigma / \sqrt q} where \eqn{\sigma} is the
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
#'   \code{\link{TrajStepLengths}} for the step lengths.
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
  1.18 * stats::sd(TrajAngles(trj, compass.direction = compass.direction)) / sqrt(segLen)
}


# E-MAX ########################################################################################

#' Trajectory straightness index, E-max
#'
#' Emax is a single-valued measure of straightness defined by (Cheung, Zhang,
#' Stricker, & Srinivasan, 2007). Emax-a is a dimensionless, scale-independent
#' measure of the maximum possible expected displacement. Emax-b is Emax-a *
#' mean step length, and gives the maximum possible expected displacement in
#' spatial units. Values closer to 0 are more sinuous, while larger values
#' (approaching infinity) are straighter.
#'
#' @param trj Trajectory to be analysed.
#' @param eMaxB If TRUE, calculates and returns Emax-b, otherwise returns Emax-a.
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
  b <- mean(cos(TrajAngles(trj, compass.direction = compass.direction)))

  # If it's E max b, multiply by mean step length
  f <- ifelse(eMaxB, mean(TrajStepLengths(trj)), 1)

  f * b / (1 - b)
}
