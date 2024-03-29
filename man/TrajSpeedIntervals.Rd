% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speed.R
\name{TrajSpeedIntervals}
\alias{TrajSpeedIntervals}
\title{Calculate speed time intervals}
\usage{
TrajSpeedIntervals(
  trj,
  fasterThan = NULL,
  slowerThan = NULL,
  interpolateTimes = TRUE,
  diff = c("backward", "central", "forward")
)
}
\arguments{
\item{trj}{Trajectory to be analysed.}

\item{fasterThan, slowerThan}{If not \code{NULL}, intervals will cover time
periods where speed exceeds/is lower than this value. At least one of
\code{fasterThan}, or \code{slowerThan} must be specified.}

\item{interpolateTimes}{If \code{TRUE}, times will be linearly interpolated
between frames.}

\item{diff}{Method used to calculate speed, see \code{\link{TrajVelocity}}
for details. The default is \code{"backward"} to maintain backwards
compatibility; in general, \code{"central"} provides a more accurate
estimate of velocity.}
}
\value{
A data frame of class "TrajSpeedIntervals", each row is an interval,
  columns are: \item{startFrame}{Indices of frames at the start of each
  interval.} \item{stopFrame}{Indices of frames at the end of each interval.}
  \item{startTime}{Time since start of trajectory at the start of each
  interval.} \item{stopTime}{Time since start of trajectory at the end of
  each interval} \item{duration}{Duration of each interval.}

  The data frame will also have non-standard attributes:
  \item{trajectory}{Value of the \code{trj} argument.}
  \item{slowerThan}{Value of the \code{slowerThan} argument.}
  \item{fasterThan}{Value of the \code{fasterThan} argument.}
  \item{speed}{Data frame with columns \code{speed} and \code{time}.}
  \item{derivs}{Value returned by calling \code{TrajDerivatives(trj)}.
  Provided for backwards-compatibility; use of the \code{speed} attribute
  is now preferred to the \code{derivs} attribute.}
}
\description{
Calculates and returns a list of time intervals during which speed is slower
and/or faster than specified values. Speed is calculated by taking the
modulus of velocity (\code{\link{TrajVelocity}}).
}
\examples{
# Plot speed, highlighting intervals where speed drops below 50 units/sec
set.seed(4)
trj <- TrajGenerate(200, random = TRUE)
smoothed <- TrajSmoothSG(trj, 3, 101)
intervals <- TrajSpeedIntervals(smoothed, diff = "central", slowerThan = 50, fasterThan = NULL)
plot(intervals)

# Report the duration of the longest period of low speed
cat(sprintf("Duration of the longest low-speed interval was \%g secs\n", max(intervals$duration)))

}
\seealso{
\code{\link{TrajVelocity}} for calculating trajectory velocity,
  \code{\link{plot.TrajSpeedIntervals}} for plotting speed over time with
  intervals highlighted.
}
