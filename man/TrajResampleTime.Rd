% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resample.R
\name{TrajResampleTime}
\alias{TrajResampleTime}
\title{Resample a trajectory to a constant time interval.}
\usage{
TrajResampleTime(trj, stepTime, newFps = NULL)
}
\arguments{
\item{trj}{The trajectory to be resampled.}

\item{stepTime}{The resampled trajectory step time. Each step in the new
trajectory will have this duration.}

\item{newFps}{Value to be stored as the FPS value in the new trajectory (see
\code{\link{TrajGetFPS}}). It is not otherwise used by this function.}
}
\value{
A new trajectory with a constant time interval for each step. Points
  in the new trajectory are calculated by linearly interpolating along
  \code{trj}.
}
\description{
Constructs a new trajectory by resampling the input trajectory to a fixed
time interval. Points are linearly interpolated along the trajectory. Spatial
and time units are preserved.
}
\details{
Resampling with the original step time (i.e., 1/fps) can be used to fill in
missing frames in a trajectory.
}
\examples{
# Simulate a trajectory with steps every 5 hours
set.seed(46)
trj <- TrajGenerate(10, stepLength = 5, fps = 1/5, timeUnits = "hours", linearErrorSd = .8)

# Resample to 1 hour steps
resampled <- TrajResampleTime(trj, 1)

par(mar = c(5, 4, .5, .5))
plot(trj, lwd = 2)
points(trj, pch = 16)
points(resampled, col = "red", draw.start.pt = FALSE)

}
