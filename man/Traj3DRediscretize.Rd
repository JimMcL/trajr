% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3D.R
\name{Traj3DRediscretize}
\alias{Traj3DRediscretize}
\title{Resample a 3D trajectory to a constant step length}
\usage{
Traj3DRediscretize(trj3d, R, simConstantSpeed = FALSE)
}
\arguments{
\item{trj3d}{The 3-dimensional trajectory to be resampled.}

\item{R}{rediscretization step length, in the spatial units of \code{trj}.}

\item{simConstantSpeed}{If TRUE, speeds are interpolated along the new
trajectory so that average speed is approximately the same as that of
\code{trj3d}.}
}
\value{
A new 3-dimensional trajectory with a constant segment length which
  follows \code{trj3d}.
}
\description{
Constructs a new 3-dimensional trajectory by resampling the input
3-dimensional trajectory to a fixed step (or segment) length. By
default, timing of frames is lost, so speed and acceleration cannot
be calculated on a rediscretized trajectory. However, a constant
speed may be applied to the rediscretized trajectory
(\code{simConstantSpeed = TRUE}), in which case the returned
trajectory will have (almost) constant speed, with average speed
approximately equal to the average speed of \code{trj3d}.
}
\details{
Unfortunately this operation is slow for large trajectories.

Based on the appendix in Bovet and Benhamou, (1988), extended to 3
dimensions.
}
\references{
Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals'
  movements using a correlated random walk model. Journal of Theoretical
  Biology, 131(4), 419-433. doi:10.1016/S0022-5193(88)80038-9
}
\seealso{
\code{\link{Traj3DFromCoords}}, \code{\link{TrajRediscretize}}
}
