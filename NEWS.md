# Changes to the `trajr` package

## trajr 1.5.1.9000

* Add more 3d functionality: `Traj3DVelocity`, `Traj3DAcceleration` and `Traj3DSpeed`.
* Fix bug in simulation vignette.

## trajr 1.5.1

* Minor documentation changes: mention 3D in package overview and use new roxygen method to generate package documentation (https://github.com/r-lib/roxygen2/issues/1491).

## trajr 1.5.0

* Added some basic 3D trajectory functionality; see `Traj3DFromCoords` and all functions prefixed with `Traj3D`.
* `TrajAngles` now returns `NA` for the angles before and after zero-length segments. Thanks to Valentin Baumann for identifying this bug. This change could potential cause problems with existing code, but such a problem would indicate a problem with the old analysis. The bug occurred because in base R, `Arg(0) == 0`, which is not correct. Various functions that take means of TrajAngles now remove NA values.
* Added function `TrajFromTrjPoints` to simplify operations such as removing zero-length segments from a trajectory.
* Fix documentation of `TrajDirectionalChange` and `TrajAngles`. The documentation incorrectly stated that the number of angles returned was number of points - 1. It now states that the number of angles is number of points - 2 (which is number of segments - 1).
* Fix bug in `TrajFromCoords` that incorrectly treated NA values in columns other than coordinates or time as an error.
* Fix bug in `TrajScale`; `yScale` was being ignored.
* Add some more named parameters to `plot.Trajectory` to avoid passing them to `graphics::lines`. 
* Add arg `simConstantSpeed` to `TrajRediscretize` and `Traj3DRediscretize`.
* Add arg `check.names` to `TrajsMergeStats`.
* Updated github continuous integration actions (using usethis::use_github_action("test-coverage") and usethis::use_github_action_check_release()). This change does not affect users

## trajr 1.4.0

* Allow the creation of trajectories with 0 points in `TrajFromCoords`.
* Fixed bug in TrajSpeedIntervals: no intervals were being returned if the entire trajectory qualified as an interval.
* Fix bug in plot.TrajSpeedIntervals: time (x-axis) didn't handle trajectories with a non-zero starting time.
* Added functions to assist in simulating trajectories that are bounded or vary spatially in some way. The new functions are `TrajInPolygon`, `TrajSplit`, `TrajMerge` and `TrajSplitAtFirstCrossing`.
* Added functions to calculate velocity and acceleration: `TrajVelocity` and `TrajAcceleration`. The documentation for `TrajDerivatives` has been updated to state that the `acceleration` returned is not really acceleration, rather it is change in speed over time.
* Added argument `diff` to `TrajSpeedIntervals` to control how speed is calculated. Default is "backward" so the function is backwards-compatible, although "central" is usually a better choice.

## trajr 1.3.0

* Detect and report some invalid parameter values.
* Altered handling of parameter `readcsvFn` to `TrajsBuild` to make it 
  possible to use `readr::read_csv` without a wrapper function.
* Added function `TrajResampleTime` to resample a trajectory to fixed step times.
* Added parameters `start.pt.pch` and `start.pt.col` to plotting functions.
* Added parameter `dt` to `TrajTranslate`
* Fix vertical extents of rectangles in `plot.TrajSpeedIntervals` to handle non-default ylim values.
* Added optional progressbar to `TrajsMergeStats`.
* TrajsMergeStats now passes the arguments `stringsAsFactors = FALSE` to rbind. This prevents incorrect 
  behaviour and the warning "invalid factor level, NA generated" if one or more of your statistics are characters.
* Enhanced `TrajRotate` to allow absolute rotation and arbitrary origin of rotation.

## trajr 1.2.0

* Added start.pt.cex parameter to function `lines.Trajectory`.
* Added function `TrajConvertTime`.

## trajr 1.1.0

* Added correct citation.
* Fixed: `plot.TrajSpeedIntervals` was not passing additional arguments (`...`) to `plot`.
* Added: functions `TrajDuration`, `TrajMeanVelocity`, `TrajTranslate`.
* Added `translateToOrigin` parameter to function `TrajsBuild`.
