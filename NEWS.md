# Changes to the `trajr` package

## trajr 1.2.0.9001

Current development version, intended to be released as version 1.3.0.

* Altered handling of parameter `readcsvFn` to `TrajsBuild` to make it 
  possible to use `readr::read_csv` without a wrapper function.
* Added function `TrajResampleTime` to resample a trajectory to fixed step times.

## trajr 1.2.0

* Added start.pt.cex parameter to function `lines.Trajectory`.
* Added function `TrajConvertTime`.

## trajr 1.1.0

* Added correct citation.
* Fixed: `plot.TrajSpeedIntervals` was not passing additional arguments (`...`) to `plot`.
* Added: functions `TrajDuration`, `TrajMeanVelocity`, `TrajTranslate`.
* Added `translateToOrigin` parameter to function `TrajsBuild`.
