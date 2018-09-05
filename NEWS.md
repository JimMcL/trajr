# Changes to the `trajr` package

## trajr 1.2.1

* Altered handling of parameter `readcsvFn` to `TrajsBuild` to make it 
  possible to use `readr::read_csv` without a wrapper function.

## trajr 1.2.0

* Added start.pt.cex parameter to function `lines.Trajectory`.
* Added function `TrajConvertTime`.

## trajr 1.1.0

* Added correct citation.
* Fixed: `plot.TrajSpeedIntervals` was not passing additional arguments (`...`) to `plot`.
* Added: functions `TrajDuration`, `TrajMeanVelocity`, `TrajTranslate`.
* Added `translateToOrigin` parameter to function `TrajsBuild`.
