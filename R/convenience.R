# Convenience functions

# Private functions ####

.isBlank <- function(str) is.null(str) || is.na(str) || !is.character(str) || grepl("^\\s*$", str)

.locateFiles <- function(fileNames, rootDir) {
  sapply(fileNames, function(fn) {
    file <- ""
    # Ignore empty/NULL file names
    if (!.isBlank(fn)) {
      file <- list.files(rootDir, pattern = fn, full.names=TRUE, recursive=TRUE)
      if (length(file) == 0) {
        stop(sprintf("Unable to locate trajectory file %s (within folder %s)", fn, rootDir))
      }
      if (length(file) > 1) {
        stop(sprintf("There are %d files named '%s' within folder %s, file names must be unique", length(file), fn, rootDir))
      }
    }
    file
  })
}

# Reads a file and returns a list of trajectory coordinates
.readAndCheckCoords <- function(fileName, csvReadFn, ...) {
  coordList <- csvReadFn(fileName, ...)
  # Note a data.frame is a type of list, so we can't use is.list() here
  if (!inherits(coordList, "list"))
    coordList <- list(coordList)
  for (i in seq(1, length.out = length(coordList))) {
    coords <- coordList[[i]]
    # Check that there are at least 2 columns
    if (ncol(coords) < 2) {
      stop(sprintf("Invalid trajectory CSV file %s, contains %d column(s) but requires at least 2", fileName, ncol(coords)))
    }
  }
  coordList
}

# Public functions ####

#' Construct multiple trajectories
#'
#' Reads multiple trajectories from files, performs some basic sanity checks on
#' them, and optionally smooths and scales them. Attempts to collect and report
#' errors for multiple trajectories in a single call.
#'
#' If \code{rootDir} is not null, it should be the name of a directory which is
#' searched for the files in \code{fileNames}. The found files are then used as
#' the list of files to be read in. This may be useful when the names of the
#' files are known, but their exact location within a directory structure is not
#' known.
#'
#' For each file name in \code{fileNames}, reads the file by calling
#' \code{csvReadFn} to obtain a set of coordinates and optionally times. A
#' Trajectory is then constructed by passing the coordinates to
#' \code{\link{TrajFromCoords}}, passing in the appropriate \code{fps} value,
#' and x, y and time column names/indices from \code{csvStruct}. If \code{scale}
#' is not \code{NULL}, the trajectory is scaled by calling
#' \code{\link{TrajScale}}. If \code{smoothP} and \code{smoothN} are not
#' \code{NULL}, the trajectory is smoothed by calling
#' \code{\link{TrajSmoothSG}}.
#'
#' @param fileNames Vector of the names of CSV files containing trajectory
#'   coordinates. All of the files must have the same columns. All file names
#'   must be unique. If \code{rootDir} is not \code{NULL}, then the file names
#'   are treated as regular expressions.
#' @param fps Vector of frames-per-second values corresponding to the
#'   trajectories in \code{fileNames}. If length is 1, it is repeated to
#'   length(fileNames).
#' @param scale Vector of scale values corresponding to the trajectories in
#'   \code{fileNames}. May be specified as character expressions (e.g. "1 /
#'   1200") rather than numeric values. If NULL, the trajectories will not be
#'   scaled. If length is 1, it is repeated to length(fileNames).
#' @param spatialUnits Abbreviated name of spatial coordinate units after
#'   scaling, e.g. "m".
#' @param timeUnits Abbreviated name of temporal units, e.g. "s".
#' @param csvStruct A list which identifies the columns in each CSV file which
#'   contain x-, y-, and optionally time-values.
#' @param smoothP Filter order to be used for Savitzky-Golay smoothing (see
#'   \code{\link{TrajSmoothSG}}). If \code{NA}, no smoothing is performed.
#' @param smoothN Filter length to be used for Savitzky-Golay smoothing (must be
#'   odd, see \code{\link{TrajSmoothSG}}). If \code{NA}, no smoothing is
#'   performed.
#' @param translateToOrigin If TRUE, each trajectory is translated so that its
#'   starting point is at (0, 0).
#' @param rootDir Optional name of a top level directory which contains the CSV
#'   files. If \code{rootDir} is not NULL, the CSV files may be located anywhere
#'   within \code{rootDir} or its sub-directories.
#' @param csvReadFn Function used to read the CSV files. Required to accept
#'   arguments \code{filename, ...}, and return a data frame of coordinates, or
#'   a list of multiple data frames (see \code{\link[utils]{read.csv}},
#'   \code{\link[utils]{read.csv2}}). The default function calls
#'   \code{\link[utils]{read.csv}} with argument \code{stringsAsFactors =
#'   FALSE}.
#' @param ... Additional arguments passed to \code{csvReadFn}.
#'
#' @return A list of trajectories.
#'
#' @seealso \code{\link[utils]{read.csv}}, \code{\link{TrajFromCoords}},
#'   \code{\link{TrajScale}}, \code{\link{TrajSmoothSG}},
#'   \code{\link{TrajTranslate}}
#'
#' @examples
#' \dontrun{
#' # Names of CSV files containing trajectory coordinates
#' fileNames <- c('xy001.csv', 'xy003.csv', 'xy004.csv')
#' # The files are all located somewhere under this directory
#' rootDir <- '.'
#' # Scale for these files is 1 / pixels per metre
#' scale <- c('1/1200', '1/1350', '1/1300')
#' # Files have columns y, x
#' csvStruct <- list(x = 2, y = 1)
#' # Apply default smoothing, and the files are formatted as conventional CSV,
#' # so there's no need to specify csvReadFn
#' trjs <- TrajsBuild(fileNames, fps = 50, scale = scale, units = "m",
#'                    csvStruct = csvStruct, rootDir = rootDir)
#' }
#'
#' @export
TrajsBuild <- function(fileNames, fps = NULL, scale = NULL,
                       spatialUnits = NULL, timeUnits = NULL,
                       csvStruct = list(x = 1, y = 2, time = NULL),
                       smoothP = 3, smoothN = 41,
                       translateToOrigin = FALSE,
                       rootDir = NULL,
                       csvReadFn = function (filename, ...) utils::read.csv(filename, stringsAsFactors = FALSE, ...),
                       ...) {
  # I hate factors!
  fileNames <- as.character(fileNames)
  # Check that file names are unique
  if (any(table(fileNames) > 1)) {
    stop(sprintf("List of file names contains duplicates: %s\n",
                 paste(sprintf("'%s'", names(which(table(fileNames) > 1))), collapse = ", ")))
  }

  # Fail with a (hopefully) helpful message on empty file name
  blanks <- which(sapply(fileNames, .isBlank))
  if (length(blanks) > 0) {
    stop(sprintf("Trajectory input file name is blank or NULL (index %s)", paste(blanks, collapse = ", ")))
  }

  # Maybe locate files within rootDir
  if (!is.null(rootDir)) {
    fileNames <- .locateFiles(fileNames, rootDir)
  }

  # Extend fps and scale if required
  if (!is.null(fps) && length(fps) == 1)
    fps <- rep(fps, length(fileNames))
  if (!is.null(scale) && length(scale) == 1)
    scale <- rep(scale, length(fileNames))

  result <- list()

  # Build up a list of errors so they can all be reported at once,
  # rather than annoying report one, fix one, report next one...
  errors <- character(0)

  # For each file...
  for (i in 1:length(fileNames)) {

    fn <- as.character(fileNames[i])

    # Read the trajectory coordinates from the file
    coordList <- .readAndCheckCoords(fn, csvReadFn, ...)
    for (j in seq(1, length.out = length(coordList))) {
      # Convert to a trajectory
      trj <- tryCatch(
        TrajFromCoords(coordList[[j]], fps = fps[i],
                       xCol = csvStruct$x, yCol = csvStruct$y, timeCol = csvStruct$time,
                       spatialUnits = spatialUnits, timeUnits = timeUnits),
        error = function (e) {
          errors <<- c(errors, sprintf("Trajectory file %s (index %d): %s", fn, i, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(trj)) {
        # Scale
        if (is.character(scale[i]) || is.numeric(scale[i])) {
          # Allow scale to be specified as a character expression such as "1 / 250"
          sc <- ifelse(is.character(scale[i]), eval(parse(text=scale[i])), scale[i])
          trj <- TrajScale(trj, sc, spatialUnits)
        }

        # Smooth
        if (is.numeric(smoothP) && is.numeric(smoothN)) {
          trj <- TrajSmoothSG(trj, p = smoothP, n = smoothN)
        }

        # Translate to origin
        if (translateToOrigin) {
          trj <- TrajTranslate(trj, -trj$x[1], -trj$y[1])
        }
      }

      result[[length(result) + 1]] <- trj
    }
  }

  # Check for errors
  if (length(errors) > 0)
    stop(paste(errors, collapse = "\n"))

  result
}

#' Merge trajectory characteristics
#'
#' Builds a data frame by combining rows of statistical values for multiple
#' trajectories. The statistics for each trajectory are defined by the caller in
#' a user defined function - see the example for one way to achieve this.
#'
#' @section Note: Any NULL valued statistics are converted to NAs.
#'
#' @param trjs List of trajectories to be characterised.
#' @param statsFn Function to calculate statistics of interest for a single
#'   trajectory.
#' @param progressBar Displays an optional progressbar, which may be helpful if
#'   processing is very slow. The progressbar is displayed by printing to the
#'   console, by using \code{winProgressBar} or
#'   \code{\link[tcltk]{tkProgressBar}}, if \code{progressBar} is \code{"text"},
#'   \code{"win"} or \code{"tk"} respectively. The default is no progressbar
#'   (value \code{"none"}). The \code{"win"} progressbar is only available on
#'   Windows.
#' @param check.names Passed to \code{as.data.frame}. Set to \code{FALSE} if
#'   statistic names are not syntactically valid variable names.
#' @param ... Additional arguments passed to \code{statsFn}.
#'
#' @examples
#' \dontrun{
#'
#' # Define a function which calculates some statistics
#' # of interest for a single trajectory
#' characteriseTrajectory <- function(trj) {
#'   # Measures of speed
#'   derivs <- TrajDerivatives(trj)
#'   mean_speed <- mean(derivs$speed)
#'   sd_speed <- sd(derivs$speed)
#'
#'   # Resample to constant step length.
#'   # Step length must be appropriate for the trajectory
#'   resampled <- TrajRediscretize(trj, 2)
#'
#'   # Measures of straightness
#'   sinuosity <- TrajSinuosity2(resampled)
#'   Emax <- TrajEmax(resampled)
#'
#'   # Periodicity
#'   resampled <- TrajRediscretize(trj, .001)
#'   corr <- TrajDirectionAutocorrelations(resampled, round(nrow(resampled) / 4))
#'   first_min <- TrajDAFindFirstMinimum(corr)
#'
#'   # Return a list with all of the statistics for this trajectory
#'   list(mean_speed = mean_speed,
#'        sd_speed = sd_speed,
#'        sinuosity = sinuosity,
#'        Emax = Emax,
#'        first_min_deltaS = first_min[1],
#'        first_min_C = first_min[2])
#' }
#'
#' trjs <- TrajsBuild(filenames)
#' stats <- TrajsMergeStats(trjs, characteriseTrajectory)
#' }
#'
#' @export
TrajsMergeStats <- function(trjs, statsFn, progressBar = c("none", "text", "win", "tk"), check.names = TRUE, ...) {
  # TODO rename this function? perhaps TrajsCombineIndices (or TrajsRBindStats) would be more meaningful

  # Setup the progress bar
  progressBar <- match.arg(progressBar)
  pb <- switch(progressBar,
               none = function(close){},
               text = ElapsedTimeProgressBarFn(length(trjs), buildTxtReportFn()),
               win = ElapsedTimeProgressBarFn(length(trjs), buildWinReportFn("TrajsMergeStats")),
               tk = ElapsedTimeProgressBarFn(length(trjs), buildTkReportFn("TrajsMergeStats"))
  )
  tryCatch({
    result <- data.frame()
    nc <- NA

    for (trj in trjs) {
      pb()
      row <- statsFn(trj, ...)
      # Replace NULL with NA because NULL values are removed by rbind, resulting
      # in: invalid list argument: all variables should have the same length
      row[sapply(row, is.null)] <- NA

      if (is.na(nc))
        nc <- length(row)
      else if (nc != length(row))
        stop(sprintf("Statistics for trajectory %d contains %d values, but expected %d", nrow(result), length(row), nc))
      dfrow <- as.data.frame(row, stringsAsFactors = FALSE, check.names = check.names)
      result <- rbind(result, dfrow)
    }
    result
  },
  finally = pb(close = TRUE)
  )
}

#' Step lengths of multiple trajectories
#'
#' Returns the lengths of all of the steps in a list of trajectories
#'
#' @param trjs A list of \code{Trajectory} objects.
#' @return A numeric vector containing the lengths of every step in every trajectory.
#'
#' @examples
#' \dontrun{
#' trjs <- TrajsBuild(fileNames, scale = scale, units = "m")
#' # Print a summary about the step sizes across all trajectories
#' print(summary(TrajsStepLengths(trjs)))
#' }
#'
#' @export
TrajsStepLengths <- function(trjs) {
  # First displacement is 0, not a real displacement, so don't return it
  unlist(lapply(trjs, TrajStepLengths))
}

#' Replace NAs in a data frame
#'
#' Replaces NAs in a single column of a data frame with an imputed uninformative
#' numeric replacement value, so that a principal component analysis can be
#' applied without discarding data. Optionally adds a new "flag" column which
#' contains \code{1} for each row which originally contained NA, otherwise
#' \code{0}.
#'
#' @param df Data frame to be adjusted.
#' @param column Name or index of the column to be adjusted.
#' @param replacementValue Numeric value to use instead of NA.
#' @param flagColumn If not NULL, specifies the name of a new column to be added
#'   to the data frame, with value 0 for non-NA rows, 1 for NA rows. The column
#'   is added regardless of whether there are any NAs in the data.
#' @return A copy of \code{df} with NAs replaced in \code{column}.
#'
#' @seealso \code{\link[stats]{prcomp}}
#'
#' @examples
#' df <- data.frame(x = c(1, 2, 3), y = c(NA, 5, 6), z = c(NA, NA, 9))
#' # Eliminate NAs in y, add a flag column, ignore other NAs
#' df <- TrajsStatsReplaceNAs(df, "y", flagColumn = "y.was.NA")
#' print(df)
#'
#' @export
TrajsStatsReplaceNAs <- function(df, column, replacementValue = mean(df[,column], na.rm = TRUE), flagColumn = NULL) {
  # Are there any NAs in the column?
  col <- df[,column]
  nas <- which(is.na(col))
  df[nas,column] <- replacementValue
  if (!is.null(flagColumn)) {
    df[,flagColumn] <- as.numeric(is.na(col))
  }
  df
}

# Time conversion ####

#' Converts a delimited time string to a numeric value
#'
#' Time values may be imported in a format which is not immediately usable by
#' `trajr`. This function converts times that are specified as a number of
#' delimited fields to a single numeric value. The default parameter values
#' handle a value with 4 colon-separated values, which are hours, minutes,
#' seconds and milliseconds, eg: "0:01:04:108" represents 1 minute, 4 seconds
#' and 108 milliseconds, or 64.108 seconds.
#'
#' Note that the base R strptime can be used to convert time values in more
#' complex date/time formats, but it does not handle millisecond fields.
#'
#' @param time A character string containing the time value to be converted.
#' @param sep Field separator.
#' @param factors Vector of numeric factors to be applied to each field, in the
#'   order they occur within `time`. The default assumes 4 fields containing
#'   numeric hours, minutes, seconds and milliseconds.
#' @return `time` converted to a numeric value.
#'
#' @seealso \code{\link[base]{strptime}}
#'
#' @examples
#' time <- c("0:00:00:029", "0:01:00:216", "0:02:01:062", "1:00:02:195", "1:06:03:949", "1:42:04:087")
#' seconds <- TrajConvertTime(time)
#'
#' @export
TrajConvertTime <- function(time, sep = ":", factors = c(60 * 60, 60, 1, .001)) {
  # Split time on separator, and use matrix (inner) multiplication to
  # convert each component to seconds then sum all components
  c(as.matrix(utils::read.table(text = time, sep = sep)) %*% factors)
}
