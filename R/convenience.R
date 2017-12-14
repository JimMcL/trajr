# Convenience functions

# Private functions ####

.locateFiles <- function(fileNames, rootDir) {
  sapply(fileNames, function(fn) {
    file <- list.files(rootDir, fn, full.names=TRUE, recursive=TRUE)
    if (length(file) == 0) {
      stop(sprintf("Unable to locate trajectory file %s (within folder %s)", fn, rootDir))
    }
    if (length(file) > 1) {
      stop(sprintf("There are %d files named %s within folder %s, file names must be unique", length(file), fn, rootDir))
    }
    file
  })
}

.readAndCheckCoords <- function(fileName, csvReadFn) {
  coords <- csvReadFn(fileName, stringsAsFactors = FALSE)
  # Check that there are at least 2 columns
  if (ncol(coords) < 2) {
    stop(sprintf("Invalid trajectory CSV file %s, contains %d column but requires at least 2", fileName, ncol(coords)))
  }
  coords
}

# Public functions ####

#' Construct multiple trajectories
#'
#' Reads multiple trajectories from files, performs some basic sanity checks on
#' them, and optionally smooths and scales them. Attempts to collect and report
#' errors for multiple trajectories in a single call.
#'
#' For each file name in \code{fileNames}, searches through the folder
#' \code{rootDir} (unless it's \code{NULL}) to find the file, then reads the
#' file by calling \code{csvReadFn} to obtain a set of coordinates and
#' optionally times. A Trajectory is constructed by passing the coordinates to
#' \code{\link{TrajFromCoords}}, passing in the appropriate \code{fps} value,
#' and x, y and time column names/indices from \code{csvStruct}. If \code{scale}
#' is not \code{NULL}, the trajectory is then scaled by calling
#' \code{\link{TrajScale}}. If \code{smoothP} and \code{smoothN} are not
#' \code{NULL}, the trajectory is smoothed by calling
#' \code{\link{TrajSmoothSG}}.
#'
#' @param fileNames Vector of the names of CSV files containing trajectory
#'   coordinates. All of the files must have the same columns. All file names
#'   must be unique.
#' @param fps Vector of frames-per-second values corresponding to the
#'   trajectories in \code{fileNames}.
#' @param scale Vector of scale values corresponding to the trajectories in
#'   \code{fileNames}. May be specified as character expressions (e.g. "1 /
#'   1200") rather than numeric values. If NULL, the trajectories will not be
#'   scaled.
#' @param units Name of spatial coordinate units after scaling, e.g. "m".
#' @param csvStruct A list which identifies the columns in each CSV file which
#'   contain x-, y-, and optionally time-values.
#' @param smoothP Filter order to be used for Savitzky-Golay smoothing (see
#'   \code{\link{TrajSmoothSG}})
#' @param smoothN Filter length to be used for Savitzky-Golay smoothing (must be
#'   odd, see \code{\link{TrajSmoothSG}})
#' @param rootDir Optional name of a top level directory which contains the CSV
#'   files. If \code{rootDir} is not NULL, the CSV files may be located anywhere
#'   within \code{rootDir} or its sub-directories, .
#' @param csvReadFn Function used to read the CSV files (see
#'   \code{\link[utils]{read.csv}}, \code{\link[utils]{read.csv2}}).
#'
#' @return A list of trajectories.
#'
#' @seealso \code{\link[utils]{read.csv}}, \code{\link{TrajFromCoords}},
#'   \code{\link{TrajScale}}, \code{\link{TrajSmoothSG}}
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
#' trjs <- TrajsBuild(fileNames, fps = rep(50, length(fileNames)),
#'                    scale = scale, units = "m",
#'                    csvStruct = csvStruct, rootDir = rootDir)
#' }
#'
#' @export
TrajsBuild <- function(fileNames, fps, scale = NULL, units = NULL, csvStruct = list(x = 1, y = 2, time = NULL),
                       smoothP = 3, smoothN = 41,
                       rootDir = NULL,
                       csvReadFn = utils::read.csv) {
  # Check that file names are unique
  if (any(table(fileNames) > 1)) {
    stop(sprintf("List of file names contains duplicates: %s\n", paste(names(which(table(fileNames) > 1)), collapse = ", ")))
  }
  # Maybe locate files within rootDir
  if (!is.null(rootDir)) {
    fileNames <- .locateFiles(fileNames, rootDir)
  }

  result <- list(length(fileNames))

  # Build up a list of errors so they can all be reported at once,
  # rather than annoying report one, fix one, report next one...
  errors <- character(0)

  # For each file...
  for (i in 1:length(fileNames)) {
    # Read the trajectory coordinates from the file
    coords <- .readAndCheckCoords(fileNames[i], csvReadFn)
    # Convert to a trajectory
    trj <- tryCatch(
      TrajFromCoords(coords, fps = fps[i], xCol = csvStruct$x, yCol = csvStruct$y, timeCol = csvStruct$time),
      error = function (e) {
        errors <<- c(errors, sprintf("Trajectory file %s (index %d): %s", fileNames[i], i, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(trj)) {
      # Scale
      if (!is.null(scale[i])) {
        # Allow scale to be specified as a character expression such as "1 / 250"
        sc <- ifelse(is.character(scale[i]), eval(parse(text=scale[i])), scale[i])
        trj <- TrajScale(trj, sc, units)
      }

      # Smooth
      if (is.numeric(smoothP) && is.numeric(smoothN)) {
        trj <- TrajSmoothSG(trj, 3, 101)
      }
    }

    result[[i]] <- trj
  }

  # Check for errors
  if (length(errors) > 0)
    stop(paste(errors, collapse = "\n"))

  result
}

#' Merge trajectory characteristics
#'
#' Builds a data frame by combining rows of statistical values for multiple trajectories.
#'
#' @section Note:
#' Any NULL valued statistics are converted to NAs.
#'
#' @param trjs List of trajectories to be characterised.
#' @param statsFn Function to calculate statistics of interest for a single trajectory.
#'
#' @export
TrajsMergeStats <- function(trjs, statsFn) {
  result <- data.frame()
  nc <- NA
  rowNum <- 1
  for (trj in trjs) {
    row <- statsFn(trj)
    # Replace NULL with NA because NULL values are removed by rbind, resulting
    # in: invalid list argument: all variables should have the same length
    row[sapply(row, is.null)] <- NA

    if(is.na(nc))
      nc <- length(row)
    else if (nc != length(row))
      stop(sprintf("Statistics for trajectory %d contain %d values, but expected %d", rowNum, length(row), nc))
    result <- rbind(result, row)

    rowNum <- rowNum + 1
  }
  result
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
#' print(summary(TrajsStepLength(trjs)))
#' }

#'
#' @export
TrajsStepLengths <- function(trjs) {
  # First displacement is 0, not a real displacement, so don't return it
  unlist(lapply(trjs, TrajStepLengths))
}

#' Replace NAs in a data frame
#'
#' Replaces NAs in a single column of a data frame with an uninformative numeric
#' replacement value, so that a principal components analysis can be applied
#' without discarding data. Optionally adds a new "flag" column which contains
#' \code{1} for each row which originally contained NA, otherwise \code{0}.
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
