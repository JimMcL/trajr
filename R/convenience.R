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
#' them, and optionally smooths and scales them.
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
#' # Names of CSV files containing trajectory coordinates
#' fileNames <- c('xy001.csv', 'xy003.csv', 'xy004.csv')
#' # The files are all located under this directory
#' rootDir <- 'c:/Jim/data'
#' # Scale is 1 / pixels per metre
#' scale <- c('1/1200', '1/1350', '1/1300')
#' # Files have columns y, x
#' csvStruct <- list(x = 2, y = 1)
#' # Apply default smoothing, and the files are conventional CSV, so no need to specify csvReadFn
#' trjs <- TrajsBuild(fileNames, fps = rep(50, length(fileNames)),
#'                    scale = scale, units = "m",
#'                    csvStruct = csvStruct, rootDir = rootDir)
TrajsBuild <- function(fileNames, fps, scale = NULL, units = NULL, csvStruct = list(x = 1, y = 2, time = NULL),
                       smoothP = 3, smoothN = 41,
                       rootDir = NULL,
                       csvReadFn = read.csv) {
  # Check that file names are unique
  if (any(table(fileNames) > 1)) {
    stop(sprintf("List of file names contains duplicates: %s\n", paste(names(which(table(fileNames) > 1)), collapse = ", ")))
  }
  # Maybe locate files within rootDir
  if (!is.null(rootDir)) {
    fileNames <- .locateFiles(fileNames, rootDir)
  }

  result <- list(length(fileNames))

  # For each file...
  for (i in 1:length(fileNames)) {
    # Read the trajectory coordinates from the file
    coords <- .readAndCheckCoords(fileNames[i], csvReadFn)
    # Convert to a trajectory
    trj <- withCallingHandlers(
      TrajFromCoords(coords, fps = fps[i], xCol = csvStruct$x, yCol = csvStruct$y, timeCol = csvStruct$time),
      error = function (e) stop(sprintf("Trajectory file %s:\n%s", fileNames[i], e))
    )

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

    result[[i]] <- trj
  }

  result
}

#' Statistically characterise a list of trajectories
#'
#' Builds a data frame containing statistical values for multiple trajectories.
#' Each row contains the statistics for a single trajectory.
#' 
#' @param trjs List of trajectories to be characterised.
TrajsCharacterise <- function(trjs, statsFn) {
  
}
