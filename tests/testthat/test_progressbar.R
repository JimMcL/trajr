context("progressbar tests")

test_that("progress bar", {
  # Just check that the progress bar - text output - doesn't crash. This check
  # exists mainly to unskew the coverage statistics because I can't test it
  # properly without lots of work

  expect_error(capture.output({
    n <- 20
    pb <- trajr:::ElapsedTimeProgressBarFn(n, trajr:::buildTxtReportFn("Progress"))
    for (i in 1:(n - 1)) {
      pb()
    }
    pb(close = TRUE)
    pb()
  }), NA)
})
