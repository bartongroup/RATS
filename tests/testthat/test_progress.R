context("Progress reporting")

#==============================================================================
test_that("The progress text output is not correct", {
  
  # Set up progress updates
  progress_steps <- data.frame(c(10,20,30),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=TRUE)
  
  expect_output(progress <- update_progress(progress), "Update1.....10%")
  expect_output(progress <- update_progress(progress), "Update2.....20%")
  expect_output(progress <- update_progress(progress), "Update3.....30%")
})

#==============================================================================
test_that("Setting 'on' to false does not silence progress text output", {
  
  # Set up progress bar
  progress_steps <- data.frame(c(10,20,30),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=FALSE)
  
  expect_silent(progress <- update_progress(progress))
  
})

#==============================================================================
test_that("Setting 'on' to false does not silence progress bar output", {
  
  # Set up progress bar
  progress_steps <- data.frame(c(10,20,30),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- BarProgressUpdate(steps=progress_steps, on=FALSE)
  
  expect_silent(progress <- update_progress(progress))
  
})
