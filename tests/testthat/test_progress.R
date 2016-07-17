context("Progress reporting")

#==============================================================================
test_that("The progress text output is correct", {
  
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
test_that("The progress text output can have new steps inserted", {
  
  # Set up progress updates
  progress_steps <- data.frame(c(10,20,100),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=TRUE)
  
  expect_output(progress <- update_progress(progress), "Update1.....10%")
  expect_output(progress <- update_progress(progress), "Update2.....20%")
  
  # Set up inserted progress updates
  new_steps <- data.frame(c(40,60),
                          c("Inserted Update1","Inserted Update2"),
                          stringsAsFactors = FALSE)
  progress <- insert_steps(progress, new_steps)
  
  expect_output(progress <- update_progress(progress), "Inserted Update1.....40%")
  expect_output(progress <- update_progress(progress), "Inserted Update2.....60%")
  expect_output(progress <- update_progress(progress), "Update3.....100%")
})

#==============================================================================
test_that("Setting 'on' to false silences progress text output", {
  
  # Set up progress bar
  progress_steps <- data.frame(c(10,20,30),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- TxtProgressUpdate(steps=progress_steps, on=FALSE)
  
  expect_silent(progress <- update_progress(progress))
  
})

#==============================================================================
test_that("Setting 'on' to false silences progress bar output", {
  
  # Set up progress bar
  progress_steps <- data.frame(c(10,20,30),
                               c("Update1","Update2","Update3"),
                               stringsAsFactors = FALSE)
  progress <- BarProgressUpdate(steps=progress_steps, on=FALSE)
  
  expect_silent(progress <- update_progress(progress))
  
})
