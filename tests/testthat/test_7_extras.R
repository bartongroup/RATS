#==============================================================================
#==============================================================================
context("Non-core functionality")

#==============================================================================
test_that("GTF annotations are parsed correctly", {
  expected_values <- data.table(target_id = c("GENE2.1","GENE2.2","GENE3.1","GENE3.2","GENE4.1","GENE5.1","GENE5.2","GENE5.3","GENE7.1","TRANSCRIPT8.1","GENE9.1","TRANSCRIPT9.2"),
                                parent_id = c("GENE2","GENE2","GENE3","GENE3","GENE4","GENE5","GENE5","GENE5","GENE7","GENE8","GENE9","GENE9"))
  
  expect_silent(result <- annot2ids("./test.gtf"))
  expect_identical(dim(result), dim(expected_values))
  expect_identical(names(result), names(expected_values))
  expect_false("ORPHAN1" %in% result$parent_id)
  expect_false("ORPHAN6.1" %in% result$target_id)
  expect_identical(result, expected_values)
  
  names(expected_values) <- c("foo", "bar")
  expect_silent(result <- annot2ids("./test.gtf", transc_header = "foo", gene_header = "bar"))
  expect_identical(result, expected_values)
})
