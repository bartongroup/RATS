library(rats)
context("DTU results")

#==============================================================================
test_that("The output structure is correct", {
  # All the expectations below should be general, regardless of the exact test data.
  dtu <- calculate_DTU(pseudo_sleuth, mini_anno, "Col", "Vir")

  expect_type(dtu, "list")
  expect_equal(length(dtu), 3)
  expect_identical(names(dtu), c("Comparison", "Genes", "Transcripts"))

  expect_true(is.vector(dtu$Comparison))
  expect_type(dtu$Comparison, "character")
  expect_equal(length(dtu$Comparison), 3)
  expect_identical(names(dtu$Comparison), c("variable_name", "reference", "compared"))

  expect_true(is.data.frame(dtu$Genes))
  expect_equal(dim(dtu$Genes)[2], 4)
  expect_identical(names(dtu$Genes), c("considered", "parent_id", "dtu", "p_value"))
  expect_true(is.logical(dtu$Genes$considered))
  expect_true(is.factor(dtu$Genes$parent_id))
  expect_true(is.logical(dtu$Genes$dtu))
  expect_true(is.numeric(dtu$Genes$p_value))

  expect_true(is.data.frame(dtu$Transcripts))
  expect_equal(dim(dtu$Transcripts)[2], 11)
  expect_identical(names(dtu$Transcripts), c("considered", "target_id", "parent_id",
                                             "ref_proportion", "comp_proportion", "ref_sum", "comp_sum",
                                             "ref_mean", "ref_variance", "comp_mean", "comp_variance" ))
  expect_true(is.logical(dtu$Transcripts$considered))
  expect_true(is.factor(dtu$Transcripts$target_id))
  expect_true(is.factor(dtu$Transcripts$parent_id))
  expect_true(is.numeric(dtu$Transcripts$ref_proportion))
  expect_true(is.numeric(dtu$Transcripts$comp_proportion))
  expect_true(is.numeric(dtu$Transcripts$ref_sum))
  expect_true(is.numeric(dtu$Transcripts$comp_sum))
  expect_true(is.numeric(dtu$Transcripts$ref_mean))
  expect_true(is.numeric(dtu$Transcripts$ref_variance))
  expect_true(is.numeric(dtu$Transcripts$comp_mean))
  expect_true(is.numeric(dtu$Transcripts$comp_variance))
})

#==============================================================================
test_that("The data munging is correct", {
  sl <- pseudo_sleuth
  ids <- mini_anno
  counts_col <- "est_counts"
  varname <- "condition"
  # Check that the intermediate values are correct.

  targets_by_parent <- parent_to_targets(ids)
  # all the parents are here:
  expect_identical(ordered(levels(as.factor(ids$parent_id))), ordered(names(targets_by_parent)))
  # all the targets are under their parent:
  for (target in ids$target_id) {
    expect_true(any(targets_by_parent[[ids[ids$target_id == target, "parent_id"]]] == target))
  }

  tx_filter <- mark_sibling_targets3(ids, targets_by_parent)
  # all the parents are here:
  expect_identical(ordered(levels(as.factor(ids$parent_id))), ordered(levels(as.factor(tx_filter$parent_id))))
  # all the targets are here:
  expect_identical(ordered(levels(as.factor(ids$target_id))), ordered(levels(as.factor(tx_filter$target_id))))
  # count parent duplicates and make sure flag was assigned correctly
  parent_counts <- table(ids$parent_id)
  for (parent in parent_counts) {
    expect_true(all((parent_counts[[parent]] > 1) == (tx_filter$has_siblings[tx_filter$parent_id == parent])))  # lengths may differ
  }

#   samples_by_condition <- group_samples(sl$sample_to_covariates)
#   # all the variables are here:
#   expect_identical(names(sl$sample_to_covariates), names(samples_by_condition))
#   # all the variable values are here:
#   for (cond in sl$sample_to_covariates){
#     expect_true(any( ordered(names(samples_by_condition[[cond]])) == ordered(levels(as.factor(sl$sample_to_covariates[[cond]]))) ))
#     # the correct number of samples have been assigned:
#     cond_counts <- table(sl$sample_to_covariates[[cond]])
#     for (subcond in samples_by_condition[[cond]]) {
#       expect_equal(cond_counts[sub_cont], length(samples_by_condition[[cond]][[sub_cond]]))
#       # the samples assigned are the correct ones:
#       for (sample in samples_by_condition[[cond]][[sub_cond]]) {
#         expect_equal(sl$sample_to_covariates[sample, cond], sub_cond)
#       }
#     }
#   }
#   samples_by_condition <- samples_by_condition[[varname]]  # focus on one covariate for the rest of the test.
#
#   count_data <- lapply(samples_by_condition, function(condition) make_filtered_bootstraps(sl, condition, tx_filter, counts_col))
#   for (condition in count_data) {
#     # all the targets are here:
#     expect_identical(ordered(rownames(condition)), ordered(tx_filter$target_id[tx_filter$has_siblings]))
#     # the correct number of counts have been gathered:
#     tot_boots <- 0
#     # the correct counts have been assigned:
#     for (sample in samples_by_condition[[condition]]) {
#       for (boot in sl$kal[[sample]]) {
#         for (target in rownames(condition)) {
#           expect_true(boot[boot$target_id == target, counts_col] %in% condition[target, ])
#         }
#         tot_boots <- tot_boots + 1
#       }
#     }
#     expect_equal(tot_boots, dim(condition)[2])
#   }
})

# #==============================================================================
# check_equal <- function(result1, result2) {
#
#   expect_equal(sort(result1$Col$mean), sort(result2$Col$mean), tolerance = 1e-6)
#   expect_equal(sort(result1$Col$variance), sort(result2$Col$variance), tolerance = 1e-6)
#   expect_equal(sort(result1$Col$proportion), sort(result2$Col$proportion), tolerance = 1e-6)
#   expect_equal(sort(result1$Vir$mean), sort(result2$Vir$mean), tolerance = 1e-6)
#   expect_equal(sort(result1$Vir$variance), sort(result2$Vir$variance), tolerance = 1e-6)
#   expect_equal(sort(result1$Vir$proportion), sort(result2$Vir$proportion), tolerance = 1e-6)
# }
#
# test_that("Mixed order bootstraps give same results as unmixed", {
#
#   # make a pseudo sleuth object with mixed up bootstrap rows
#   mixed_pseudo_sleuth <- pseudo_sleuth
#   mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]][c(1,3,2,5,4),]
#   mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]][c(5,4,3,2,1),]
#   mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]] <- mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]][c(3,5,1,4,2),]
#
#   mixed_stats <- calculate_DTU(mixed_pseudo_sleuth, mini_anno)
#   unmixed_stats <- calculate_DTU(pseudo_sleuth, mini_anno)
#
#   check_equal(mixed_stats, unmixed_stats)
# })
#
# #==============================================================================
# test_that("Bootstraps with all 0 / NA entries are discarded", {
#
#   # add a new transcript which will have zero entries
#   TRANSCRIPT <- "AT1G01020.3"
#   GENE <- "AT1G01020"
#   z_mini_anno <- rbind(mini_anno, c(TRANSCRIPT, GENE,""))
#
#   # make a new sleuth with this transcript added, and give it 0 entries
#   zero_entry_sleuth <- pseudo_sleuth
#   zero_entry_sleuth$kal <- lapply(zero_entry_sleuth$kal, function(kal)
#     lapply(kal, function(bs) lapply(bs, rbind, list(TRANSCRIPT, 0, 0))))
#
#   # make a new sleuth with one of the zero entries set to NA
#   na_entry_sleuth <- zero_entry_sleuth
#   na_entry_sleuth$kal[[2]]$bootstrap[[1]][6,2] <- NA
#
#   # calculate DTU for new sleuths and old one
#   z_result <- calculate_DTU(zero_entry_sleuth, z_mini_anno)
#   n_result <- calculate_DTU(na_entry_sleuth, z_mini_anno)
#   result <- calculate_DTU(pseudo_sleuth, mini_anno)
#
#   # check results are all equal
#   check_equal(z_result, result)
#   check_equal(n_result, result)
# })
