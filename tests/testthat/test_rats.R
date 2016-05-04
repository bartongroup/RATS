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

  targets_by_parent <- split(as.matrix(ids[TARGET_ID]), ids[[PARENT_ID]])
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
  for (parent in names(parent_counts)) {
    expect_true(all((parent_counts[[parent]] > 1) == (tx_filter$has_siblings[tx_filter$parent_id == parent])))  # lengths may differ
  }

  samples_by_condition <- group_samples(sl$sample_to_covariates)
  # all the variables are here:
  expect_identical(names(sl$sample_to_covariates), names(samples_by_condition))
  # all the variable values are here:
  for (condname in names(sl$sample_to_covariates)){
    expect_true(any( ordered(names(samples_by_condition[[condname]])) == ordered(levels(as.factor(sl$sample_to_covariates[[condname]]))) ))
    # the correct number of samples have been assigned:
    cond_counts <- table(sl$sample_to_covariates[[condname]])
    for (sub_cond in names(samples_by_condition[[condname]])) {
      expect_equal(cond_counts[[sub_cond]], length(samples_by_condition[[condname]][[sub_cond]]))
      # the samples assigned are the correct ones:
      for (smpl in samples_by_condition[[condname]][[sub_cond]]) {
        expect_equal(sl$sample_to_covariates[smpl, condname], sub_cond)
      }
    }
  }
  samples_by_condition <- samples_by_condition[[varname]]  # focus on one covariate for the rest of the test.

  count_data <- lapply(samples_by_condition, function(condition) make_filtered_bootstraps(sl, condition, tx_filter, counts_col))
  for (cond_name in names(count_data)) {
    # all the targets are here:
    expect_identical(ordered(rownames(count_data[[cond_name]])), ordered(tx_filter$target_id[tx_filter$has_siblings]))
    # the correct number of counts have been gathered:
    tot_boots <- 0
    # the correct counts have been assigned:
    for (smpl in samples_by_condition[[cond_name]]) {
      for (boot in sl$kal[[smpl]]$bootstrap) {
        for (target in rownames(count_data[[cond_name]])) {
          expect_true(boot[boot$target_id == target, counts_col] %in% count_data[[cond_name]][target, ])
        }
        tot_boots <- tot_boots + 1
      }
    }
    expect_equal(tot_boots, dim(count_data[[cond_name]])[2])
  }
})

#==============================================================================
test_that("Mixed order bootstraps give same results as unmixed", {

  # make a pseudo sleuth object with mixed up bootstrap rows
  mixed_pseudo_sleuth <- pseudo_sleuth
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]][c(1,3,2,5,4),]
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[4]][c(5,4,3,2,1),]
  mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]] <- mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]][c(3,5,1,4,2),]

  expect_equal(calculate_DTU(mixed_pseudo_sleuth, mini_anno, "Col", "Vir"),
               calculate_DTU(pseudo_sleuth, mini_anno, "Col", "Vir"))
})

#==============================================================================
test_that("Bootstraps with all 0 / NA entries are discarded", {

  # add a new transcript which will have zero entries
  TRANSCRIPT <- "AT1G01020.3"
  GENE <- "AT1G01020"
  z_mini_anno <- rbind(mini_anno, c(TRANSCRIPT, GENE,""))

  # make a new sleuth with this transcript added, and give it 0 entries
  zero_entry_sleuth <- pseudo_sleuth
  zero_entry_sleuth$kal <- lapply(zero_entry_sleuth$kal, function(kal)
    lapply(kal, function(bs) lapply(bs, rbind, list(TRANSCRIPT, 0, 0))))

  # make a new sleuth with one of the zero entries set to NA
  na_entry_sleuth <- zero_entry_sleuth
  na_entry_sleuth$kal[[2]]$bootstrap[[1]][6,2] <- NA

  # calculate DTU for new sleuths and old one
  z_result <- calculate_DTU(zero_entry_sleuth, z_mini_anno, "Col", "Vir")
  n_result <- calculate_DTU(na_entry_sleuth, z_mini_anno, "Col", "Vir")
  result <- calculate_DTU(pseudo_sleuth, mini_anno, "Col", "Vir")

  # check results are all equal
  expect_equal(z_result, result)
  expect_equal(n_result, result)
})
