library(rats)

context("DTU input controls.")

#==============================================================================
test_that("The input checks work", {
  sl <- pseudo_sleuth
  ids <- mini_anno
  counts_col <- "est_counts"
  varname <- "condition"
  TARGET_ID <- "target_id"
  BS_TARGET_ID <- "target_id"
  PARENT_ID <- "parent_id"
  ref <- "Col"
  comp <- "Vir"
  wrong_name <- c("JUNK_A", "JUNK_B", "JUNK_C", "JUNK_D", "JUNK_E", "JUNK_F")

  # No false alarms.
  expect_silent(calculate_DTU(sl, ids, ref, comp))
  expect_silent(calculate_DTU(sl, ids, ref, comp, counts_col=counts_col, varname=varname,
                              TARGET_ID=TARGET_ID, PARENT_ID=PARENT_ID, BS_TARGET_ID=BS_TARGET_ID))

  # Transcript ids are not a dataframe.
  expect_error(calculate_DTU(sl, c("not", "a", "dataframe"), ref, comp), "transcripts is not a data.frame.")

  # Transcript field names.
  expect_error(calculate_DTU(sl, ids, ref, comp, TARGET_ID=wrong_name[1]),
               "The specified target and parent IDs field-names do not exist in transcripts.", fixed=TRUE)
  expect_error(calculate_DTU(sl, ids, ref, comp, PARENT_ID=wrong_name[2]),
               "The specified target and parent IDs field-names do not exist in transcripts.", fixed=TRUE)
  expect_error(calculate_DTU(sl, ids, ref, comp, TARGET_ID=wrong_name[1], PARENT_ID=wrong_name[2]),
               "The specified target and parent IDs field-names do not exist in transcripts.", fixed=TRUE)

  # Bootstrap field names.
  expect_error(calculate_DTU(sl, ids, ref, comp, BS_TARGET_ID=wrong_name[1]),
               "The specified target IDs field-name does not exist in the bootstraps.", fixed=TRUE)
  expect_error(calculate_DTU(sl, ids, ref, comp, counts_col=wrong_name[3]),
               "The specified counts field-name does not exist.", fixed=TRUE)

  # Correction method.
  expect_error(calculate_DTU(sl, ids, ref, comp, correction=wrong_name[4]),
               "Invalid p-value correction method name. Refer to stats::p.adjust.methods.", fixed=TRUE)

  # Covariate name.
  expect_error(calculate_DTU(sl, ids, ref, comp, varname=wrong_name[5]),
               "The specified covariate name does not exist.", fixed=TRUE)

  # Condition names.
  expect_error(calculate_DTU(sl, ids, wrong_name[6], comp),
               "One or both of the specified conditions do not exist.", fixed=TRUE)
  expect_error(calculate_DTU(sl, ids, ref, wrong_name[7]),
               "One or both of the specified conditions do not exist.", fixed=TRUE)
  expect_error(calculate_DTU(sl, ids, wrong_name[6], wrong_name[7]),
               "One or both of the specified conditions do not exist.", fixed=TRUE)

  # Verbose is bool.
  expect_error(calculate_DTU(sl, ids, ref, comp, verbose="yes"),
               "verbose must be a logical value.", fixed=TRUE)
})

context("DTU results")

#==============================================================================
test_that("The output structure is correct", {
  dtu <- calculate_DTU(pseudo_sleuth, mini_anno, "Col", "Vir")

  expect_type(dtu, "list")
  expect_equal(length(dtu), 3)
  expect_named(dtu, c("Parameters", "Genes", "Transcripts"))

  expect_type(dtu$Parameters, "list")
  expect_length(dtu$Parameters, 6)
  expect_named(dtu$Parameters, c("var_name", "cond_A", "cond_B", "replicates_A", "replicates_B", "p_thresh"))

  expect_true(is.data.frame(dtu$Genes))
  expect_equal(dim(dtu$Genes)[2], 10)
  expect_named(dtu$Genes, c("parent_id", "known_transc", "usable_transc", "pval_AB", "pval_BA",
                            "pval_AB_corr", "pval_BA_corr", "dtu_AB", "dtu_BA", "dtu"))
  expect_true(is.numeric(dtu$Genes$known_transc))
  expect_true(is.numeric(dtu$Genes$usable_transc))
  expect_true(is.numeric(dtu$Genes$pval_AB))
  expect_true(is.numeric(dtu$Genes$pval_BA))
  expect_true(is.numeric(dtu$Genes$pval_AB_corr))
  expect_true(is.numeric(dtu$Genes$pval_BA_corr))
  expect_true(is.logical(dtu$Genes$dtu_AB))
  expect_true(is.logical(dtu$Genes$dtu_BA))
  expect_true(is.logical(dtu$Genes$dtu))

  expect_true(is.data.frame(dtu$Transcripts))
  expect_equal(dim(dtu$Transcripts)[2], 10)
  expect_named(dtu$Transcripts, c("target_id", "parent_id",
                                  "prop_A", "prop_B", "sum_A", "sum_B",
                                  "mean_A", "mean_B", "stdev_A", "stdev_B"))
  expect_true(is.numeric(dtu$Transcripts$prop_A))
  expect_true(is.numeric(dtu$Transcripts$prop_B))
  expect_true(is.numeric(dtu$Transcripts$sum_A))
  expect_true(is.numeric(dtu$Transcripts$sum_B))
  expect_true(is.numeric(dtu$Transcripts$mean_A))
  expect_true(is.numeric(dtu$Transcripts$mean_B))
  expect_true(is.numeric(dtu$Transcripts$stdev_A))
  expect_true(is.numeric(dtu$Transcripts$stdev_B))
})

#==============================================================================
test_that("The data munging is correct", {
  sl <- pseudo_sleuth
  transcripts <- mini_anno
  counts_col <- "est_counts"
  varname <- "condition"
  TARGET_ID <- "target_id"
  BS_TARGET_ID <- "target_id"
  PARENT_ID <- "parent_id"
  # Check that the intermediate values are correct.

  targets <- data.table("target_id"=transcripts[[TARGET_ID]], "parent_id"=transcripts[[PARENT_ID]])
  setkey(targets, target_id)
  parents <- data.table(targets)
  setkey(parents, parent_id)
  
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

  count_data <- lapply(samples_by_condition, function(condition) make_filtered_bootstraps(sl, condition, targets, counts_col, BS_TARGET_ID))
  for (cond_name in names(count_data)) {
    # all the targets are here:
    expect_identical(ordered(rownames(count_data[[cond_name]])), ordered(targets$target_id))
    # the number of replicates is correct
    expect_equal(dim(count_data[[cond_name]])[2], length(sl$kal[samples_by_condition[[cond_name]]]))
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
  result <- calculate_DTU(pseudo_sleuth, z_mini_anno, "Col", "Vir")

  # check results are all equal
  expect_equal(z_result, result)
  expect_equal(n_result, result)
})
