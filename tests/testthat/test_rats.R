library(rats)

context("DTU input controls.")

#==============================================================================
test_that("The input checks work", {
  # Sensible values.
  slo <- pseudo_sleuth
  annot <- mini_anno
  name_A <- "Col"
  name_B <- "Vir"
  varname="condition"
  p_thresh = 0.05
  count_thresh = 5
  testmode = "both"
  correction = "BH"
  verbose = FALSE
  boots = "none"
  bootnum = 10L
  threads = parallel::detectCores()
  COUNTS_COL="est_counts" 
  TARGET_COL="target_id" 
  PARENT_COL="parent_id"
  BS_TARGET_COL="target_id"
  
  wrong_name <- "JUNK_TEST_NAME"

  # No false alarms.
  expect_silent(calculate_DTU(slo, annot, name_A, name_B))
  expect_silent(calculate_DTU(slo, annot, name_A, name_B, varname = varname, p_thresh = p_thresh, count_thresh = count_thresh,
                              testmode = testmode, correction = correction, verbose = verbose, boots = boots,
                              bootnum = bootnum, threads = threads, COUNTS_COL = COUNTS_COL, TARGET_COL = TARGET_COL, 
                              PARENT_COL = PARENT_COL, BS_TARGET_COL = BS_TARGET_COL))

  # Annottaion is not a dataframe.
  expect_error(calculate_DTU(slo, c("not", "a", "dataframe"), name_A, name_B), "annot is not a data.frame.")
  # Annotation field names.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, TARGET_COL=wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed=TRUE)
  expect_error(calculate_DTU(slo, annot, name_A, name_B, PARENT_COL=wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed=TRUE)
  
  # Bootstrap field names.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, BS_TARGET_COL=wrong_name),
               "target IDs field-name does not exist in the bootstraps", fixed=TRUE)
  expect_error(calculate_DTU(slo, annot, name_A, name_B, COUNTS_COL=wrong_name),
               "counts field-name does not exist in the bootstraps", fixed=TRUE)

  # Correction method.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, correction=wrong_name),
               "Invalid p-value correction method name", fixed=TRUE)

  # Covariate name.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, varname=wrong_name),
               "covariate name does not exist", fixed=TRUE)

  # Condition names.
  expect_error(calculate_DTU(slo, annot, wrong_name, name_B),
               "conditions do not exist", fixed=TRUE)
  expect_error(calculate_DTU(slo, annot, name_A, wrong_name),
               "conditions do not exist", fixed=TRUE)
  
  # Verbose is bool.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, verbose="yes"),
               "verbose must be a logical", fixed=TRUE)
  
  # Probability threshold.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, p_thresh = 666),
               "Invalid p-value threshold", fixed=TRUE)
  expect_error(calculate_DTU(slo, annot, name_A, name_B, p_thresh = -0.05),
               "Invalid p-value threshold", fixed=TRUE)
  
  # Read counts threshold.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, count_thresh = -5),
               "Invalid read-count threshold", fixed=TRUE)
  
  # Tests.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, testmode="GCSE"),
               "Unrecognized value for testmode", fixed=TRUE)
  expect_silent(calculate_DTU(slo, annot, name_A, name_B, testmode="g-test"))
  expect_silent(calculate_DTU(slo, annot, name_A, name_B, testmode="prop-test"))
  
  expect_error(calculate_DTU(slo, annot, name_A, name_B, boots="GCSE"),
               "Unrecognized value for boots", fixed=TRUE)
  expect_silent(calculate_DTU(slo, annot, name_A, name_B, boots="g-test"))
  expect_silent(calculate_DTU(slo, annot, name_A, name_B, boots="prop-test"))
  
  # Number of bootstraps.
  expect_error(calculate_DTU(slo, annot, name_A, name_B, bootnum = -5),
               "Invalid number of bootstraps", fixed=TRUE)
})



context("DTU results")

#==============================================================================
test_that("The output structure is correct", {
  dtu <- calculate_DTU(pseudo_sleuth, mini_anno, "Col", "Vir")

  expect_type(dtu, "list")
  expect_equal(length(dtu), 3)
  expect_named(dtu, c("Parameters", "Genes", "Transcripts"))

  expect_type(dtu$Parameters, "list")
  expect_length(dtu$Parameters, 11)
  expect_named(dtu$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                 "count_thresh", "tests", "bootstrap", "bootnum", "threads"))

  expect_true(is.data.frame(dtu$Genes))
  expect_equal(dim(dtu$Genes)[2], 23)
  expect_named(dtu$Genes, c("parent_id", "known_transc", "detect_transc", "eligible", "Pt_DTU", "Gt_DTU", 
                            "Gt_dtuAB", "Gt_dtuBA", "Gt_pvalAB", "Gt_pvalBA", "Gt_pvalAB_corr", "Gt_pvalBA_corr", 
                            "Gt_boot_dtuAB", "Gt_boot_dtuBA", "Gt_boot_meanAB", "Gt_boot_meanBA", "Gt_boot_stdevAB", 
                            "Gt_boot_stdevBA", "Gt_boot_minAB", "Gt_boot_minBA", "Gt_boot_maxAB", "Gt_boot_maxBA", 
                            "Gt_boot_na"))
  expect_true(is.numeric(dtu$Genes$known_transc))
  expect_true(is.numeric(dtu$Genes$detect_transc))
  expect_true(is.logical(dtu$Genes$eligible))
  expect_true(is.logical(dtu$Genes$Pt_DTU))
  expect_true(is.logical(dtu$Genes$Gt_DTU))
  expect_true(is.logical(dtu$Genes$Gt_dtuAB))
  expect_true(is.logical(dtu$Genes$Gt_dtuBA))
  expect_true(is.numeric(dtu$Genes$Gt_pvalAB))
  expect_true(is.numeric(dtu$Genes$Gt_pvalBA))
  expect_true(is.numeric(dtu$Genes$Gt_pvalAB_corr))
  expect_true(is.numeric(dtu$Genes$Gt_pvalBA_corr))
  expect_true(is.numeric(dtu$Genes$Gt_boot_dtuAB))
  expect_true(is.numeric(dtu$Genes$Gt_boot_dtuBA))
  expect_true(is.numeric(dtu$Genes$Gt_boot_meanAB))
  expect_true(is.numeric(dtu$Genes$Gt_boot_meanBA))
  expect_true(is.numeric(dtu$Genes$Gt_boot_stdevAB))
  expect_true(is.numeric(dtu$Genes$Gt_boot_stdevBA))
  expect_true(is.numeric(dtu$Genes$Gt_boot_minAB))
  expect_true(is.numeric(dtu$Genes$Gt_boot_minBA))
  expect_true(is.numeric(dtu$Genes$Gt_boot_maxAB))
  expect_true(is.numeric(dtu$Genes$Gt_boot_maxBA))
  expect_true(is.numeric(dtu$Genes$Gt_boot_na))
  
  expect_true(is.data.frame(dtu$Transcripts))
  expect_equal(dim(dtu$Transcripts)[2], 24)
  expect_named(dtu$Transcripts, c("target_id", "parent_id", "propA", "propB", "Dprop", "eligible", "Gt_DTU", "Pt_DTU", 
                                  "Pt_pval", "Pt_pval_corr", "Pt_boot_dtu", "Pt_boot_mean", "Pt_boot_stdev", 
                                  "Pt_boot_min", "Pt_boot_max", "Pt_boot_na",
                                  "sumA", "sumB", "meanA", "meanB", "stdevA", "stdevB", "totalA", "totalB"))
  expect_true(is.numeric(dtu$Transcripts$prop_A))
  expect_true(is.numeric(dtu$Transcripts$prop_B))
  expect_true(is.numeric(dtu$Transcripts$Dprop))
  expect_true(is.logical(dtu$Transcripts$eligible))
  expect_true(is.logical(dtu$Transcripts$Gt_DTU))
  expect_true(is.logical(dtu$Transcripts$Pt_DTU))
  expect_true(is.numeric(dtu$Transcripts$Pt_pval))
  expect_true(is.numeric(dtu$Transcripts$Pt_pva_corrl))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_dtu))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_mean))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_stdev))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_min))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_max))
  expect_true(is.numeric(dtu$Transcripts$Pt_boot_na))
  expect_true(is.numeric(dtu$Transcripts$sum_A))
  expect_true(is.numeric(dtu$Transcripts$sum_B))
  expect_true(is.numeric(dtu$Transcripts$mean_A))
  expect_true(is.numeric(dtu$Transcripts$mean_B))
  expect_true(is.numeric(dtu$Transcripts$stdev_A))
  expect_true(is.numeric(dtu$Transcripts$stdev_B))
  expect_true(is.numeric(dtu$Transcripts$total_A))
  expect_true(is.numeric(dtu$Transcripts$total_B))
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
