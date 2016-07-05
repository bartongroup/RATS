sim_sleuth_data <- function(varnames=c("condition","foobar"), COUNTS_COL="est_counts", TARGET_COL="target_id" , PARENT_COL="parent_id", BS_TARGET_COL="target_id", cnames=c("one","two")) {
  tx <- data.frame("target_id" <- c("a","b","c","d","e"), "parent_id"=c("AC","B","AC","DE","DE"))
  names(tx) <- c(TARGET_COL, PARENT_COL)
  
  sl <- list()
  sl[["sample_to_covariates"]] <- data.frame("foo"=c(cnames[1],cnames[2],cnames[1],cnames[2]), "bar"=c("ba","ba", "bb","bb"))
  names(sl[["sample_to_covariates"]]) <- varnames
  
  sl[["kal"]] <- list()
  sl$kal[[1]] <- list()
  sl$kal[[1]]["bootstrap"] <- list()
  sl$kal[[1]]$bootstrap[[1]] <- data.frame("target"=c("e","c","a","b","d"), "my_counts"=c(15,13,11,12,14))
  sl$kal[[1]]$bootstrap[[2]] <- data.frame("target"=c("a","c","e","d","b"), "my_counts"=c(21,23,25,24,22))
  sl$kal[[1]]$bootstrap[[3]] <- data.frame("target"=c("e","d","c","b","a"), "my_counts"=c(35,34,33,32,31))
  
  sl$kal[[2]] <- list()
  sl$kal[[2]]["bootstrap"] <- list()
  sl$kal[[2]]$bootstrap[[1]] <- data.frame("target"=c("e","c","a","b","d"), "my_counts"=c(500,300,100,200,400))
  sl$kal[[2]]$bootstrap[[2]] <- data.frame("target"=c("a","c","e","d","b"), "my_counts"=c(1100,1300,1500,1400,1200))
  sl$kal[[2]]$bootstrap[[3]] <- data.frame("target"=c("e","d","c","b","a"), "my_counts"=c(2500,2400,2300,2200,2100))
  
  sl$kal[[3]] <- list()
  sl$kal[[3]]["bootstrap"] <- list()
  sl$kal[[3]]$bootstrap[[1]] <- data.frame("target"=c("c","a","d","b","e"), "my_counts"=c(103,101,104,102,105))
  sl$kal[[3]]$bootstrap[[2]] <- data.frame("target"=c("e","b","a","d","c"), "my_counts"=c(115,112,111,114,113))
  sl$kal[[3]]$bootstrap[[3]] <- data.frame("target"=c("b","a","c","e","d"), "my_counts"=c(122,121,123,125,124))
  
  sl$kal[[4]] <- list()
  sl$kal[[4]]["bootstrap"] <- list()
  sl$kal[[4]]$bootstrap[[1]] <- data.frame("target"=c("e","c","a","b","d"), "my_counts"=c(50000,30000,10000,20000,4000))
  sl$kal[[4]]$bootstrap[[2]] <- data.frame("target"=c("a","c","e","d","b"), "my_counts"=c(11000,13000,15000,14000,12000))
  sl$kal[[4]]$bootstrap[[3]] <- data.frame("target"=c("e","d","c","b","a"), "my_counts"=c(25000,24000,23000,22000,21000))
  
  names(sl$kal[[1]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[1]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[2]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[3]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[1]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[2]]) <- c(BS_TARGET_COL, COUNTS_COL)
  names(sl$kal[[4]]$bootstrap[[3]]) <- c(BS_TARGET_COL, COUNTS_COL)
  
  return(list("annot" = tx, "slo" = sl))
}


context("DTU Input checks.")

#==============================================================================
test_that("The input checks don't work", {
  name_A <- "one"
  name_B <- "two"
  wrong_name <- "RUBBISH_COLUMN_NAME"
  
  # No false alarms.
  sim <- sim_sleuth_data(varnames=c("foo","bar"), COUNTS_COL="counts", TARGET_COL="target" , PARENT_COL="parent", BS_TARGET_COL="id", cnames=c("AAAA","BBBB"))
  expect_silent(calculate_DTU(sim$slo, sim$annot, "AAAA", "BBBB", varname = "foo", p_thresh = 0.01, count_thresh = 10,
                              testmode = "prop-test", correction = "bonferroni", verbose = FALSE, boots = "g-test",
                              bootnum = 2, threads = 1, COUNTS_COL = "counts", TARGET_COL = "target", 
                              PARENT_COL = "parent", BS_TARGET_COL = "id"))
  sim <- sim_sleuth_data(cnames=c(name_A, name_B))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B))

  # Annottaion is not a dataframe.
  expect_error(calculate_DTU(sim$slo, c("not", "a", "dataframe"), name_A, name_B), "annot is not a data.frame.")
  # Annotation field names.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, TARGET_COL=wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed=TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, PARENT_COL=wrong_name),
               "target and/or parent IDs field-names do not exist in annot", fixed=TRUE)
  
  # Bootstrap field names.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, BS_TARGET_COL=wrong_name),
               "target IDs field-name does not exist in the bootstraps", fixed=TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, COUNTS_COL=wrong_name),
               "counts field-name does not exist", fixed=TRUE)

  # Correction method.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, correction=wrong_name),
               "Invalid p-value correction method name", fixed=TRUE)

  # Covariate name.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, varname=wrong_name),
               "covariate name does not exist", fixed=TRUE)

  # Condition names.
  expect_error(calculate_DTU(sim$slo, sim$annot, wrong_name, name_B),
               "conditions do not exist", fixed=TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, wrong_name),
               "conditions do not exist", fixed=TRUE)
  
  # Verbose is bool.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, verbose="yes"),
               "verbose must be a logical", fixed=TRUE)
  
  # Probability threshold.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = 666),
               "Invalid p-value threshold", fixed=TRUE)
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, p_thresh = -0.05),
               "Invalid p-value threshold", fixed=TRUE)
  
  # Read counts threshold.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, count_thresh = -5),
               "Invalid read-count threshold", fixed=TRUE)
  
  # Tests.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="GCSE"),
               "Unrecognized value for testmode", fixed=TRUE)
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="g-test"))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, testmode="prop-test"))
  
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="GCSE"),
               "Unrecognized value for boots", fixed=TRUE)
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="g-test", bootnum = 2))
  expect_silent(calculate_DTU(sim$slo, sim$annot, name_A, name_B, boots="prop-test", bootnum = 2))
  
  # Number of bootstraps.
  expect_error(calculate_DTU(sim$slo, sim$annot, name_A, name_B, bootnum = -5),
               "Invalid number of bootstraps", fixed=TRUE)
})


context("DTU Output")

#==============================================================================
test_that("The output structure is not correct", {
  sim <- sim_sleuth_data(cnames=c("ONE","TWO"))
  full <- calculate_DTU(sim$slo, sim$annot, "ONE", "TWO", boots="both", bootnum=2)

  expect_type(full, "list")
  expect_equal(length(full), 3)
  expect_named(full, c("Parameters", "Genes", "Transcripts"))

  expect_type(full$Parameters, "list")
  expect_length(full$Parameters, 11)
  expect_named(full$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "tests", "bootstrap", "bootnum", "threads"))
  
  expect_true(is.data.frame(full$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_named(full$Genes, c("parent_id", "known_transc", "detect_transc", "eligible", "Pt_DTU", "Gt_DTU", 
                             "Gt_dtuAB", "Gt_dtuBA", "Gt_pvalAB", "Gt_pvalBA", "Gt_pvalAB_corr", "Gt_pvalBA_corr", 
                             "Gt_boot_dtuAB", "Gt_boot_dtuBA", "Gt_boot_meanAB", "Gt_boot_meanBA", "Gt_boot_stdevAB", 
                             "Gt_boot_stdevBA", "Gt_boot_minAB", "Gt_boot_minBA", "Gt_boot_maxAB", "Gt_boot_maxBA", 
                             "Gt_boot_na"))
  expect_true(is.numeric(full$Genes[["known_transc"]]))
  expect_true(is.numeric(full$Genes[["detect_transc"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalAB_corr"]]))
  expect_true(is.numeric(full$Genes[["Gt_pvalBA_corr"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_dtuAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_dtuBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_meanAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_meanBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_stdevAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_stdevBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_minAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_minBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_maxAB"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_maxBA"]]))
  expect_true(is.numeric(full$Genes[["Gt_boot_na"]]))
  expect_true(is.logical(full$Genes[["eligible"]]))
  expect_true(is.logical(full$Genes[["Pt_DTU"]]))
  expect_true(is.logical(full$Genes[["Gt_DTU"]]))
  expect_true(is.logical(full$Genes[["Gt_dtuAB"]]))
  expect_true(is.logical(full$Genes[["Gt_dtuBA"]]))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 24)
  expect_named(full$Transcripts, c("target_id", "parent_id", "propA", "propB", "Dprop", "eligible", "Gt_DTU", "Pt_DTU", 
                                   "Pt_pval", "Pt_pval_corr", "Pt_boot_dtu", "Pt_boot_mean", "Pt_boot_stdev", 
                                   "Pt_boot_min", "Pt_boot_max", "Pt_boot_na",
                                   "sumA", "sumB", "meanA", "meanB", "stdevA", "stdevB", "totalA", "totalB"))
  expect_true(is.logical(full$Transcripts[["eligible"]]))
  expect_true(is.logical(full$Transcripts[["Gt_DTU"]]))
  expect_true(is.logical(full$Transcripts[["Pt_DTU"]]))
  expect_true(is.numeric(full$Transcripts[["propA"]]))
  expect_true(is.numeric(full$Transcripts[["propB"]]))
  expect_true(is.numeric(full$Transcripts[["Dprop"]]))
  expect_true(is.numeric(full$Transcripts[["sumA"]]))
  expect_true(is.numeric(full$Transcripts[["sumB"]]))
  expect_true(is.numeric(full$Transcripts[["meanA"]]))
  expect_true(is.numeric(full$Transcripts[["meanB"]]))
  expect_true(is.numeric(full$Transcripts[["stdevA"]]))
  expect_true(is.numeric(full$Transcripts[["stdevB"]]))
  expect_true(is.numeric(full$Transcripts[["totalA"]]))
  expect_true(is.numeric(full$Transcripts[["totalB"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_pval"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_pval_corr"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_dtu"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_mean"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_stdev"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_min"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_max"]]))
  expect_true(is.numeric(full$Transcripts[["Pt_boot_na"]]))
})


context("DTU internal")

#==============================================================================
test_that("The reporting structures are not created correctly", {
  full <- alloc_out(mini_anno, "full")
  short <- alloc_out(mini_anno, "short")
  
  expect_type(full, "list")
  expect_type(short, "list")
  expect_equal(length(full), 3)
  expect_equal(length(full), length(short))
  expect_named(full, c("Parameters", "Genes", "Transcripts"))
  expect_true(all(names(full) == names(short)))
  
  expect_type(full$Parameters, "list")
  expect_true(typeof(full$Parameters) == typeof(short$Parameters))
  expect_length(full$Parameters, 11)
  expect_length(short$Parameters, 2)
  expect_named(full$Parameters, c("var_name", "cond_A", "cond_B", "num_replic_A", "num_replic_B", "p_thresh", 
                                  "count_thresh", "tests", "bootstrap", "bootnum", "threads"))
  expect_true(all(names(short$Parameters) %in% names(full$Parameters)))
  
  expect_true(is.data.frame(full$Genes))
  expect_true(is.data.frame(short$Genes))
  expect_equal(dim(full$Genes)[2], 23)
  expect_equal(dim(short$Genes)[2], 7)
  expect_named(full$Genes, c("parent_id", "known_transc", "detect_transc", "eligible", "Pt_DTU", "Gt_DTU", 
                             "Gt_dtuAB", "Gt_dtuBA", "Gt_pvalAB", "Gt_pvalBA", "Gt_pvalAB_corr", "Gt_pvalBA_corr", 
                             "Gt_boot_dtuAB", "Gt_boot_dtuBA", "Gt_boot_meanAB", "Gt_boot_meanBA", "Gt_boot_stdevAB", 
                             "Gt_boot_stdevBA", "Gt_boot_minAB", "Gt_boot_minBA", "Gt_boot_maxAB", "Gt_boot_maxBA", 
                             "Gt_boot_na"))
  expect_true(all(names(short$Genes) %in% names(full$Genes)))
  
  expect_true(is.data.frame(full$Transcripts))
  expect_true(is.data.frame(short$Transcripts))
  expect_equal(dim(full$Transcripts)[2], 24)
  expect_equal(dim(short$Transcripts)[2], 12)
  expect_named(full$Transcripts, c("target_id", "parent_id", "propA", "propB", "Dprop", "eligible", "Gt_DTU", "Pt_DTU", 
                                   "Pt_pval", "Pt_pval_corr", "Pt_boot_dtu", "Pt_boot_mean", "Pt_boot_stdev", 
                                   "Pt_boot_min", "Pt_boot_max", "Pt_boot_na",
                                   "sumA", "sumB", "meanA", "meanB", "stdevA", "stdevB", "totalA", "totalB"))
  expect_true(all(names(short$Transcripts) %in% names(full$Transcripts)))
})


#==============================================================================
test_that("Samples are not grouped correctly", {
  sim <- sim_sleuth_data()
  r <- group_samples(sim$slo$sample_to_covariates)
  
  expect_equal(length(r), length(sim$slo$sample_to_covariates))  # number of covariates
  expect_named(r, names(sim$slo$sample_to_covariates))  # names of covariates
  expect_equal( length(r[[1]]), length(levels(as.factor(sim$slo$sample_to_covariates[[1]]))) )  # number of values of each covariate
  expect_equal( length(r[[2]]), length(levels(as.factor(sim$slo$sample_to_covariates[[2]]))) )
  expect_equal( sum(sapply(r[[1]], length)), length(sim$slo$sample_to_covariates[[1]]) )  # total number of samples
  expect_equal( sum(sapply(r[[2]], length)), length(sim$slo$sample_to_covariates[[2]]) )
})


#==============================================================================
test_that("Bootstrapped counts are not extracted correctly", {
  samples <- c(1,3)
  sim <- sim_sleuth_data(COUNTS_COL = "counts", BS_TARGET_COL = "id")
  lr <- denest_boots(sim$slo, sim$annot[[1]], samples, "counts", "id")
  
  expect_equal(length(lr), length(samples))
  # browser()
})


#==============================================================================
test_that("Mixed order bootstraps don't give same results as unmixed", {
  sim <- sim_sleuth_data(cnames=c("one","two"))
  mixed_pseudo_sleuth <- sim$slo
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[3]][c(1,3,2,5,4),]
  mixed_pseudo_sleuth$kal[[1]]$bootstrap[[2]] <- mixed_pseudo_sleuth$kal[[1]]$bootstrap[[2]][c(5,4,3,2,1),]
  mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]] <- mixed_pseudo_sleuth$kal[[3]]$bootstrap[[1]][c(3,5,1,4,2),]

  expect_equal(calculate_DTU(mixed_pseudo_sleuth, sim$annot, "one", "two"),
               calculate_DTU(sim$slo, sim$annot, "one", "two"))
})


#==============================================================================
test_that("Bootstraps with all 0 / NA entries are not filtered out", {
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


