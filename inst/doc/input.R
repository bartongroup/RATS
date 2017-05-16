## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(rats)

## ---- echo=FALSE---------------------------------------------------------
# Show the first rows of the table corresponding to one sample, from simulated data.
head(sim_boot_data()[[2]][[1]])

## ---- echo=FALSE---------------------------------------------------------
# Show the first rows of the table corresponding to the annotation, from simulated data.
head(sim_count_data()[[1]])

## ----eval=FALSE----------------------------------------------------------
#  # Extract transcript ID to gene ID index from a GTF annotation.
#  myannot <- annot2ids("my_annotation_file.gtf")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_count_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]       # Simulated abundances for one condition.
mycond_B <- simdat[[3]]       # Simulated abundances for other condition.
myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.

## ------------------------------------------------------------------------
# Find DTU between the simulated datasets.
mydtu <- call_DTU(annot= myannot, count_data_A= mycond_A, count_data_B= mycond_B, 
                  verbose= FALSE,
                  name_A= "healthy", name_B= "patients", varname= "My phenotype",
                  description="Comparison of two simulated counts datasets for the tutorial. Simulated using built-in functionality of RATs.")

## ------------------------------------------------------------------------
# Simulate some data. (Notice it is a different function than before.)
simdat <- sim_boot_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]       # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]       # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.

## ------------------------------------------------------------------------
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot= myannot, boot_data_A= mycond_A, boot_data_B= mycond_B, 
                  verbose= FALSE, 
                  name_A= "wildtype", name_B= "some mutant", varname = "My phenotype",
                  description="Comparison of two simulated datasets of bootstrapped counts for the tutorial. Simulated using built-in functionality of RATs.")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# controls and patients are arbitrary names to use as conditions.

# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo       # Simulated minimal sleuth object.
myannot <- simdat$annot   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot= myannot, slo= myslo, name_A= "controls", name_B= "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object for the purposes of the tutorial. Simulated using built-in functionality of RATs.")

## ------------------------------------------------------------------------
# 1. Restructure sample_to_covariates as covariate to samples.
samples_by_covariate <- group_samples(myslo$sample_to_covariates)

# There are two covariates in the simulated data:
print( names(samples_by_covariate) )

# Covariate "condition" in our simulated data has two values:
print( names( samples_by_covariate$condition ) )

# 2. Extract the bootstrapped abundance tables for each of the two conditions:
condA_boots <- denest_sleuth_boots(myslo, myannot, 
                                   samples_by_covariate$condition[["controls"]])
condB_boots <- denest_sleuth_boots(myslo, myannot, 
                                   samples_by_covariate$condition[["patients"]])

# 3. Remove the sleuth object from memory.
# Make sure you have saved a copy of it to file before doing this!
rm(myslo)

# 4. Run RATs with the generic bootstrapped format:
mydtu <- call_DTU(annot= myannot, boot_data_A= condA_boots, boot_data_B= condB_boots, 
                  verbose= FALSE, 
                  name_A= "controls", name_B= "patients", varname = "condition",
                  description="Comparison of two sets of bootstrapped counts for the tutorial, extracted from a simulated sleuth object. Simulated using built-in functionality of RATs.")

## ---- eval=FALSE---------------------------------------------------------
#  # 1. Collect your outputs into vectors. The end of each path should be a
#  # directory with a unique name/identifier for one sample.
#  samples_A <- c("your/path/SAMPLE1", "your/path/SAMPLE4","your/path/SAMPLE5", etc)
#  samples_B <- c("your/path/SAMPLE2", "your/path/SAMPLE3","your/path/SAMPLE7", etc)
#  
#  # 2a. Convert, import, and extract from Salmon.
#  boots <- fish4rodents(A_paths= samples_A, B_paths= samples_B, half_cooked= FALSE)
#  
#  # 2b. OR if it is already in Kallisto format.
#  boots <- fish4rodents(A_paths= samples_A, B_paths= samples_B, half_cooked= TRUE)
#  
#  # 3. You might want to save boots to file, in case you want to re-run RATs on it
#  # later with different parameters or in case you make a mistake, etc...
#  
#  # 4. Run RATs with the generic bootstrapped format:
#  mydtu <- call_DTU(annot= myannot, boot_data_A= boots$boot_data_A, boot_data_B= boots$boot_data_B,
#                    verbose= FALSE,
#                    name_A= "controls", name_B= "patients", varname = "condition",
#                    description="Comparison of two sets of bootstrapped counts imported from the quantification output.")

## ----eval=FALSE----------------------------------------------------------
#  # Calling DTU with custom thresholds.
#  mydtu <- call_DTU(annot= myannot, slo= myslo, name_A= "controls", name_B= "patients",
#                    p_thresh= 0.01, abund_thresh= 10, dprop_thres = 0.25)

## ----eval=FALSE----------------------------------------------------------
#  # Bootstrap (default). Do 100 iterations.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", qboot = TRUE, qbootnum = 100, qrep_thresh= 0.95)
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", qboot = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Bootstrap (default).
#  # Note that for few replicates, the reproducibility values are highly discrete.
#  # NOTE: The number of iterations for 3 samples per condition is 3*3=9.
#  # So the minimum possible error rate is 1/9=0.111.
#  # The corresponding threshold is 1-0.111=0.888.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", rboot = TRUE, qrep_thresh= 0.85)
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", rboot = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Using 8 threads/cores for parallel computing.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", threads = 8)

## ----eval=FALSE----------------------------------------------------------
#  # Transcripts only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", testmode="transc")
#  # Genes only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", testmode="genes")

## ----eval=FALSE----------------------------------------------------------
#  # Bonferroni correction.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", correction = "bonferroni")

## ------------------------------------------------------------------------
# Lets emulate some input with custom field names. 
sim <- sim_sleuth_data(varname="mouse", cnames=c("Splinter", "Mickey"), 
                       COUNTS_COL="the-counts", TARGET_COL="transcript", 
                       PARENT_COL="gene", BS_TARGET_COL = "trscr")
myslo <- sim$slo
myannot <- sim$annot

# Sleuth covariates table.
print( sim$slo$sample_to_covariates )
# Sleuth bootstrapped quantifications.
print( head(sim$slo$kal[[1]]$bootstrap[[1]]) )
# Annotation.
print( head(sim$annot) )

## ---- eval=FALSE---------------------------------------------------------
#  # Call DTU on data with custom field names.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey",
#                    varname="mouse", COUNTS_COL="the-counts", BS_TARGET_COL="trscr",
#                    verbose = FALSE)
#  
#  #! In our example data here, the annotation fields have been modified as well,
#  #! so this command will NOT run as it is shown. You will also need the parameters
#  #! shown in teh section below.

## ---- eval=FALSE---------------------------------------------------------
#  # Call DTU using annotation with custom field names.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey",
#                    varname="mouse", TARGET_COL="transcript", PARENT_COL="gene", verbose = FALSE)
#  
#  #! In our example data here, the Sleuth fields have been modified as well,
#  #! so this command will NOT run as it is shown. You will also need the parameters
#  #! shown in the section above.

