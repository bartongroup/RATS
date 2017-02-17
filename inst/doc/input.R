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
                  description="Comparison of two simulated counts datasets for the
                    tutorial. Simulated using built-in functionality of RATs.")

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
                  description="Comparison of two simulated datasets of bootstrapped 
                    counts for the tutorial. Simulated using built-in functionality 
                    of RATs.")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# controls and patients are arbitrary names to use as conditions.

# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo       # Simulated minimal sleuth object.
myannot <- simdat$annot   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object 
                    for the purposes of the tutorial. Simulated using built-in functionality 
                    of RATs.")

## ----eval=FALSE----------------------------------------------------------
#  # Calling DTU with custom thresholds.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients",
#                    p_thresh = 0.01, abund_thresh = 10, dprop_thresh = 0.25)

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
#  # NOTE: The number of iterations for 3 samples per condition is 3*3=9.
#  # So the minimum possible error rate is 1/9=0.111. The corresponding threshold is 1-0.111=0.888.
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
#                    varname="mouse", COUNTS_COL="the-counts", BS_TARGET_COL="trscr", verbose = FALSE)
#  
#  #! In our example data here, the annotation fields have been modified as well, so this command
#  #! will NOT run as it is shown. You will also need the parameters shown in teh section below.

## ---- eval=FALSE---------------------------------------------------------
#  # Call DTU using annotation with custom field names.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey",
#                    varname="mouse", TARGET_COL="transcript", PARENT_COL="gene", verbose = FALSE)
#  
#  #! In our example data here, the Sleuth fields have been modified as well, so this command
#  #! will NOT run as it is shown. You will also need the parameters shown in the section above.

