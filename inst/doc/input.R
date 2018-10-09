## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE------------------------------------------------------
library(rats)

## ---- eval=FALSE---------------------------------------------------------
#  library(rats)

## ---- echo=FALSE---------------------------------------------------------
# Show the first rows of the table corresponding to one sample, from emulated data.
head(sim_boot_data()[[2]][[1]])

## ---- echo=FALSE---------------------------------------------------------
# Show the first rows of the table corresponding to the annotation, from simulated data.
head(sim_count_data()[[1]])

## ---- eval=FALSE---------------------------------------------------------
#  # Extract transcript ID to gene ID index from a GTF annotation.
#  myannot <- annot2ids("my_annotation_file.gtf")

## ---- eval=FALSE---------------------------------------------------------
#  # Simulate some data.
#  simdat <- sim_count_data(clean=TRUE)
#  # For convenience let's assign the contents of the list to separate variables
#  mycond_A <- simdat[[2]]       # Simulated abundances for one condition.
#  mycond_B <- simdat[[3]]       # Simulated abundances for other condition.
#  myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.

## ---- eval=FALSE---------------------------------------------------------
#  # Find DTU between the simulated datasets.
#  mydtu <- call_DTU(annot= myannot, count_data_A= mycond_A,
#                    count_data_B= mycond_B, verbose= FALSE,
#                    name_A= "healthy", name_B= "patients",
#                    varname= "My phenotype",
#                    description="Comparison of two simulated counts datasets
#                    for the tutorial. Simulated using built-in functionality
#                    of RATs.")

## ------------------------------------------------------------------------
# Simulate some data. (Notice it is a different function than before.)
simdat <- sim_boot_data(clean=TRUE)

# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

## ---- eval=FALSE---------------------------------------------------------
#  # Find DTU between conditions.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mycond_A,
#                    boot_data_B= mycond_B, verbose= FALSE,
#                    name_A= "wildtype", name_B= "some mutant",
#                    varname = "My phenotype", description="Comparison of
#                    two simulated datasets of bootstrapped counts for the
#                    tutorial. Simulated using built-in functionality
#                    of RATs.")

## ---- eval=FALSE---------------------------------------------------------
#  # 1. Collect your outputs into vectors. The end of each path should be a
#  #    directory with a unique name/identifier for one sample.
#  samples_A <- c("your/path/SAMPLE1", "your/path/SAMPLE4","your/path/SAMPLE5")
#  samples_B <- c("your/path/SAMPLE2", "your/path/SAMPLE3","your/path/SAMPLE7",
#                 "your/path/SAMPLE10")
#  
#  # 2. Calculate length- & library-normalised abundances.
#  #    Scale them to 1M reads for TPM values.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot, scaleto=100000000)

## ---- eval=FALSE---------------------------------------------------------
#  # 3. Run RATs with the bootstrapped table data format.
#  #    Scale the TPM abundances to the respective library sizes.
#  #    Scaling factors are: (target size) / (current size)
#  #    Here the current size is 1M (ie. TPMs), as per the fish4rodents() step.
#  #    For the example, assume library sizes of 25M, 26M, 23M, 50M, 45M,
#  #    48M and 52M for the 7 samples, in the same order.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B, verbose= FALSE,
#                    scaling=c(25, 26, 23, 50, 45, 48, 52))

## ---- eval=FALSE---------------------------------------------------------
#  # Calling DTU with custom thresholds.
#  mydtu <- call_DTU(annot= myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    p_thresh= 0.01, dprop_thres = 0.15, abund_thresh= 10)

## ---- eval=FALSE---------------------------------------------------------
#  # Bootstrap (default). Do 100 iterations.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    qboot = TRUE, qbootnum = 100, qrep_thresh= 0.95)
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    qboot = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  # Bootstrap (default).
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    rboot = TRUE, qrep_thresh= 0.85)
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    rboot = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  # Lean run (default).
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    lean =TRUE)
#  
#  # Extra info on variance across iterations.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    lean = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Using 8 threads/cores for parallel computing.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    threads = 8)

## ---- eval=FALSE---------------------------------------------------------
#  # Transcripts only.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    testmode="transc")
#  # Genes only.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    testmode="genes")

## ---- eval=FALSE---------------------------------------------------------
#  # Bonferroni correction.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    correction = "bonferroni")

## ---- eval=FALSE---------------------------------------------------------
#  # Call DTU using annotation with custom field names.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    TARGET_COL="transcript", PARENT_COL="gene")

## ---- eval=FALSE---------------------------------------------------------
#  # The following are equivalent.
#  
#  # 1:
#  # Scale directly to library sizes at the import step.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto=c(25123456, 2665431, 23131313,
#                                   5000000, 45123132, 48456654, 52363636))
#  # No additional scaling needed.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling=1)  # default
#  
#  # 2:
#  # Normalise quantifications but do not scale them.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto=1)
#  # Scale directly to the smallest library size at the run step.
#  libsiz <- min(25123456, 2665431, 23131313, 5000000, 45123132, 48456654, 52363636)
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling=libsiz)
#  
#  # 3:
#  # Scale Kallisto/Salmon quantifications to TPMs.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto=10000000)  # default
#  
#  # Scale TPMs to actual library sizes.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling=c(25.123456, 26.65431, 23.131313, 50.0, 45.123132, 48.456654, 52.363636))

## ---- eval=FALSE---------------------------------------------------------
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    reckless=TRUE, verbose=TRUE)

