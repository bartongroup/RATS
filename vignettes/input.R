## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE-----------------------------------------------------------
library(rats)

## ---- echo=FALSE--------------------------------------------------------------
# Show the first rows of the table corresponding to one sample, from emulated data.
head(sim_boot_data()[[2]][[1]])

## ---- echo=FALSE--------------------------------------------------------------
# Show the first rows of the table corresponding to the annotation, from simulated data.
head(sim_count_data()[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  # Extract transcript ID to gene ID index from a GTF annotation.
#  myannot <- gtf2ids("my_annotation_file.gtf")
#  
#  # Extract transcript ID and gene ID from a GRanges object.
#  # It must have GTF-style metadata columns "gene_id" and "transcript_id".
#  myannot <- granges2ids(mygranges)

## ----eval=FALSE---------------------------------------------------------------
#  myannot <- select(mytxdb, keys(mytxdb), "GENEID", "TXNAME")
#  # Rename the columns to match what RATs expects.
#  names(myannot) <- c('gene_id', 'target_id')

## -----------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_count_data(clean=TRUE)
# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]       # Simulated abundances for one condition.
mycond_B <- simdat[[3]]       # Simulated abundances for other condition.
myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.

## -----------------------------------------------------------------------------
print( class(mycond_A) )

## -----------------------------------------------------------------------------
# Find DTU between the simulated datasets.
mydtu <- call_DTU(annot= myannot, count_data_A= mycond_A, 
                  count_data_B= mycond_B, verbose= FALSE,
                  scaling= 1,
                  name_A= "healthy", name_B= "patients", 
                  varname= "My phenotype",
                  description="Comparison of two simulated counts datasets 
                  for the tutorial. Simulated using built-in functionality 
                  of RATs.")

## -----------------------------------------------------------------------------
# Simulate some data. (Notice it is a different function than before.)
simdat <- sim_boot_data(clean=TRUE)

# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

## -----------------------------------------------------------------------------
print( class(mycond_A) )
print( class(mycond_A[[1]]) )

## -----------------------------------------------------------------------------
# Find DTU between conditions.
mydtu <- call_DTU(annot= myannot, boot_data_A= mycond_A, 
                  boot_data_B= mycond_B, verbose= FALSE, 
                  scaling= 1,
                  name_A= "wildtype", name_B= "some mutant", 
                  varname = "My phenotype", description="Comparison of 
                  two simulated datasets of bootstrapped counts for the 
                  tutorial. Simulated using built-in functionality 
                  of RATs.")

## ---- eval=FALSE--------------------------------------------------------------
#  # Mock-up code, does not run.
#  
#  # 1. Collect your outputs into vectors. The end of each path should be a
#  #    directory with a unique name/identifier for one sample and containing
#  #    the respective output of Kallisto.
#  samples_A <- c("your/path/SAMPLE1", "your/path/SAMPLE4","your/path/SAMPLE5")
#  samples_B <- c("your/path/SAMPLE2", "your/path/SAMPLE3","your/path/SAMPLE7", "your/path/SAMPLE10")
#  
#  # 2. Calculate length- & library-normalised abundances.
#  #    Scale them to 1M reads for TPM values.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot, scaleto=1e6) # scaling to TPM

## ---- echo=FALSE--------------------------------------------------------------
simdat <- sim_boot_data(clean=TRUE)
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

## -----------------------------------------------------------------------------
# Calling DTU with custom thresholds.
mydtu <- call_DTU(annot= myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  p_thresh= 0.01, dprop_thres = 0.15, abund_thresh= 10,
                  verbose= FALSE)

## -----------------------------------------------------------------------------
# Do bootstrap (default). Do 100 iterations.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  qboot = TRUE, qbootnum = 100, qrep_thresh= 0.95,
                  verbose= FALSE)

## -----------------------------------------------------------------------------
# Bootstrap (default).
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  rboot = TRUE, qrep_thresh= 0.85, verbose= FALSE)

## -----------------------------------------------------------------------------
# Extra info on variance across iterations.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  lean = FALSE, verbose= FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  # The following code is for demonstration only and won't run
#  # without valid paths supplied to fish4rodents().
#  
#  # The following are equivalent.
#  
#  # 1:
#  # Scale to individual library sizes directly at the import step.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto= c(25123456, 2665431, 23131313,
#                                   5000000, 45123132, 48456654, 52363636),
#                         verbose= FALSE)
#  # No additional scaling needed.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling= 1,  # default
#                    verbose= FALSE)
#  
#  # 2:
#  # Normalise quantifications but do not scale them at all.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto=1)
#  # Scale library fractions to the library sizes.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling= c(25123456, 2665431, 23131313, 5000000,
#                              45123132, 48456654, 52363636),
#                    verbose= FALSE)
#  
#  # 3:
#  # Scale Kallisto quantifications to TPMs.
#  mydata <- fish4rodents(A_paths= samples_A, B_paths= samples_B,
#                         annot= myannot,
#                         scaleto= 1000000)  # default
#  # Scale TPMs to individual library sizes.
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    scaling=c(25.123456, 26.65431, 23.131313, 50.0, 45.123132, 48.456654, 52.363636),
#                    verbose= FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  # Using 4 threads for parallel computing.
#  mydtu <- call_DTU(annot = myannot,
#                    boot_data_A= mycond_A, boot_data_B= mycond_B,
#                    threads = 4, verbose= FALSE)

## -----------------------------------------------------------------------------
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  use_sums = TRUE, verbose= FALSE)

## -----------------------------------------------------------------------------
# Bonferroni correction.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  correction = "bonferroni", verbose= FALSE)

## -----------------------------------------------------------------------------
# Transcripts only.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  testmode="transc", verbose= FALSE)
# Genes only.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  testmode="genes", verbose= FALSE)

## ---- echo=FALSE--------------------------------------------------------------
simdat <- sim_boot_data(clean=TRUE, PARENT_COL='gene', TARGET_COL='transcript')
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

## -----------------------------------------------------------------------------
# Call DTU using annotation with custom field names.
mydtu <- call_DTU(annot = myannot, 
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  TARGET_COL="transcript", PARENT_COL="gene",
                  verbose= FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  mydtu <- call_DTU(annot= myannot, boot_data_A= mydata$boot_data_A,
#                    boot_data_B= mydata$boot_data_B,
#                    reckless=TRUE, verbose=TRUE)

