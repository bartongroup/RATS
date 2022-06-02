## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(rats)
library(data.table)

## -----------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_boot_data(clean=TRUE) 
# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

# Call DTU
mydtu <- call_DTU(annot= myannot, verbose= FALSE,
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  dprop_thresh=0.1, qboot=TRUE, rboot=FALSE)

## -----------------------------------------------------------------------------
print( names(mydtu) )

## -----------------------------------------------------------------------------
# Parameter list's elements.
print( names(mydtu$Parameters) )

## -----------------------------------------------------------------------------
# Genes table's fields.
print( names(mydtu$Genes) )

## -----------------------------------------------------------------------------
# Transcripts table's fields.
print( names(mydtu$Transcripts) )

## -----------------------------------------------------------------------------
# Elements of Abundances.
print( names(mydtu$Abundances) )

## -----------------------------------------------------------------------------
# Abundance table for first condition.
print( head(mydtu$Abundances[[1]]) )

## -----------------------------------------------------------------------------
# A tally of the outcome.
print( dtu_summary(mydtu) )

## -----------------------------------------------------------------------------
# Gene and transcript IDs corresponding to the tally above.
ids <- get_dtu_ids(mydtu)
print( names(ids) )
print( ids )

## -----------------------------------------------------------------------------
# A tally of genes switching isoform ranks.
print( dtu_switch_summary(mydtu) )

# The gene IDs displaying isoform switching.
ids <- get_switch_ids(mydtu)
print( names(ids) )

## -----------------------------------------------------------------------------
# A tally of genes switching isoform ranks.
print( dtu_plurality_summary(mydtu) )

# The gene IDs displaying isoform switching.
ids <- get_plurality_ids(mydtu)

