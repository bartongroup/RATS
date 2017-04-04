## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(rats)

# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo
myannot <- simdat$annot

# Call DTU
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object 
                    for the purposes of the tutorial. Simulated using built-in functionality 
                    of RATs.")

## ------------------------------------------------------------------------
# A really simple tally of the outcome.
print( dtu_summary(mydtu) )

## ------------------------------------------------------------------------
# Gene and transcript IDs corresponding to the tally above.
ids <- get_dtu_ids(mydtu)

# Contents
print( names(ids) )

# DTU positive genes according to the transcript-level test.
print( ids[[4]] )


## ------------------------------------------------------------------------
print( names(mydtu) )

## ------------------------------------------------------------------------
# Parameter list's elements.
print( names(mydtu$Parameters) )

## ------------------------------------------------------------------------
# Genes table's fields.
print( names(mydtu$Genes) )

## ------------------------------------------------------------------------
# Transcripts table's fields.
print( names(mydtu$Transcripts) )

## ------------------------------------------------------------------------
# Elements of ReplicateData
print( names(mydtu$Abundances) )

