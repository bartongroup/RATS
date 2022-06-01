## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE-----------------------------------------------------------
library(rats)

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
# Grouping by condition (DEFAULT):
plot_gene(mydtu, "MIX", style="bycondition")

## -----------------------------------------------------------------------------
# Grouping by isoform:
plot_gene(mydtu, "MIX", style="byisoform")

## ----eval=FALSE---------------------------------------------------------------
#  models <- annot2models('/my/annotation/file.gtf')
#  library(ggbio)
#  # This will plot the structure of all isoforms for the given gene ID.
#  autoplot(models[['mygeneID']])

## ---- eval=FALSE--------------------------------------------------------------
#  # Proportion change VS transcript-level significance. Each point is a transcript
#  plot_overview(mydtu, type="tvolcano")
#  
#  # This can also be plotted for genes, by using the largest isoform effect size as proxy.
#  plot_overview(mydtu, type="gvolcano")

## ---- eval=FALSE--------------------------------------------------------------
#  # Distribution of proportion change.
#  plot_overview(mydtu, type="dprop")
#  
#  # Distribution of largest isoform proportion change per gene.
#  plot_overview(mydtu, type="maxdprop")

## ---- eval=FALSE--------------------------------------------------------------
#  # Proportion change VS transcript-level significance. Each point is a transcript
#  plot_overview(mydtu, type="fcvolcano")
#  
#  # This can also be plotted for genes, by using the largest isoform effect size as proxy.
#  plot_overview(mydtu, type="fcVSdprop")

## -----------------------------------------------------------------------------
# Pairwise Pearson's correlations among samples.
plot_diagnostics(mydtu, type='cormat') # Default type.

## ---- eval=FALSE--------------------------------------------------------------
#  # Start the interactive volcano plot.
#  plot_shiny_volcano(mydtu)

