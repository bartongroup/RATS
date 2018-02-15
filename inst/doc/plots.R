## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE------------------------------------------------------
library(rats)

## ---- eval=FALSE---------------------------------------------------------
#  library(rats)

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_boot_data() 
# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

# Call DTU
mydtu <- call_DTU(annot= myannot, verbose= FALSE,
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  dprop_thresh=0.1, qboot=TRUE, rboot=FALSE)

## ------------------------------------------------------------------------
# Grouping by condition (DEFAULT):
plot_gene(mydtu, "MIX6", style="bycondition")

## ---- eval=FALSE---------------------------------------------------------
#  # Grouping by condition (minimalist):
#  plot_gene(mydtu, "MIX6", style="lines")

## ------------------------------------------------------------------------
# Grouping by isoform:
plot_gene(mydtu, "MIX6", style="byisoform")

## ---- eval=FALSE---------------------------------------------------------
#  # Proportion change VS transcript-level significance. Each point is a transcript
#  plot_overview(mydtu, type="tvolcano")
#  
#  # This can also be plotted for genes, by using the largest isoform effect size as proxy.
#  plot_overview(mydtu, type="gvolcano")

## ---- eval=FALSE---------------------------------------------------------
#  # Distribution of proportion change.
#  plot_overview(mydtu, type="dprop")
#  
#  # Distribution of largest isoform proportion change per gene.
#  plot_overview(mydtu, type="maxdprop")

## ---- eval=FALSE---------------------------------------------------------
#  # Proportion change VS transcript-level significance. Each point is a transcript
#  plot_overview(mydtu, type="fcvolcano")
#  
#  # This can also be plotted for genes, by using the largest isoform effect size as proxy.
#  plot_overview(mydtu, type="fcVSdprop")

## ------------------------------------------------------------------------
# Matrix of pairwise Pearson's correlations among samples.
plot_diagnostics(mydtu, type='cormat') # Default type.

## ---- eval=FALSE---------------------------------------------------------
#  # Start the interactive volcano plot.
#  plot_shiny_volcano(mydtu)

## ------------------------------------------------------------------------
# For a less busy look, any of the information layers can be disabled.
plot_gene(mydtu, "MIX6", style="byisoform", colourby="none", shapeby="none")

## ------------------------------------------------------------------------
plot_gene(mydtu, "MIX6", fillby="DTU", shapeby="none")

## ------------------------------------------------------------------------
# Colour codes can be customised by specifying new values to the
# corresponding parameters.
plot_gene(mydtu, "MIX6", style="bycondition", fillby="condition", 
          condcolvec=c("magenta", "cyan"),
          replcolvec=c("green", "darkgreen"))

