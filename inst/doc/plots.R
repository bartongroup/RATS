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
mydtu <- call_DTU(annot= myannot, slo= myslo, name_A= "controls", name_B= "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object for the purposes of the tutorial. Simulated using built-in functionality of RATs.")

## ------------------------------------------------------------------------
# Grouping by condition (DEFAULT):
#   plot_gene(mydtu, "MIX6")
plot_gene(mydtu, "MIX6", style="bycondition")

## ------------------------------------------------------------------------
# Grouping by condition (minimalist):
plot_gene(mydtu, "MIX6", style="lines")

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

