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
# Split by condition for easier view of the overall gene profile.
plot_gene(mydtu, "MIX6", style="plain")

## ------------------------------------------------------------------------
# Isoforms paired up for easier individual comparisons.
plot_gene(mydtu, "MIX6", style="paired")

## ------------------------------------------------------------------------
# Split by condition.
plot_gene(mydtu, "MIX6", style="points")

## ------------------------------------------------------------------------
# Paired by isoform.
plot_gene(mydtu, "MIX6", style="pairedpnt")

## ------------------------------------------------------------------------
# Split by condition.
# This is the DEFAULT view if the style is omitted, as it is the most informative.
plot_gene(mydtu, "MIX6", style="lines")

## ------------------------------------------------------------------------
# A cleaner version, although it no longer shows which isoforms are DTU.
plot_gene(mydtu, "MIX6", style="linesonly")

## ------------------------------------------------------------------------
# You can change the information that is colour-coded.
plot_gene(mydtu, "MIX6", style="plain", fillby="DTU")
plot_gene(mydtu, "MIX6", style="points", fillby="isoform", colourby="replicate")
plot_gene(mydtu, "MIX6", style="pairedpnt", colourby="isoform", shapeby="replicate")

# For a less colourful look, the layered information can be disabled.
plot_gene(mydtu, "MIX6", style="points", fillby="none", colourby="none", shapeby="none")

## ------------------------------------------------------------------------
# You can also customise the colours used by specifying new values for
# condcolvec, replcolvec, isofcolvec, dtucolvec and nonecol.
plot_gene(mydtu, "MIX6", style="lines", fillby="condition", condcolvec=c("magenta", "cyan"))

## ------------------------------------------------------------------------
# Proportion change VS significance.
plot_overview(mydtu, type="volcano")

## ------------------------------------------------------------------------
# Distribution of maximum proportion change.
plot_overview(mydtu, type="maxdprop")

## ---- eval=FALSE---------------------------------------------------------
#  # Start the interactive volcano plot.
#  plot_shiny_volcano(mydtu)

## ------------------------------------------------------------------------
library(ggplot2)

myplot <- plot_overview(mydtu, "volcano")
myplot  # display

# Change title. 
myplot2 <- myplot + ggtitle("MY EPIC TITLE")
myplot2

