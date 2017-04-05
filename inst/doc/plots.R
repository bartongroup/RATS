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
# Grouping by condition (DEAFULT):
#   plot_gene(mydtu, "MIX6")
plot_gene(mydtu, "MIX6", style="bycondition")

## ------------------------------------------------------------------------
# Grouping by condition (minimalist):
plot_gene(mydtu, "MIX6", style="linesonly")

## ------------------------------------------------------------------------
# Grouping by isoform:
plot_gene(mydtu, "MIX6", style="byisoform")

## ------------------------------------------------------------------------
# Change the encoded information.
plot_gene(mydtu, "MIX6", style="bycondition", fillby="isoform")
plot_gene(mydtu, "MIX6", style="byisoform", colourby="DTU", shapeby="replicate")

# For a less busy look, any of the information layers can be disabled.
plot_gene(mydtu, "MIX6", style="byisoform", colourby="none", shapeby="none")

## ------------------------------------------------------------------------
# Colour codes can be customised by specifying new values for
# condcolvec, replcolvec, isofcolvec, dtucolvec and nonecol.
plot_gene(mydtu, "MIX6", style="bycondition", fillby="condition", condcolvec=c("magenta", "cyan"))

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

