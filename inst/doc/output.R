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

# DTU positive genes.
print( ids[["dtu-genes"]] )

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

## ------------------------------------------------------------------------
# Proportion and count changes for all the transcripts of the "MIX6" gene.
plot_gene(mydtu, "MIX6", style="lines")  # default

## ----eval=FALSE----------------------------------------------------------
#  plot_gene(mydtu, "MIX6", style="points")
#  plot_gene(mydtu, "MIX6", style="rainbow")
#  plot_gene(mydtu, "MIX6", style="merged")
#  plot_gene(mydtu, "MIX6", style="dashed")

## ----eval=FALSE----------------------------------------------------------
#  # Proportion change VS significance.
#  plot_overview(mydtu, type="volcano")

## ----eval=FALSE----------------------------------------------------------
#  # Distribution of maximum proportion change.
#  plot_overview(mydtu, type="maxdprop")

## ----eval=FALSE----------------------------------------------------------
#  # Start the interactive volcano plot.
#  plot_shiny_volcano(mydtu)

## ------------------------------------------------------------------------
library(ggplot2)

myplot <- plot_overview(mydtu, "volcano")
myplot  # display

# Change title. 
myplot2 <- myplot + ggtitle("My epic title")
myplot2

