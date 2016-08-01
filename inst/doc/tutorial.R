## ------------------------------------------------------------------------
library(rats)
library(data.table)

## ------------------------------------------------------------------------
simdat <- sim_sleuth_data(cnames = c("foo", "bar")) # foo and bar are the names to use for the conditions.

# For convenience let's assign the contents of the list to separate variables.
slo <- simdat$slo       # Simulated minimal sleuth object.
annot <- simdat$annot  # Annotation for the above data.


## ------------------------------------------------------------------------
# Let's have a look at the simulated data.
summary(slo)

# Names of the samples and conditions.
print( slo$sample_to_covariates )

# Estimated fragment counts in bootstrap 1 of sample 1.
print( slo$kal[[1]]$bootstrap[[1]] )

## ------------------------------------------------------------------------
# And this is what the annotation looks like.
print(annot)

## ------------------------------------------------------------------------
mydtu <- call_DTU(slo, annot, "foo", "bar")

## ------------------------------------------------------------------------
mydtu <- call_DTU(slo, annot, "foo", "bar", verbose = FALSE)

## ------------------------------------------------------------------------
# A really simple tally of the outcome.
print( dtu_summary(mydtu) )

# Gene and transcript IDs corresponding to the tally above.
ids <- get_dtu_ids(mydtu)
print( ids[["dtu-genes"]] )  # Also "dtu-transc", "ndtu-genes", "ndtu-transc", "na-genes", "na-transc"

## ------------------------------------------------------------------------
print( summary(mydtu) )

## ------------------------------------------------------------------------
# See the elements.
print( names(mydtu$Parameters) )

# See the value of one of the elements (i.e. the P-value threshold).
print( mydtu$Parameters[["p_thresh"]] )

## ------------------------------------------------------------------------
# See the fields in the table.
print( names(mydtu$Genes) )

# See details of genes that show DTU.
print( mydtu$Genes[(DTU), ] )

# Or if you only want their codes:
print( mydtu$Genes[(DTU), parent_id] )

# NA are valid values for genes that didn't meet some pre-requisites:
# (more on this in the next two examples)
print( mydtu$Genes[, .(parent_id, elig, elig_fx, DTU)] )

## ------------------------------------------------------------------------
# See the fields in the table.
print( names(mydtu$Transcripts) )

# See the details of Transcripts that show relative change.
print( mydtu$Transcripts[(DTU), ] )

# Or if you only want their codes, without details:
print( mydtu$Transcripts[(DTU), .(target_id, parent_id)] )

# Not all transcripts meet the prerequisites:
# (more on this in the next two examples)
print( mydtu$Transcripts[, .(target_id, elig_xp, elig, elig_fx)] )

## ------------------------------------------------------------------------
print( mydtu$Genes["MIX6", ] )
print( mydtu$Transcripts["MIX6", ] )

# To search by any other criteria (like transcript code), you do it like this:
print( mydtu$Transcripts[target_id == "MIX6.c4", ] )

## ------------------------------------------------------------------------
library(data.table)
library(rats)

# Same data simulation as previous example.
sim <- sim_sleuth_data(cnames = c("foo", "bar"))

# Calling and bootstrapping DTU.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", boots = "both", bootnum = 100)

## ------------------------------------------------------------------------
# Types of list's elements.
print( summary(mydtu$Parameters) )

## ------------------------------------------------------------------------
# Types of table's elements.
print( lapply(mydtu$Genes, typeof) )

## ------------------------------------------------------------------------
# Types of table's elements.
print( lapply(mydtu$Transcripts, typeof) )

## ------------------------------------------------------------------------
# Let's check the info and settings.
print( mydtu$Parameters )

## ------------------------------------------------------------------------
# Gene-level calls.
print( mydtu$Genes )

## ------------------------------------------------------------------------
# Transcript-level calls.
print( mydtu$Transcripts )

## ------------------------------------------------------------------------
library(rats)
library(data.table)

sim <- sim_sleuth_data(cnames = c("foo", "bar"))

## ------------------------------------------------------------------------
# Compare by a different variable. In this case "batch".
mydtu <- call_DTU(sim$slo, sim$annot, "ba", "bb", varname= "batch", verbose = FALSE)

## ------------------------------------------------------------------------
# Bootstrap everything.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", boots = "both", bootnum = 100, verbose = FALSE)

# Only bootstrap transcript calls.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", boots = "transc", verbose = FALSE)

# Only bootstrap gene calls.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", boots = "genes", verbose = FALSE)

## ------------------------------------------------------------------------
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", p_thresh = 0.01, 
                  count_thresh = 10, dprop_thresh = 0.25, verbose = FALSE)
# You can view if/how this changes the result by using the commands shown in the 
# previous two examples.

## ------------------------------------------------------------------------
# Inspect higher statistical significance threshold.
new_thresh = 0.01
newgenecalls <- mydtu$Genes[, .(parent_id, 
                                pvalAB_corr < new_thresh  & 
                                  pvalBA_corr < new_thresh  & 
                                  elig_fx)]  # DTU requires both stats and biol significance
newtransccalls <- mydtu$Transcripts[, .(target_id, parent_id, 
                                        pval_corr < new_thresh  &  elig_fx)]  # same
print( newgenecalls )

## ------------------------------------------------------------------------
# Transcripts only.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", testmode="transc", verbose = FALSE)
# Genes only.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", testmode="genes", verbose = FALSE)

## ------------------------------------------------------------------------
# Bonferroni correction.
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", correction = "bonferroni", verbose = FALSE)

## ------------------------------------------------------------------------
# Lets create some input with custom field names. The data is exactly the same as before.
sim <- sim_sleuth_data(varname="mouse", cnames=c("Splinter", "Mickey"), COUNTS_COL="the-counts", 
                       TARGET_COL="transcript", PARENT_COL="gene", BS_TARGET_COL = "trscr")
print( sim$slo$sample_to_covariates )
print( sim$slo$kal[[1]]$bootstrap[[1]] )
print( sim$annot )

## ------------------------------------------------------------------------
mydtu <- call_DTU(sim$slo, sim$annot, "Splinter", "Mickey", varname="mouse", 
                  TARGET_COL="transcript", PARENT_COL="gene", 
                  COUNTS_COL="the-counts", BS_TARGET_COL="trscr", verbose = FALSE)

# The output structure will always use the same field names, regardless of 
# what the input names are.
print( names(mydtu$Transcripts) )

## ------------------------------------------------------------------------
library(rats)

sim <- sim_sleuth_data(cnames = c("foo", "bar"))
mydtu <- call_DTU(sim$slo, sim$annot, "foo", "bar", boots="both", bootnum=100, verbose = FALSE)

## ------------------------------------------------------------------------
# Proportion changes for all the transcripts of the "MIX6" gene.
# !!! The ERROR BARS represent the standard error for the estimated proportions.
plot_gene(mydtu, "MIX6")

## ------------------------------------------------------------------------
# Absolute expression changes for all the transcripts of the "MIX6" gene.
# !!! The ERROR BARS represent 2 standard deviations from the mean count across replicates.
plot_gene(mydtu, "MIX6", vals="counts")

## ------------------------------------------------------------------------
# Proportion change VS significance. Probably the most common plot. Similar to volcano plot.
plot_overview(mydtu, type="dpropVsig")

# Proportion change VS transcript expression.
plot_overview(mydtu, type="dpropVcount")

# Proportion VS transcript expression.
plot_overview(mydtu, type="propVcount")

# Distribution of maximum proportion change.
plot_overview(mydtu, type="maxdprop")


# Distribution of bootstrapped confidence of transcript-level DTU.
plot_overview(mydtu, type="transc-conf")

# Distribution of bootstrapped confidence of gene-level DTU.
plot_overview(mydtu, type="gene-conf")

# Transcript-level confidence threshold VS. number of DTU positive calls.
plot_overview(mydtu, type="trconfVdtu")

# Gene-level confidence threshold VS. number of DTU positive calls.
plot_overview(mydtu, type="gconfVdtu")

## ------------------------------------------------------------------------
library(ggplot2)

myplot <- plot_overview(mydtu, type="dpropVsig")
myplot  # display

# Change scale of y axis.
myplot2 <- myplot + coord_trans(y="sqrt")  # Or "log10", "log10", etc...
myplot2

# Revert to linear scale.
myplot2 <- myplot + coord_trans(y="identity")
myplot2

