## ----eval=FALSE----------------------------------------------------------
#  # 1. Build latest developmental version from Github:
#  devtools::install_github("bartongroup/rats")
#  
#  # 2. Load into R session.
#  library{rats}
#  
#  # 3. Call DTU on a sleuth object, using default settings.
#  mydtu <- call_DTU(annot = my_identifiers_table, slo = my_sleuth_object, name_A = "My_condition", name_B = "My_other_condition")
#  
#  # 4. Tally of results.
#  dtu_summary(mydtu)
#  
#  # 5. Get all gene and transcript identifiers that correspond to the above tally.
#  myids <- get_dtu_ids(mydtu)
#  # A list of vectors, one for each category of the tally.
#  
#  # 6. Plot significance VS effect size.
#  plot_overview(mydtu)

## ------------------------------------------------------------------------
library(rats)

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("foo", "bar")) # foo and bar are arbitrary names 
                                                    # to use as conditions.

# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo       # Simulated minimal sleuth object.
myannot <- simdat$annot   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# This table is important. It assigns samples to variables. We can only compare data based on the 
# variables and values listed in this table.
# This one has two variables: "condition" and "batch". Yours may have more/fewer/different ones.
# Notice that "foo" and "bar" (the names we gave to our simulated data) are under "condition".
print( myslo$sample_to_covariates )

# This is what the estimated counts tables look like. One such table per bootstrap, per sample.
# head() shows the first few rows only.
print( head(myslo$kal[[1]]$bootstrap[[1]]) )

## ------------------------------------------------------------------------
# This is what the annotation table should look like.
# head() shows the first few rows only.
print( head(myannot) )

## ------------------------------------------------------------------------
# Find DTU between conditions "foo" and "bar" in the simulated data.
# Be warned that the vignette format does not display progress bars properly!
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar")

## ----eval=FALSE----------------------------------------------------------
#  # Comparing samples by a different variable.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "ba", name_B = "bb", varname="batch")

## ------------------------------------------------------------------------
# This prints out nothing at all.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", verbose = FALSE)

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
# Parameter list's elements.
print( names(mydtu$Parameters) )

## ------------------------------------------------------------------------
# Genes table's fields.
print( names(mydtu$Genes) )

## ------------------------------------------------------------------------
# Transcripts table's fields.
print( names(mydtu$Transcripts) )

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
# Proportion changes for all the transcripts of the "MIX6" gene.
plot_gene(mydtu, "MIX6", vals="proportions")

## ------------------------------------------------------------------------
# Absolute expression changes for all the transcripts of the "MIX6" gene.
# The ERROR BARS represent 2 standard deviations from the mean count across replicates.
plot_gene(mydtu, "MIX6", vals="counts")

## ------------------------------------------------------------------------
# Proportion change VS significance.
plot_overview(mydtu, type="volcano")

## ------------------------------------------------------------------------
# Distribution of maximum proportion change.
plot_overview(mydtu, type="maxdprop")

## ------------------------------------------------------------------------
# Transcript-level confidence threshold VS. number of DTU positive calls.
plot_overview(mydtu, type="transc_conf")

## ------------------------------------------------------------------------
# Gene-level confidence threshold VS. number of DTU positive calls.
plot_overview(mydtu, type="gene_conf")

## ------------------------------------------------------------------------
library(ggplot2)

myplot <- plot_overview(mydtu, "volcano")
myplot  # display

# Change scale of y axis to linear. 
myplot2 <- myplot + scale_y_continuous(trans = "identity")
myplot2

## ----eval=FALSE----------------------------------------------------------
#  # Calling DTU with custom thresholds.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", p_thresh = 0.01,
#                    count_thresh = 10, dprop_thresh = 0.25)

## ------------------------------------------------------------------------
# Compare by a different variable. In this case "batch".
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "ba", name_B = "bb", varname= "batch", verbose = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Bootstrap everything. (default)
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", boots = "both", bootnum = 100)
#  
#  # Only bootstrap transcript calls.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", boots = "transc")
#  
#  # Only bootstrap gene calls.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", boots = "genes")
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", boots = "none")

## ----eval=FALSE----------------------------------------------------------
#  # Transcripts only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", testmode="transc")
#  # Genes only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", testmode="genes")

## ----eval=FALSE----------------------------------------------------------
#  # Bonferroni correction.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "foo", name_B = "bar", correction = "bonferroni")

## ------------------------------------------------------------------------
# Lets create some input with custom field names. The data is exactly the same as before.
sim <- sim_sleuth_data(varname="mouse", cnames=c("Splinter", "Mickey"), COUNTS_COL="the-counts", 
                       TARGET_COL="transcript", PARENT_COL="gene", BS_TARGET_COL = "trscr")
myslo <- sim$slo
myannot <- sim$annot

print( sim$slo$sample_to_covariates )
print( head(sim$slo$kal[[1]]$bootstrap[[1]]) )
print( head(sim$annot) )

## ------------------------------------------------------------------------
# Call DTU on data with custom field names.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey", varname="mouse", 
                  TARGET_COL="transcript", PARENT_COL="gene", 
                  COUNTS_COL="the-counts", BS_TARGET_COL="trscr", verbose = FALSE)

# The output structure will always use the same field names, regardless of 
# what the input field names are.
print( names(mydtu$Transcripts) )

