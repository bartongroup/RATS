## ----eval=FALSE----------------------------------------------------------
#  # 1. Build latest developmental version from Github:
#  devtools::install_github("bartongroup/rats")
#  
#  # 2. Load into R session.
#  library{rats}
#  
#  # 3. Specify transcript grouping:
#  my_identifiers_table <- annot2ids("my_annotation.gtf")
#  
#  # 4a. Call DTU on a sleuth object, using default settings:
#  mydtu <- call_DTU(annot= my_identifiers_table, slo= my_sleuth_object,
#                    name_A= "My_condition", name_B= "My_other_condition")
#  # 4b. Call DTU on generic bootstrapped abundance estimates:
#  mydtu <- call_DTU(annot= my_identifiers_table, boot_data_A= my_list_data_tables_A,
#                    boot_data_B= my_list_data_tables_A)
#  # 4c. Call DTU on generic abundance estimates:
#  mydtu <- call_DTU(annot= my_identifiers_table, count_data_A= my_data_table_A,
#                    count_data_B= my_data_table_B, boots= "none")
#  
#  # 5. Get all gene and transcript identifiers per category (significant DTU,
#  # no DTU, Not Applicable):
#  myids <- get_dtu_ids(mydtu)
#  
#  # 6. Plot significance VS effect size:
#  plot_overview(mydtu)

## ------------------------------------------------------------------------
library(rats)

## ------------------------------------------------------------------------
# Show the first rows of the table corresponding to one sample, from simulated data.
head(sim_boot_data()[[2]][[1]])

## ------------------------------------------------------------------------
# Show the first rows of the table corresponding to the annotation, from simulated data.
head(sim_count_data()[[1]])

## ----eval=FALSE----------------------------------------------------------
#  # Extract transcript ID to gene ID index from a GTF annotation.
#  myannot <- annot2ids("my_annotation_file.gtf")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_count_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]          # Simulated abundances for one condition.
mycond_B <- simdat[[3]]          # Simulated abundances for other condition.
myannot <- simdat[[1]]   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# Find DTU between the simulated datasets.
mydtu <- call_DTU(annot= myannot, count_data_A= mycond_A, count_data_B= mycond_B, 
                  boots= "none", verbose= FALSE,
                  name_A= "healthy", name_B= "patients", varname= "My phenotype",
                  description="Comparison of two simulated counts datasets for the
                               tutorial. Simulated using built-in functionality of 
                               `rats`.")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_boot_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]          # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]          # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot= myannot, boot_data_A= mycond_A, boot_data_B= mycond_B, 
                  name_A= "wildtype", name_B= "some mutant", varname = "My phenotype",
                  verbose= FALSE, description="Comparison of two simulated datasets 
                    of bootstrapped counts for the tutorial. Simulated using built-in 
                    functionality of `rats`.")

## ------------------------------------------------------------------------
# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# controls and patients are arbitrary names to use as conditions.

# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo       # Simulated minimal sleuth object.
myannot <- simdat$annot   # Transcript and gene Identifiers for the above data.

## ------------------------------------------------------------------------
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Using a simulated sleuth object for the purposes of the tutorial.
                               Simulated using built-in functionality of `rats`.")

## ------------------------------------------------------------------------
# See available variables and values.
print( myslo$sample_to_covariates )

## ----eval=FALSE----------------------------------------------------------
#  # Compare samples by a non-default variable.
#  mydtu <- call_DTU(annot= myannot, slo= myslo, name_A= "ba", name_B= "bb", varname= "batch")

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
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients",
#                    p_thresh = 0.01, count_thresh = 10, dprop_thresh = 0.25)

## ------------------------------------------------------------------------
# Compare by a different variable. In this case "batch".
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "ba", name_B = "bb", 
                  varname= "batch", verbose = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Bootstrap both types of DTU calls (default), for 100 iterations (default).
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", boots = "both", bootnum = 100)
#  
#  # Only bootstrap transcript calls.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", boots = "transc")
#  
#  # Only bootstrap gene calls.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", boots = "genes")
#  
#  # Skip bootstraps.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", boots = "none")

## ----eval=FALSE----------------------------------------------------------
#  # Transcripts only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", testmode="transc")
#  # Genes only.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", testmode="genes")

## ----eval=FALSE----------------------------------------------------------
#  # Bonferroni correction.
#  mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls",
#                    name_B = "patients", correction = "bonferroni")

## ------------------------------------------------------------------------
# Lets emulate some input with custom field names. 
sim <- sim_sleuth_data(varname="mouse", cnames=c("Splinter", "Mickey"), 
                       COUNTS_COL="the-counts", TARGET_COL="transcript", 
                       PARENT_COL="gene", BS_TARGET_COL = "trscr")
myslo <- sim$slo
myannot <- sim$annot

print( sim$slo$sample_to_covariates )
print( head(sim$slo$kal[[1]]$bootstrap[[1]]) )
print( head(sim$annot) )

## ------------------------------------------------------------------------
# Call DTU on data with custom field names.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey", 
                  varname="mouse", TARGET_COL="transcript", PARENT_COL="gene", 
                  COUNTS_COL="the-counts", BS_TARGET_COL="trscr", verbose = FALSE)

# The output structure will always use the same field names, regardless of 
# what the input field names are.
print( names(mydtu$Transcripts) )

