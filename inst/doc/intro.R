## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=FALSE)

## ------------------------------------------------------------------------
#  install.packages(c("data.table", "matrixStats"), dependencies=TRUE)

## ------------------------------------------------------------------------
#  install.packages("ggplot2", dependencies=TRUE)

## ------------------------------------------------------------------------
#  # Devtools (available on CRAN), needed by the other two.
#  install.packages("devtools", dependencies=TRUE)
#  
#  # Sleuth
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("rhdf5")
#  devtools::install_github("pachterlab/sleuth")
#  
#  # Wasabi
#  biocLite("COMBINE-lab/wasabi")

## ------------------------------------------------------------------------
#  install.packages("shiny", dependencies=TRUE)

## ------------------------------------------------------------------------
#  install.packages("<path/to/dowloaded/package>", repos = NULL, type="source")

## ------------------------------------------------------------------------
#  devtools::install_github("bartongroup/rats", ref="master")

## ------------------------------------------------------------------------
#  devtools::install_github("bartongroup/rats", ref="development")

## ------------------------------------------------------------------------
#  # 1. Load into R session.
#  library{rats}
#  
#  # 2. Specify transcript grouping:
#  my_identifiers_table <- annot2ids("my_annotation.gtf")
#  
#  # 3a. Call DTU on a sleuth object, using default settings:
#  mydtu <- call_DTU(annot= my_identifiers_table, slo= my_sleuth_object,
#                    name_A= "My_condition", name_B= "My_other_condition")
#  # 3b. Call DTU on generic bootstrapped abundance estimates:
#  mydtu <- call_DTU(annot= my_identifiers_table, boot_data_A= my_list_data_tables_A,
#                    boot_data_B= my_list_data_tables_A)
#  # 3c. Call DTU on generic abundance estimates:
#  mydtu <- call_DTU(annot= my_identifiers_table, count_data_A= my_data_table_A,
#                    count_data_B= my_data_table_B, qboots= FALSE)
#  
#  # 4. Plot significance VS effect size:
#  plot_overview(mydtu)
#  
#  # 5a. Get all gene and transcript identifiers per category
#  # (significant DTU, no DTU, Not Applicable):
#  myids <- get_dtu_ids(mydtu)
#  
#  # 5b. Get all gene and transcript identifiers implicated in isoform switching:
#  myids <- get_switch_ids(mydtu)
#  
#  # 6. Plot isoform changes for a given gene.
#  plot_gene(mydtu, "my_awesome_gene_ID")

