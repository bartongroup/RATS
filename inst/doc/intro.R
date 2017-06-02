## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval= FALSE--------------------------------------------------------
#  install.packages(c("data.table", "matrixStats"), dependencies=TRUE)

## ---- eval= FALSE--------------------------------------------------------
#  install.packages("ggplot2", dependencies=TRUE)

## ---- eval= FALSE--------------------------------------------------------
#  # Devtools (available on both CRAN and Bioconductor), needed by wasabi apparently?
#  install.packages("devtools", dependencies=TRUE)
#  
#  source("http://bioconductor.org/biocLite.R")
#  # Wasabi converter from Salmon/Sailfish to Kallisto.
#  biocLite("COMBINE-lab/wasabi")
#  #Kallisto parser
#  biocLite("rhdf5")

## ---- eval= FALSE--------------------------------------------------------
#  install.packages("shiny", dependencies=TRUE)

## ---- eval= FALSE--------------------------------------------------------
#  install.packages("<path/to/dowloaded/package>", repos = NULL, type="source")

## ---- eval= FALSE--------------------------------------------------------
#  devtools::install_github("bartongroup/rats", ref="master")

## ---- eval= FALSE--------------------------------------------------------
#  devtools::install_github("bartongroup/rats", ref="development")

## ---- eval= FALSE--------------------------------------------------------
#  # 1. Load into R session.
#  library(rats)
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
#                    count_data_B= my_data_table_B, qboot= FALSE)
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

