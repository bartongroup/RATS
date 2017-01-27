## ----eval=FALSE----------------------------------------------------------
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
#                    count_data_B= my_data_table_B, boots= "none")
#  
#  # 4. Plot significance VS effect size:
#  plot_overview(mydtu)
#  
#  # 5. Get all gene and transcript identifiers per category (significant DTU,
#  # no DTU, Not Applicable):
#  myids <- get_dtu_ids(mydtu)
#  
#  # 6. Plot isoform changes for the top DTU gene returned in previous step.
#  plot_gene(mydtu, myids[["dtu-genes"]][1], style = "lines"")

