# Contents:

1 Analysis Scripts
  - commands.txt : commands and parameters used to achieve the results presented.
  - dexseq_dtu.R : DEXSeq DTU for the simulated datasets
  - drimseq_human.R : DRIMSeq DTU for the Deng et al 2013 dataset
  - fileutilities.py : General housekeeping and helper script used to avoid error-prone copy-paste-editing of commands
  - gen_extended_data_1.Rmd : Analysis worksheet for the A.thaliana dataset DTU results. Produces Extended_Data_1.html
  - gen_extended_data_2.Rmd : Analysis worksheet for the simulated datasets DTU results. Produces Extended_Data_2.html
  - gen_extended_data_5.Rmd : Analysis worksheet for the Deng et al 2103 dataset DTU results. Produces Extended_Data_5.html
  - mylogs.py : Command-logging functions. Irrelevant to the analysis, but required by fileutilities.py.
  - ratsDTU_FDR_count-FP.R : Count False Positives from the A.thaliana dataset DTU bootstraps.
  - ratsDTU_FDR_single.R : Execute a RATs DTU call for a single iteration of the A.thaliana dataset bootstraps.
  - ratsDTU-human_60.R : RATs DTU for the Deng et al 2013 dataset, using the Ensembl 60 assembly/annotation
  - ratsDTU-human_87.R : RATs DTU for the Deng et al 2013 dataset, using the Ensembl 87 assembly/annotation
  - ratsDTU_sonesonDm70-k_.R : RATs DTU for the simulated fruitfly dataset quantified with Kallisto.
  - ratsDTU_sonesonDm70-s_.R : RATs DTU for the simulated fruitfly dataset quantified with Salmon.
  - ratsDTU_sonesonHs71-k_.R : RATs DTU for the simulated human dataset quantified with Kallisto.
  - ratsDTU_sonesonHs71-s_.R : RATs DTU for the simulated human dataset quantified with Salmon.
  - scale_table.R : Scale TPM abundance values in a text table by a given factor.
  - sequtilities.py : Used to generate transcript-gene look-up pairs from a GTF annotation file.
2 Extended Data
  - Extended_Data_1.html : Analysis of A.thaliana dataset results.
  - Extended_Data_2.html : Analysis of simulated datasets results.
  - Extended_Data_3.txt : RATs DTU genes for the Deng et al 2013 dataset, using Ensembl 60 annotation.
  - Extended_Data_4.txt : RATs DTU genes for the Deng et al 2013 dataset, using Ensembl 87 annotation.
  - Extended_Data_5.html : Analysis of the Deng et al 2013 dataset results.
  - Extended_Data_6.docx : Mapping of Deng et al 2013 qRT-PCR primers to Ensembl 87 assembly/annotation.
3 f1000researchpreprint.pdf : This manuscript
4 Figure1.tiff : manuscript figure 1
5 Figure2.tiff : manuscript figure 2
6 Figure3.tiff : manuscript figure 3
7 Figure4.tiff : manuscript figure 4
8 Figure5.tiff : manuscript figure 5
