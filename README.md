# Relative Abundance of Transcripts (rats)
[![DOI](https://zenodo.org/badge/55973542.svg)](https://zenodo.org/badge/latestdoi/55973542)

Publication: https://doi.org/10.12688/f1000research.17916.1
***


## Description


#### Who it is for

Anyone working in transcriptomics, analysing gene expression and transcript abundances.


#### What it does

1. It provides a method to detect changes in the abundance ratios of transcript isoforms within a gene.
This is called **Differential Transcript Usage (DTU)**. 

2. **RATs is workflow-agnostic**. Quantification quality details are left to the quantification tools; 
RATs uses only the transcript abundances, which you can obtain using any tool you like. This makes it 
*suitable for use with alignment-free quantification tools* like [Kallisto](http://pachterlab.github.io/kallisto/)
or [Salmon](https://github.com/COMBINE-lab/salmon), as well as with traditional alignment-based quantification methods.

3. RATs is able to take advantage of the bootstrapped quantifications provided by the alignment-free tools. These bootstrapped
data are used by RATs to assess how much the technical variability of the heuristic quantifications affects differential transcript usage
and thus provide a measure of confidence in the DTU calls. 


#### What it needs

1. This is an R source package, and will run on any platform with a reasonably up-to-date R environment. A few third-party R packages are also required (see below).

2. As input, RATs requires transcript abundance estimates with or without bootstrapping. The format either way is tables with the samples (or iterations) as columns and the transcripts as rows. The first column holds the transcript IDs. Some functionality to create these from Salmon or Kallisto quantification files is provided by RATs, but there are many other ways to set your data in the required format.

3. RATs also requires a look-up table matching the transcript identifiers to the respective gene identifiers. This can be obtained through various means, such as extracting this info from a GTF file using functionality provided by RATs.

***

## How to use RATs

### Documentation

The package comes with full documentation: The vignettes are available locally by calling `browseVignettes("rats")` after the package is installed.
We recommend studying the vignettes before using RATs.


### Dependencies

The following instructions assume Bioconductor >=3.5 syntax. Consult [Bioconductor](https://bioconductor.org/install/) for the old syntax.

* Mandatory dependencies

```
# From CRAN

install.packages(c("data.table", "matrixStats", "ggplot2"))

# From Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicRanges")
```

* Optional dependencies

```
# Gene models
BiocManager::install("ggbio")

# Interactive volcano plot.
install.packages("shiny")
```


### Installation

The recommended option is to obtain the latest platform-independent package release from the [releases section](https://github.com/bartongroup/Rats/releases) here on **Github**.
Download the package file and then install it using:

`install.packages("<path/to/downloaded/package>", repos = NULL, type="source")`

The latest release can also be installed directly from Github, using the `devtools` package:

`devtools::install_github("bartongroup/rats", ref="master")`

For testing purposes (bug resolutions, new features), you can install the on-going development version from Github:

`devtools::install_github("bartongroup/rats", ref="development")`

Development versions are not archived and are **not suitable** for reproducible/publishable analyses. They may also temporarily not work correctly or not work at all. For important analyses you should use the latest release version.

Eventually, we aim to make the package also available through Bioconductor.


### Differential Transcript Usage

A minimal command to call DTU from bootstrapped quantifications looks like this:

`results <- call_DTU(annot= my_identifiers_table, boot_data_A= my_list_of_tables_A, boot_data_B= my_list_of_tables_B)`

and for plain quantifications like this:

`results <- call_DTU(annot= my_identifiers_table, count_data_A= my_table_A, count_data_B= my_table_B)`

In most cases, additional parameters will be needed. Please consult the input vignette before running RATs, to find the appropriate usage case for your needs.

***

## Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To report **problems** or ask for **assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *small reproducible working example* that demonstrates your issue is strongly encouraged. 
Also, be sure to **read the vignette(s)**, and browse/search the support forum before posting a new issue.

Enjoy!

![](./vignettes/figs/rats_logo.png)


