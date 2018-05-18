# Relative Abundance of Transcripts (rats)

***

## Description


#### Who it is for

Anyone working in transcriptomics, analysing gene expression and transcript abundances.


#### What it does

1. It provides a method to detect changes in the abundance ratios of transcript isoforms of a gene.
This is called **Differential Transcript Usage (DTU)**. 

2. **RATs is workflow-agnostic**. Quantification quality details are left to the quantification tools; RATs uses only the
transcript abundances, which you can obtain using any tool you like. This makes it *suitable for use with alignment-free quantification tools* 
like [Kallisto](http://pachterlab.github.io/kallisto/) or [Salmon](https://github.com/COMBINE-lab/salmon). 

3. RATs is able to take advantage of the bootstrapped quantifications provided by the alignment-free tools. These bootstrapped
data are used by RATs to assess how much the technical variability of the heuristic quantifications affects differential transcript usage
and thus provide a measure of confidence in the DTU calls. 


#### What it needs

1. This is an R source package, and will run on any platform with a reasonably up-to-date R environment. A few third-party R packages are also required (see below).

2. As input, RATs requires transcript abundance estimates with or without bootstrapping. The format either way is tables with the samples as columns and the transcripts as rows. An extra column holds the transcript IDs. Some functionality to create these from Salmon or Kallisto quantification files is provided by RATs.

3. RATs also requires a look-up table matching the transcript identifiers to the respective gene identifiers. This can be obtained through various means,
one of them being extracting this info from a GTF file using functionality provided by RATs.

***

## How to use RATs

### Documentation

The package comes with full documentation: The vignettes are available locally by calling `browseVignettes("rats")` after the package is installed.
We recommend studying the vignettes before using RATs.


### Dependencies

The package depends on a few third-party packages, which you may need to install first, if they are not present already. 
Most of these relate to specific functionality that you may not wish to use, thus are optional:

* Packages needed for computation

```
install.packages(c("data.table", "matrixStats"))
```

* Package needed for plotting results (optional, recommended)

```
install.packages("ggplot2")
```

* Packages needed for importing abundances from Salmon/Kallisto output (optional, recommended)

```
install.packages("devtools")

source("http://bioconductor.org/biocLite.R")

# Format converter from Salmon/Sailfish to Kallisto.
biocLite("COMBINE-lab/wasabi")

# Compressed format parser, for Kallisto.
biocLite("rhdf5")
```

* Package needed for interactive visualisation feature (optional)

```
install.packages("shiny")
```

If you have trouble installing these dependencies, your system could be missing source compilers for C and/or Fortran, and possibly other libraries, which you can see by scrolling back through the installation output to look for the errors. Please refer to the R manual or the respective package documentation for help.


### Installation

The recommended option is to obtain the latest platform-independent package release from the [releases section](https://github.com/bartongroup/Rats/releases) here on **Github**.
Download the package file and then install it using:

`install.packages("<path/to/downloaded/package>", repos = NULL, type="source")`

The latest release can also be installed directly from Github, using the `devtools` package (but this option has a tendency to install exhaustive minor dependencies or update existing packages, so it can take a while to complete):

`devtools::install_github("bartongroup/rats", ref="master")`

For testing purposes (bug resolutions, new features), you can install the ongoing developmental version from Github:

`devtools::install_github("bartongroup/rats", ref="development")`

From `v0.6.0` onwards, release versions of RATs continue to have a 3-part release number, whereas developmental versions now have a 4-part version number. Despite having version numbers, developmental versions are not archived and are **not suitable** for reproducible/critical/publishable analyses. They may also temporarily not work correctly or at all. For critical analyses you should use the latest release version.

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
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged. Also, be sure to **read the vignette(s)**, and browse/search
the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./vignettes/figs/rats_logo.png)


