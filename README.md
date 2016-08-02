# Relative Abundance of Transcripts (rats)

## Description

#### Who it is for

Anyone working in transcriptomics, analyzing gene expression and transcript abundancies.

#### What it does

It provides a method to detect changes in the relative abundance of the alternative transcripts (isoforms) of genes. 
This is called **Differential Transcript Usage (DTU)**.  

Detecting DTU is supplementary to the quantification of transcripts by tools like [Salmon](http://combine-lab.github.io/salmon/), 
[Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [Kallisto](http://pachterlab.github.io/kallisto/) and the detection 
of Differential Transcript Expression (DTE) by tools such as [Sleuth](http://pachterlab.github.io/sleuth/).

#### What it needs

As input it requires an R object similar to the output of [Sleuth](http://pachterlab.github.io/sleuth/). It also requires an index 
table matching transcript identifiers to respective gene identifiers.  

The package makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and 
[matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, as well as 
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) for visualisations. All these are
available from CRAN.


## How to use rats


A full **tutorial vignette** is included in the package, explaining the input, output, commands and options. 
If you install the package from one of the binaries or through Bioconductor, the vignette should then also be 
available locally by calling:

`browseVignettes("rats")`

If all else fails, you can view a static version of the tutorial [in the wiki](https://github.com/bartongroup/Rats/wiki).


### Dependencies

The package depends on a few additional CRAN packages, which you may need to install first, 
if they are not present already or if they are not installed automatically while installing `rats`.

`install.packages(c("data.table", "matrixStats", "ggplot2"))`


### Installation

Some **precompiled binary releases** are available from the [releases section](https://github.com/bartongroup/Rats/releases) on Github.
Download the appropriate binary for your system and then install it using:

`install.packages(<pathtopackage>, repos = NULL)` or
`install.packages(<pathtopackage>, repos = NULL, type="source")`

You can also install the **latest developmental version** directly from Github, using the [devtools](https://www.rstudio.com/products/rpackages/devtools/) package:

`devtools::install_github("bartongroup/rats")`

Eventually, we aim to make `rats` available through **Bioconductor** as well.


### Differential Transcript Usage

A typical command to call DTU looks like this:

`results <- call_DTU(my_sleuth_object, my_identifiers_table, "Condition-1", "Condition-2")`

Mandatory parameters:

* a sleuth object
* a dataframe matching unique transcript identifiers to gene identifiers
* the names of two conditions recorded in the sleuth object

For more details on the parameters, please refer to the tutorial vignette.


## Contact information

The rats R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourao and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged. Also, be sure to **read the vignette(s)**, and browse/search
the support forum before posting a new issue, in case your question is already answered there.

Enjoy!