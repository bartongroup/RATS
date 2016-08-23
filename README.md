# Relative Abundance of Transcripts (rats)

## Description

#### Who it is for

Anyone working in transcriptomics, analysing gene expression and transcript abundances.

#### What it does

It provides a method to detect changes in the relative abundance of the alternative transcripts (isoforms) of genes. 
This is called **Differential Transcript Usage (DTU)**.  

Detecting DTU is supplementary to the quantification of transcripts by tools like [Salmon](http://combine-lab.github.io/salmon/), 
[Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [Kallisto](http://pachterlab.github.io/kallisto/) and the detection 
of Differential Transcript Expression (DTE) by tools such as [Sleuth](http://pachterlab.github.io/sleuth/).

#### What it needs

This is an R package, therefore access to an R terminal is required for all the commands shown here and in the tutorial.

As input, `rats` requires transcript abundance estimates with or without bootstrapping. For convenience, these can also be extracted directly
from the output of [Sleuth](http://pachterlab.github.io/sleuth/). It also requires a look-up table matching transcript identifiers to 
respective gene identifiers.  

The package makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and 
[matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, as well as 
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) for visualisations. All these are
available from CRAN.


## How to use rats


A full **tutorial vignette** is included in the package, explaining the input, output, commands and options. 
If you install the package either from a compiled source release or through Bioconductor, the vignette should be 
available locally by calling:

`browseVignettes("rats")`

If this fails, you can view the vignette on-line directly [from its Github file, here](https://github.com/bartongroup/Rats/blob/master/vignettes/tutorial.md). 

We recommend studying the vignette before using `rats`.

### Dependencies

The package depends on a few additional CRAN packages, which you may need to install first, 
if they are not present already:

`install.packages(c("data.table", "matrixStats", "ggplot2"), dependencies=TRUE)`

If you have trouble installing these dependencies, your system could be missing source compilers for C and/or Fortran, and possibly other libraries, which you can see by scrolling back through the output to look for the errors. Please refer to the R manual for help.


### Installation

Platform-independent package releases are available from the [releases section](https://github.com/bartongroup/Rats/releases) on **Github**.
Download the latest release and then install it using:

`install.packages("<path/to/dowloaded/package>", repos = NULL, type="source")`

Eventually, we aim to make `rats` available through **Bioconductor** as well.


### Differential Transcript Usage

A typical command to call DTU looks like this:

`results <- call_DTU(annot = my_identifiers_table, slo = my_sleuth_object,  name_A = "Condition-1", name_B = "Condition-2")`

Mandatory parameters:

* a data frame matching unique transcript identifiers to gene identifiers
* a sleuth object
* the names of two conditions recorded in the sleuth object

`rats` also accepts data input in generic format that does not depend on Sleuth. Please consult the vignette for syntax and specifications.

The output is a list containing two tables that list the final results as well as the intermediate calculations and decisions.
For details on the parameters, input and output, please refer to the tutorial vignette.


## Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged. Also, be sure to **read the vignette(s)**, and browse/search
the support forum before posting a new issue, in case your question is already answered there.

Enjoy!
