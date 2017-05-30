# Relative Abundance of Transcripts (rats)

***

## Description


#### Who it is for

Anyone working in transcriptomics, analysing gene expression and transcript abundances.


#### What it does

1. It provides a method to detect changes in the abundance ratios of transcript isoforms of a gene.
This is called **Differential Transcript Usage (DTU)**. 

2. **RATs is workflow-agnostic**. Quantification quality details are left to the quantification tools; RATs uses only the
transcript abundances. This makes it *suitable for use with alignment-free quantification tools* like [Kallisto](http://pachterlab.github.io/kallisto/)
or [Salmon](https://github.com/COMBINE-lab/salmon). It is also compatible with DTE output from [Sleuth](http://pachterlab.github.io/sleuth).

3. RATs is able to take advantage of the bootstrapped quantifications provided by the alignment-free tools. These bootstrapped
data are used by `RATs to assess how much the technical variability of the heuristic quantifications affects differential transcript usage
and thus provide a measure of confidence in the DTU calls. 


#### What it needs

1. This is an R source package, and will run on any platform with a reasonably up-to-date R environment.

2. As input, RATs requires transcript abundance estimates with or without bootstrapping. For convenience, these can also be extracted directly
from the output of [Sleuth](http://pachterlab.github.io/sleuth/). 

3. RATs also requires a look-up table matching transcript identifiers to respective gene identifiers. This can be obtained through various means,
one of them being extracting this info from a GTF file.

4. RATs makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and 
[matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, as well as 
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) and [shiny](https://cran.r-project.org/web/packages/shiny/shiny.pdf) for visualisations. All these are
available from CRAN.

***

## How to use RATs

### Documentation

The package comes with full documentation: The vignettes are available locally by calling `browseVignettes("rats")` after the package is installed.
We recommend studying the vignettes before using RATs.


### Dependencies

The package depends on a few third-party packages, which you may need to install first, 
if they are not present already:

* Packages needed for computation (mandatory)

```
install.packages(c("data.table", "matrixStats"), dependencies=TRUE)
```

* Package needed only for plotting (recommended)

```
install.packages("ggplot2", dependencies=TRUE)
```

* Packages needed only for importing abundances from Salmon/Kallisto output (optional) - 

```
install.packages("devtools", dependencies=TRUE)
source("http://bioconductor.org/biocLite.R")
# Wasabi converter from Salmon/Sailfish to Kallisto.
biocLite("COMBINE-lab/wasabi")

# Kallisto parser
biocLite("rhdf5")
```

* Package needed only for interactive visualisation feature (optional)

```
install.packages("shiny", dependencies=TRUE)
```

If you have trouble installing these dependencies, your system could be missing source compilers for C and/or Fortran, and possibly other libraries, 
which you can see by scrolling back through the installation output to look for the errors. Please refer to the R manual or respective package manuals for help.


### Installation

1. Platform-independent package releases are available from the [releases section](https://github.com/bartongroup/Rats/releases) on **Github**.
Download the latest release and then install it using:

`install.packages("<path/to/downloaded/package>", repos = NULL, type="source")`

The latest release can also be directly installed from github, using the `devtools` package:

`devtools::install_github("bartongroup/rats", ref="master")`

2. For testing purposes (bug resolutions, new features), you can install the ongoing developmental version from Github:
`devtools::install_github("bartongroup/rats", ref="development")`
Developmental versions are works in progress and will not be archived in snapshots. For reproducible/publishable analyses, 
**always use a release version**, NOT a developmental version.


### Differential Transcript Usage

A typical command to call DTU given a Sleuth object looks like this:

`results <- call_DTU(annot = my_identifiers_table, slo = my_sleuth_object,  name_A = "Condition-1", name_B = "Condition-2")`

RATs also accepts data input in **generic formats**. Please consult the vignettes for syntax details, format specifications and additional settings.

The output is a list containing (among other items) two tables that list the final results as well as the intermediate calculations and decisions.
Details on the output structure and visualisation options are provided in the vignettes.

***

## Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or ask for **assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged. Also, be sure to **read the vignette(s)**, and browse/search
the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./vignettes/fig/rats.png)


