# Relative Abundance of Transcripts (rats)

## Description

#### Who it is for

Anyone working in transcriptomics, analyzing gene expression and transcript abundancies.

#### What it does

It provides a method to detect changes in the relative abundance of the alternative transcripts (isoforms) of genes. 
This is called **Differential Transcript Usage (DTU)**.  

Detecting DTU is supplementary to the quantification of transcripts by tools like [salmon](http://combine-lab.github.io/salmon/), 
[sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [kallisto](http://pachterlab.github.io/kallisto/) and the detection 
of Differential Transcript Expression (DTE) by tools such as [sleuth](http://pachterlab.github.io/sleuth/).

#### What it needs

As input it requires an R object similar to the output of [sleuth](http://pachterlab.github.io/sleuth/). It also requires an index 
table matching transcript identifiers to respective gene identifiers.  

The package makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and 
[matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, as well as 
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) for visualisations. All these are
avoilable from CRAN.


## How to use rats

#### Dependencies

The package depends on a few additional CRAN packages, which you may need to install first, 
if they are not present already or if they are not installed automatically while installing `rats`.

`install.packages(c("data.table", "matrixStats", "ggplot2"))`

#### Obtaining rats

**from Github:** `rats` is hosted on GitHub [here](https://github.com/bartongroup/Rats), and you can 
dowload the zip file from that page or by [clicking here](https://github.com/bartongroup/Rats/archive/master.zip). 
Then you can install from the zipped file like so:

`install.packages("/full/path/to/downloded/rats.zip", repos=NULL)`


#### Calling Differential Transcript Usage with rats

A full tutorial vignette is included in the package, explaining the input, output, commands and options. It can be accessed with the
command below. Alternatively, it can be found by browsing to the `vignettes` directoy in the GitHub repository or the downloaded 
.zip archive.

`browseVignettes("rats")`

## Contact information

The rats R package was developed by members of [The Barton Group] (http://www.compbio.dundee.ac.uk) at [The University of Dundee] (http://www.dundee.ac.uk).
It was conceived by Dr. Nick Schurch, and implemented by Dr. Kimon Froussios and Dr. Kira Mourao.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged. Also, be sure to **read the vignette(s)**, and browse/search
the support forum before posting a new issue, in case your question is already answered there.
