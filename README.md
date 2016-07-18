# Relative Abundance of Transcripts (rats)

## Description

#### Who it is for

Anyone working in transcriptomics, analyzing gene expression and transcript abundancies.

#### What it does

It provides a method to detect changes in the relative abundance of the alternative transcripts (isomorphs) of genes. 
This is called **Differential Transcript Usage (DTU)**.  

Detecting DTU is supplementary to the quantification of transcripts by tools like [salmon](http://combine-lab.github.io/salmon/), 
[sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [kallisto](http://pachterlab.github.io/kallisto/) and the detection 
of Differential Transcript Expression (DTE) by tools such as [sleuth](http://pachterlab.github.io/sleuth/).

#### What it needs

As input it requires an R object similar to the output of [sleuth](http://pachterlab.github.io/sleuth/). It also requires an index t
able matching transcript identifiers to respective gene identifiers.  

The package makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and 
[matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, both available through CRAN.
Visalisations are done with [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html).

#### Who it is by

The rats R package was developed by members of [The Barton Group] (http://www.compbio.dundee.ac.uk) at [The University of Dundee] (http://www.dundee.ac.uk).


## How to use rats

### Obtain rats

The package depends on a few external packages, which you may need to install first, if they are not present already or if they are not i
nstalled automatically by installing rats.

```{r eval=FALSE}
install.packages("data.table")
install.packages("matrixStats")
install.packages("ggplot2")
```

#### from Github

Rats is hosted in GitHub [here](https://github.com/nickschurch/Rats.git), and you can dowload the zip file from that page or 
by [clicking here](https://github.com/nickschurch/Rats/archive/master.zip).
Then you can install the zipped file like so:

```{r eval=FALSE}
install.packages("/path/to/downloded.zip", repos=NULL)
```

#### from CRAN (not yet available)

#### from Bioconductor (not yet available)


### Loading rats

Once you have obtained rats through one of the above standard ways, you must then load it into R with the following command, 
before being able to use it:

```{r eval=FALSE}
library("rats")
```

### Calling Differential Transcript Usage with rats

A full tutorial vignette is bundled in the package, explaining the input, output, commands and options. It can be accessed with the
command below. Alternatively, it can be found by browsing the GitHub repository or the downloaded .zip archive in the `vignettes` directory.

```{r eval=FALSE}
browseVignettes("rats")
```

## Contact information

For further information or assistance with this repository please contact one of the authors (accurate as of June 2016):

* Nick Schurch: <n.schurch@dundee.ac.uk>
* Kimon Froussios: <k.froussios@dundee.ac.uk>
* Kira Mourao <k.mourao@dundee.ac.uk>
* Geoff Barton: <gjbarton@dundee.ac.uk>
