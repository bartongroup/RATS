# Relative Abundance of Transcripts (rats)

## Description

#### Who it is for

Anyone working in transcriptomics, analyzing gene expression and transcript abundancies.

#### What it does

It provides a method to detect changes in the relative abundance of the alternative transcripts (isomorphs) of genes. This is called **Differential Transcript Usage (DTU)**.  

Detecting DTU is supplementary to the quantification of transcripts by tools like [salmon](http://combine-lab.github.io/salmon/), [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [kallisto](http://pachterlab.github.io/kallisto/) and the detection of Differential Transcript Expression (DTE) by tools such as [sleuth](http://pachterlab.github.io/sleuth/).

#### What it needs

As input it requires an R object similar to the output of [sleuth](http://pachterlab.github.io/sleuth/). It also requires an index table matching transcript identifiers to respective gene identifiers.  

The package makes use of the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and [matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html) packages, both available through CRAN.

#### Who it is by

The rats R package was developed by members of [The Barton Group] (http://www.compbio.dundee.ac.uk) at [The University of Dundee] (http://www.dundee.ac.uk).


## How to use rats

### Obtain rats

#### from Github

If you use git or SVN, you can [checkout/clone](https://github.com/nickschurch/Rats.git) the repository. Otherwise you can [download as .zip](https://github.com/nickschurch/Rats/archive/master.zip) and extract it to a destination of your choice.

The package depends on two external packages, [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and [matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html), so you will also need to obtain these first:

```{r eval=FALSE}
install.packages("data.table")
install.packages("matrixStats")
```

#### from Bioconductor (not yet available)

Eventually we hope to make it available through [Bioconductor](https://bioconductor.org/). Instructions on how to install packages through Bioconductor can be found on the [Bioconductor installation guide](https://www.bioconductor.org/install/). Installing through Bioconductor should take care of dependencies automatically.

### Loading rats

Once you have obtained rats through one of the above standard ways, you must then load it into R with the following command, before being able to use it:

```{r eval=FALSE}
library("rats")
```

or 

```{r eval=FALSE}
library("path/to/rats")
```

### Detecting Differential Transcript Usage with rats

First we will need a sleuth object containing the transcript abundance estimate data. The input format recognized is the output of
[sleuth](http://pachterlab.github.io/sleuth/). See the [introduction to sleuth](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html) pages for details on how to load the transcript abundance estimate data from [kallisto](https://pachterlab.github.io/kallisto/) into a sleuth object, and see the [wasabi](https://github.com/COMBINE-lab/wasabi) tool for how to load the transcript abundance estimate data from   [Sailfish](https://github.com/kingsfordgroup/sailfish) or [Salmon](https://github.com/COMBINE-lab/salmon) into a sleuth object.

We will also need a data frame that maps the transcript IDs to their parent gene IDs. This needs to have at least two variables - `target_id` and `parent_id` - but can also contain other annotation variables as well (which will be ignored by the current implementation). It is also possible to use different column names and override the default parameter values when calling `calculate_DTU()`. Each row must represent one transcript and provide a `target_id` and a `parent_id` for it.

For our example we will call our sleuth object `so` and the identifiers table `t2g`. 
With our annotation data and our sleuth object defined, calling DTU on this data is then as easy as:

```{r eval=FALSE}
DTU <- calculate_DTU(so, t2g, "x", "y")
```

where `x` and `y` are the names of the conditions to compare, as they appear in the `sample_to_covariates` table within the `so` object. For additional parameters, please refer to the documentation of this function.

Depending on the number of transcripts in the annotation, this calculation can become quite slow. We are looking into ways to improve performance there. In the meantime, if you are in doubt, try with a reduced annotation to verify that the program works properly in your environment.

The output is a List object with three elements: 

1. `Parameters`: A list that contains the details of the comparison that was run.
2. `Genes`: A data table summarising the identified DTU at the gene level.
3. `Transcripts`: A data table containing the transcript level evidence that was used to compute the DTU.

### Examining the results of rats

With the DTU object calculated, the DTU results for the transcripts associated with a given gene can be investigated graphically with:

```{r eval=FALSE}
plotGeneDTU(DTU, "geneid", nreps=7, ptype="proportion")
```

where `geneid` is the parent_id that identifies the gene and `nreps` is the number of biological replicates in the study.

More visualization options may be made available in the future. Suggestions are welcome!

### Example data for rats

Example input and output structures are provided in the `data` subdirectory. They should be automatically loaded when requested by name (without the .rda file extension). Otherwise try loading them explicitly:

```{r eval=FALSE}
load("path/to/file.rda")
```

### Documentation for rats

Documentation for the package's functions can be obtained in the standard R way, for example:

```{r eval=FALSE}
?calculate_DTU
```


## Contact information

For further information or assistance with this repository please contact one of the authors (accurate as of June 2016):

* Nick Schurch: <n.schurch@dundee.ac.uk>
* Kimon Froussios: <k.froussios@dundee.ac.uk>
* Kira Mourao <k.mourao@dundee.ac.uk>
* Geoff Barton: <gjbarton@dundee.ac.uk>
