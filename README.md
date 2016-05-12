# Relative Abundance of Transcripts (rats)

This is the Github repository for the rats R package developed by members of [The Barton Group] (http://www.compbio.dundee.ac.uk) at [The University of Dundee] (http://www.dundee.ac.uk) to compute the relative usage of the transcripts associated with a gene from the output of transcript counting tools like [salmon](http://combine-lab.github.io/salmon/), [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) and [kallisto](http://pachterlab.github.io/kallisto/). The idea here is to identify genes that show Differential Transcript Usage (DTU, as opposed to Differential Transcript Expression; DTE) within the transcript abundance estimate data.

The aim is to make this a simple package available in [bioconductor](http://bioconductor.org/).

# How to use rats

## Installation

Checkout the rats package from github and install in R. The simplest way to do this is to use [rstudio](https://www.rstudio.com/) to chekout, build and install the package directly.

## loading and using rats

Once you ahve installed the library, load rats with:

```r
library(rats)
```

With the library loaded, you can now use it to identify the Differential Transcript Usage. 

First we will need a sleuth object containing the transcript abundance estimate data. See the [introduction to sleuth](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html) pages for details on how to load the transcript abundance estimate data from [kallisto](https://pachterlab.github.io/kallisto/) runs into a sleuth object, and see the [wasabi](https://github.com/COMBINE-lab/wasabi) tool for how to load the transcript abundance estimate data from   [Sailfish](https://github.com/kingsfordgroup/sailfish) or [Salmon](https://github.com/COMBINE-lab/salmon) into a sleuth object.

For our example we will call our sleuth object, `so`. 

Next we will need a data frame that maps the transcript IDs to their parent gene IDs. This needs to have at least two variables - `target_id` and `parent_id` - but can also contain other annotation variable as well.

With our annotation data and out sleuth object defined, calling DTU on this data is then as easy as:

```r
DTU = calculate_DTU(so, t2g, "x", "y")
```

where `x` is the name of the reference condition to use (typically the 'wild-type')  and `y` is the name of the condition to compare with the reference (for example a mutant or drug treatment).

The output is a List object with three elements: 

1. `Comparison`: This contains the details of the comparison that was run.
2. `Genes`: a data frame summarising the identified DTU at the gene level
3. `Transcripts`: a data table containing the transcript level evidence that was used to compute the DTU.

## Examining the results

With the DTU object calculated, the DTU results for the transcripts associated with a given gene can be investigated graphically with:

```r
plotGeneDTU(DTU, "geneid", nreps=7, ptype="proportion")
```

where `geneid` is the parent_id that identifies the gene and nreps is the number of biological replicates in the study.

## Contact information

For further information or assistance with this repository please contact one of:

* Nick Schurch: <nschurch@dundee.ac.uk>
* Kimon Froussios: <k.froussios@dundee.ac.uk>
* Kira Mourao <k.mourao@dundee.ac.uk>
* Geoff Barton: <gjbarton@dundee.ac.uk>
