---
title: "RATs: Input and Settings"
author: "Kimon Froussios"
date: "14 MAY 2017"
output: 
  html_document: 
    keep_md: yes
    theme: readable
    toc: yes
    toc_float: yes
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{RATs 2: Input & Settings}
  %\VignetteEngine{knitr::knitr}
---



***


```r
library(rats)
```

***

# Input formats

## Data

RATs can work with any of three input types:

1. A [Sleuth](http://pachterlab.github.io/sleuth/) object.
2. Bootstrapped abundance estimates.
3. Abundance estimates.

-> [Sleuth](http://pachterlab.github.io/sleuth/) objects contain their input data, so RATs can extract the bootstrapped abundance estimates directly from there, 
for your convenience if Sleuth is in your workflow already.

-> Bootstrapped abundance estimates can be input directly as `list`s of `data.table`s. Two lists are needed, one per condition.
Each datatable should contain the transcript identifiers in the first column, followed by columns containing the estimates from the bootstrap iterations:


```
##    target_id V1 V2 V3
## 1:     NIB.1  0  0  0
## 2:    1A1N-2 20 21 18
## 3:  1D1C:one  0  0  0
## 4:  1D1C:two 76 80 72
## 5:    1B1C.1  0  0  0
## 6:    1B1C.2 52 55 50
```

-> Abundance estimates, without bootstrapping information, can be input simply as two `data.table`s, one per condition. The first column should
contain the transcript identifiers, followed by columns listing the abundance per sample. The format of each table is identical to the one shown above.


## Annotation

Regardless of data format, RATs also needs an annotation `data.frame` or `data.table` that matches transcript identifiers to gene identifiers. This looks like this:


```
##   target_id parent_id
## 1     NIB.1       NIB
## 2    1A1N-2      1A1N
## 3  1D1C:one      1D1C
## 4  1D1C:two      1D1C
## 5    1B1C.1      1B1C
## 6    1B1C.2      1B1C
```

A function is provided to create this table, given a GTF file.
(**Note:** GFF3 is **not** supported, as the specification is too relaxed.)


```r
# Extract transcript ID to gene ID index from a GTF annotation.
myannot <- annot2ids("my_annotation_file.gtf")
```


***


# Calling DTU 

To bypass the complexity of running third-party tools in this tutorial, we will instead use **emulated data**. RATs comes
with data emulators, intended to be used for testing the code. However, they are also convenient to use for showcasing
how RATs works. If you happen to have real data available, you can use that instead if you wish.

By default, RATs reports on its progress and produces a summary report.
The progress messages and summary can be suppressed by adding the `verbose = FALSE` parameter to the call. 
To prevent cluttering this tutorial with verbose output, we will use this option in all the examples.

If you leave `verbose = TRUE` when trying out the examples below using the emulated data, you will get **warnings**
about the number of bootstraps. The warning is triggered because the emulated dataset used in the examples immitates only the structure of 
real data, not the actual volume of it, and as such contains too few bootstrap iterations.


### from abundance estimates, without bootstraps

This is the simplest usage case, provided only for completeness. We recommend using bootstrapped data whenever possible,
for reasons that will be discussed in the relevant section below.

First, let's emulate some data to work with:


```r
# Simulate some data.
simdat <- sim_count_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]       # Simulated abundances for one condition.
mycond_B <- simdat[[3]]       # Simulated abundances for other condition.
myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.
```

Now we can call DTU:


```r
# Find DTU between the simulated datasets.
mydtu <- call_DTU(annot= myannot, count_data_A= mycond_A, count_data_B= mycond_B, 
                  verbose= FALSE,
                  name_A= "healthy", name_B= "patients", varname= "My phenotype",
                  description="Comparison of two simulated counts datasets for the
                    tutorial. Simulated using built-in functionality of RATs.")
```

#### Mandatory arguents:

1. `annot` - An annotation data frame, as described in the *Input formats* section.
2. `count_data_A` and `count_data_B` - are each a `data.table` of transcript abundances as described in the Input section.

#### Optional arguments:

* `name_A`, `name_B` - A name for each conditon. (Default `NA`)
* `varname` - The name of the variable/condition. (Default `NA`)
* `description` - Free-text description. You can note experiment details, annotation source, anything you might find useful later. (Default `NA`)


### from bootstrapped abundance estimates

First, let's emulate some data, as we did before.


```r
# Simulate some data. (Notice it is a different function than before.)
simdat <- sim_boot_data()

# For convenience let's assign the contents of the list to separate variables.
mycond_A <- simdat[[2]]       # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]       # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]        # Transcript and gene IDs for the above data.
```

Now we can call DTU:


```r
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot= myannot, boot_data_A= mycond_A, boot_data_B= mycond_B, 
                  verbose= FALSE, 
                  name_A= "wildtype", name_B= "some mutant", varname = "My phenotype",
                  description="Comparison of two simulated datasets of bootstrapped 
                    counts for the tutorial. Simulated using built-in functionality 
                    of RATs.")
```

#### Mandatory arguents:

1. `annot` - An annotation data frame, as described in the *Input formats* section.
2. `boot_data_A` and `boot_data_B` - are each a list of `data.table` objects, as described in the Input section.

#### Optional arguments:

* `name_A`, `name_B` - A name for each conditon. (Default `NA`)
* `varname` - The name of the variable/condition. (Default `NA`)
* `description` - Free-text description. You can note experiment details, annotation source, anything you might find useful later. (Default `NA`)


### with a sleuth object

First, let's emulate a Sleuth object, using the bundled tools. The real Sleuth objects are very large and very complex nested lists. 
The emulated one contains only the minimum essential elements relevant to calling DTU with RATs. 


```r
# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# controls and patients are arbitrary names to use as conditions.

# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo       # Simulated minimal sleuth object.
myannot <- simdat$annot   # Transcript and gene Identifiers for the above data.
```

Now call DTU.


```r
# Find DTU between conditions "controls" and "patients" in the simulated data.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object 
                    for the purposes of the tutorial. Simulated using built-in functionality 
                    of RATs.")
```

#### Mandatory arguents:

1. `annot` - an annotation data frame, as described in the *Input formats* section.
2. `slo` - a Sleuth object.
3. `name_A` and `name_B` - the two groups of samples to compare. These are restricted to values available in `myslo$sample_to_covariates`. The values must both exist in the column specified by `varname`.
4. `varname` - the variable for which samples are compared. must be an existing column header in `myslo$sample_to_covariates`.

#### Optional arguments:

* `description` - Free-text description. You can note experiment details, annotation source, anything you might find useful later. (Default `NA`)

Please note that, unlike the other two usage cases, `name_A`, `name_B` and `varname` are **not optional**, as they specify how data is extracted from the Sleuth object.


### with a sleuth object on limited RAM

Sleuth objects can take up several GBs of memory. RATs extracts the abundances from the sleuth object but then restructures them so as to access them efficiently.
As a consequence, the bootstrapped abundances will be held in memory twice, once for internal use by RATs and once for the sleuth object. This is a waste fo resources,
and can be a problem on limited systems like personal laptops.

To allow you to free up RAM by unloading the sleuth object before calling DTU, RATs now exports some of its previously internal functions, allowing the separation of
the DTU call step from the data restructuring step.

You will need a sleuth object, the annotation table, and a vector of sample numbers.
We will be using the simulated data created in the sleuth section above:


```r
# 1. Restructure sample_to_covariates as covariate to samples.
samples_by_covariate <- group_samples(myslo$sample_to_covariates)

# There are two covariates in the simulated data:
print( names(samples_by_covariate) )
```

```
## [1] "condition" "batch"
```

```r
# Covariate "condition" in our simulated data has two values:
print( names( samples_by_covariate$condition ) )
```

```
## [1] "controls" "patients"
```

```r
# 2. Extract the bootstrapped abundance tables for each of the two conditions:
condA_boots <- denest_sleuth_boots(myslo, myannot$target_id, samples_by_covariate$condition[["controls"]])
condB_boots <- denest_sleuth_boots(myslo, myannot$target_id, samples_by_covariate$condition[["patients"]])

# 3. Remove the sleuth object from memory.
# Make sure you have saved a copy of it to file before doing this!
rm(myslo)

# 4. Run RATs with the generic bootstrapped format:
mydtu <- call_DTU(annot= myannot, boot_data_A= condA_boots, boot_data_B= condB_boots, 
                  verbose= FALSE, 
                  name_A= "controls", name_B= "patients", varname = "condition",
                  description="Comparison of two sets of bootstrapped counts for the 
                              tutorial, extracted from a simulated sleuth object. 
                              Simulated using built-in functionality of RATs.")
```

#### Mandatory arguents (group_samples):

1. `covariates` - any data frame in the format of a sleuth object's `sample_to_covariates` table. Each row is a sample, each column is a covariate,
and the values in the cells represent the category of the covariate for the sample.

#### Mandatory arguents (denest_sleuth_boots):

1. `slo` - A Sleuth object.
2. `tx` - The transcript ID column of the annotation table you will use for calling DTU, to ensure/enforce consistent number and order of transcripts across all tables. 
3. `samples` - A vector of sample numbers (their index in the sleuth object). Can be obtained by subsetting the output of `group_samples()``.


NOTE: In this example `rm(myslo)` won't free up the memory, because there are other active references to it (it is part of the `simdat` object).


***


# Additional Parameters and Settings


In summary, `rats` works as follows:

![RATs design](./fig/RATs_decisions.png)


## Thresholds

The following three main thresholds are used in RATs:


```r
# Calling DTU with custom thresholds.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  p_thresh = 0.01, abund_thresh = 10, dprop_thresh = 0.25)
```

1. `p_thresh` - Statistical significance level. P-values below this will be considered significant. Lower threshold values are stricter. (Default 0.05)
2. `abund_thresh` - Noise threshold. Transcripts with abundance below that value in both conditions are ignored. Higher threshold values are stricter. (Default 5, assumes use of estimated counts similar to those reported by [Salmon](https://github.com/COMBINE-lab/salmon))
3. `dprop_thresh` - Effect size threshold. Transcripts whose proportion changes between conditions by less than the threshold are considered non-DTU, regardless of their statistical significance. Higher threshold values are stricter. (Default 0.20)

The default values for these thresholds have been chosen such that they achieve a *median* FPR <5% for a high quality dataset from *Arabidopsis thaliana*, even with only 3 replicates per condition. Your mileage may vary and you should give some consideration to selecting appropriate values.
Depending on the settings, *additional thresholds* are available and will be discussed in their respective sections below.


## Bootstrapping 

RATs offers two types of bootstrapping:

1. Bootstrapping of `p-value` and `Dprop` against the technical fluctuations of the quantifications. This requires bootstrapped quantifications as input.
2. Bootstrapping of `p-value` and `Dprop` against the abundance fluctuations among samples.

Enabling these two procedures assesses the robustness of the DTU calls. Calls that are easily overturned depending on which subset of the data is used
may be spurious but may also be valid. Deeper sequencing and/or higher replication may resolve these ambiguities.

In both cases, what is measured is the fraction of iterations in which the significance and effect size meet their respective thresholds. 
These fractions can be used as indicators of confidence in the calls, especially when a large number of iterations is performed.

If bootstrapping is switched off, the respective output fields will not be included in the output.


### Quantification bootstraping

Three parameters control bootstrapping of DTU calls on the abundance estimates:

1. `qboot` - Whether to bootstrap the quantifications or not. (Default TRUE)
2. `qbootnum` - How many bootstrap iterations to do. Preferably at least 100. If 0, RATs will try to infer a value from the data. (Default 0)
3. `qrep_thresh` - Reproducibility threshold. What fraction of the iterations has to agree on a result to consider it confident. Higher threshold values are stricter. (Default 0.95)

In this process, one quantification iteration will be randomly selected from each sample and DTU will be called on it. This will be repeated `qbootnum` times. Because the
number of replicates remains the same, the statistical power is not compromised. Therefore, the reproducibility will be **used as a criterion** in calling DTU, along with
statistical significance and effect size. The process is **stochastic**; the quantifications are randomly sampled, so runs on the same data may yield slight differences in DTU. 
Using higher `qbootnum` improves reproducibility between runs on the same dataset.

Low reproducibility indicates that the quantification tool found it hard to distinguish these transcripts. This can be caused by high similarity of the isoforms, genes with a large
number of isoforms, and/or poor read coverage in the regions differentiating the isoforms from one another. In these cases, skepticism is required about the quantifications and any 
potential differential expression regarding these transcripts.

Warnings will be generated if `qbootnum` is too low or too high, but in most cases RATs will continue with the analysis. The warnings will not be shown if `verbose = FALSE`.


```r
# Bootstrap (default). Do 100 iterations.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", qboot = TRUE, qbootnum = 100, qrep_thresh= 0.95)

# Skip bootstraps.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", qboot = FALSE)
```


### Replicate bootstraping

Three parameters control bootstrapping of DTU calls agaisnt the samples:

1. `rboot` - Whether to bootstrap the replicates or not. (Default FALSE)
2. `rrep_thresh` - Reproducibility threshold. What fraction of the iterations has to agree on a result to consider it confident. Higher threshold values are stricter. (Default 0.85)
3. `rrep_as_crit` - Whether the replicate-wise reproducibility should be used as a criterion in calling DTU. (Default FALSE)

Unlike bootstrapping the quantifications, bootstrapping the replicates is currently **not stochastic**. RATs will do ALL the 1 vs 1 combinations of samples, one from each condition.
This behaviour may be changed if in the future there is demand for high replication that reaches the capacity limits of R matrices.

Because only one replicate is used each time, this process has reduced statistical power and only the strongest results will persist. Due to this, the sample-wise reproducibility
is **not used as a criterion** in calling DTU, unless explicitly required by setting the `rrep_as_crit` flag. For low replication levels in particular, the reproducibility may
be easily underestimated or overestimated.


```r
# Bootstrap (default). 
# NOTE: The number of iterations for 3 samples per condition is 3*3=9. 
# So the minimum possible error rate is 1/9=0.111. The corresponding threshold is 1-0.111=0.888.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", rboot = TRUE, qrep_thresh= 0.85)

# Skip bootstraps.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", rboot = FALSE)
```


## Multi-threading

RATs completion time depends on the number of expressed annotated transcripts. Single-threaded, RATs can take up to a few minutes per iteration for large annotations.
Fortunately, the task is highly parallelisable:


```r
# Using 8 threads/cores for parallel computing.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", threads = 8)
```

1. `threads` - The number of threads to use. (Default 1)

Due to core R implementation limitations, multi-threading works only in POSIX-compliant systems. In other systems, the option will be ignored and R will run serially.
Most Linux distros, Unix versions, and recent OSX versions are supposedly compatible, and there are environments that can be installed on Windows that might circumvent the issue.


## Test selection

RATs runs both gene-level calls and transcript-level calls. They are complementary and we recommend using them
together, but the option to skip either is provided for special use cases. The fields of the skipped test will be filled with `NA`.


```r
# Transcripts only.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", testmode="transc")
# Genes only.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", testmode="genes")
```

1. `testmode` - Which test(s) to run {"transc", "genes", "both""}. (Default "both")

## Correction for multiple testing

Testing multiple null hypotheses increases the chance of one being falsely rejected. To keep the overall false rate at the 
desired level, the raw p-values must be adjusted. The default adjustment method is `BH` (Benjamini-Hochberg). A full list of 
options is listed in R's `p.adjust.methods`.


```r
# Bonferroni correction.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", 
                  name_B = "patients", correction = "bonferroni")
```

1. `correction` - Type of multiple testing correction. (Default "BH")

## Input field names

RATs needs to pull information from different fields of the data and annotation and it does so based on the names of the columns. 
You can override the default names of these fields.



```r
# Lets emulate some input with custom field names. 
sim <- sim_sleuth_data(varname="mouse", cnames=c("Splinter", "Mickey"), 
                       COUNTS_COL="the-counts", TARGET_COL="transcript", 
                       PARENT_COL="gene", BS_TARGET_COL = "trscr")
myslo <- sim$slo
myannot <- sim$annot

# Sleuth covariates table.
print( sim$slo$sample_to_covariates )
```

```
##      mouse batch
## 1 Splinter    ba
## 2   Mickey    ba
## 3 Splinter    bb
## 4   Mickey    bb
```

```r
# Sleuth bootstrapped quantifications.
print( head(sim$slo$kal[[1]]$bootstrap[[1]]) )
```

```
##      trscr the-counts
## 1      LC1          3
## 2     NIA1        333
## 3     NIA2        666
## 4   1A1N-1         10
## 5   1A1N-2         20
## 6 1D1C:one          0
```

```r
# Annotation.
print( head(sim$annot) )
```

```
##   transcript gene
## 1      NIB.1  NIB
## 2     1A1N-2 1A1N
## 3   1D1C:one 1D1C
## 4   1D1C:two 1D1C
## 5     1B1C.1 1B1C
## 6     1B1C.2 1B1C
```


### Sleuth field names

Although we expect Sleuth field names to be constant with the exception of covariate names, RATs allows overriding the
expected names:


```r
# Call DTU on data with custom field names.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey", 
                  varname="mouse", COUNTS_COL="the-counts", BS_TARGET_COL="trscr", verbose = FALSE)

#! In our example data here, the annotation fields have been modified as well, so this command 
#! will NOT run as it is shown. You will also need the parameters shown in teh section below.
```

1. `varname` - The field name in `myslo$sample_to_covariates` where the desired condition names are listed. (Default "condition")
2. `COUNTS_COL` - The name of the field holding the estimated counts in the Sleuth object's bootstrap tables. (Default "est-counts")
3. `BS_TARGET_COL` - The name of the field holding the transcript identifiers in the Sleuth object's bootstrap tables. (Default "target_id")

The `COUNTS_COL` and `BS_TARGET_COL` parameters are also available for `denest_sleuth_boots()`:


### Annotation field names

Although it is easy to rename the columns of a table to comply with the expected names, this may sometimes be undesireable, so RATs
allows you to change the expected names instead.


```r
# Call DTU using annotation with custom field names.
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "Splinter", name_B = "Mickey", 
                  varname="mouse", TARGET_COL="transcript", PARENT_COL="gene", verbose = FALSE)

#! In our example data here, the Sleuth fields have been modified as well, so this command 
#! will NOT run as it is shown. You will also need the parameters shown in the section above.
```

1. `TARGET_COL` - The name of the field holding the transcript identifiers in the annotation data frame. (Default "target_id")
2. `PARENT_COL` - The name of the field holding the respective gene identifiers in the annotation data frame. (Default "parent_id")


***


# Annotation discrepancies

Different annotation versions often preserve transcript IDs, despite altering the details of the transcript models,
such as UTRs. It is important to use the same annotation throughout the workflow, especially for the quantifications.
Otherwise, the data is not comparable.

All internal operations and the output of RATs will be based on the annotation provided:

* Any transcripts present in the data but missing from the annotation will be ignored completely and 
will not show up in the output, as there is no reliable way to match them to gene IDs.
* Any transcript/gene present in the annotation but missing from the data will be included in the output as zero expression.
* If the samples appear to use different annotations from one another, RATs will abort.


***


# Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged to help us understand the problem. Also, be sure 
to **read the vignette(s)**, and browse/search the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./fig/rats.png)
