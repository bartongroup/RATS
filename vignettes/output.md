---
title: 'RATs: Raw Output'
author: "Kimon Froussios"
date: "08 JUL 2019"
output:
  html_document:
    keep_md: yes
    theme: readable
    toc: yes
    toc_float: yes
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{RATs 2: Raw Output} 
  %\VignetteEngine{knitr::rmarkdown} 
---



***


```r
library(rats)
library(data.table)
```


Let's set up an example, using RAT's data emulator (used for code testing, not suitable for data simulations).


```r
# Simulate some data.
simdat <- sim_boot_data(clean=TRUE) 
# For convenience let's assign the contents of the list to separate variables
mycond_A <- simdat[[2]]   # Simulated bootstrapped data for one condition.
mycond_B <- simdat[[3]]   # Simulated bootstrapped data for other condition.
myannot <- simdat[[1]]    # Transcript and gene IDs for the above data.

# Call DTU
mydtu <- call_DTU(annot= myannot, verbose= FALSE,
                  boot_data_A= mycond_A, boot_data_B= mycond_B,
                  dprop_thresh=0.1, qboot=TRUE, rboot=FALSE)
```


***


# Output structure

The output of RATs is a list containing 4 elements:


```r
print( names(mydtu) )
```

```
## [1] "Parameters"  "Genes"       "Transcripts" "Abundances"
```


## Parameters

`Parameters` is a list that contains information about the data and the settings for a particular run.


```r
# Parameter list's elements.
print( names(mydtu$Parameters) )
```

```
##  [1] "description"         "time"                "rats_version"       
##  [4] "R_version"           "var_name"            "cond_A"             
##  [7] "cond_B"              "data_type"           "num_replic_A"       
## [10] "num_replic_B"        "num_genes"           "num_transc"         
## [13] "tests"               "use_sums"            "correction"         
## [16] "p_thresh"            "abund_thresh"        "dprop_thresh"       
## [19] "abund_scaling"       "quant_boot"          "quant_reprod_thresh"
## [22] "quant_bootnum"       "rep_boot"            "rep_reprod_thresh"  
## [25] "rep_bootnum"         "seed"                "reckless"           
## [28] "lean"
```

* `description` - (str) Free-text description of the run. It is useful to record data sources, annotation source and version, experimental parameters, etc...
* `time` - (str) Date and time for the run. This does not represent the exact time the run started.
* `rats_version` - (str) The version of RATs.
* `R_version` - (str) The version of R (including OS architecture).
* `var_name` - (str) The value passed to the `varname` parameter.
* `cond_A` & `cond_B` - (str) The values passed to the `name_A` and `name_B` parameters.
* `data_type` - (str) The format of the input data.
* `num_replic_A` & `num_replic_B` - (int) The number of samples in each condition.
* `num_genes` - (int) The number of genes in the provided annotation.
* `num_transc` - (int) The number of transcripts in the provided annotation.
* `tests` - (str) The value passed to the `testmode` parameter.
* `use_sums` - (bool) Sum the replicates instead of averaging them (pre-0.7.0 RATs always summed).
* `p_thresh` - (num) The value passed to the `p_thresh` parameter.
* `abund_thresh` - (num) The value passed to the `abund_thresh` parameter.
* `dprop_thresh` - (num) The value passed to the `dprop_thresh` parameter.
* `correction` - (str) Multiple testing correction method.
* `abund_scaling` - (num) The value(s) passed to the `scaling` parameter. It will be either a single numeric value or a named vector of values.
* `quant_boot` - (bool) The value passed to the `qboots` parameter.
* `quant_reprod_thresh` - (num) The value passed to the `qrep_thresh` parameter, if `qboot==TRUE`.
* `quant_bootnum` - (int) The value passed to the `qbootnum` parameter, if `qboot==TRUE`.
* `rep_boot` - (bool) The value passed to the `rboot` parameter.
* `rep_reprod_thresh` - (num) The value passed to the `rrep_thresh` parameter, if `rboot==TRUE`.
* `rep_bootnum` - (int) The number of replicate bootstrapping iterations (currently  M*N, where M and N are the numbers of samples in the two conditions), if `rboot==TRUE`.
* `seed` - (int) Custom seed for the random generator. `NA` if none was specified.
* `reckless` - (bool) The value passed to the `reckless` parameter.
* `lean` - (bool) The value passed to the `lean` parameter.

**Note:** If bootstraps are disabled, the bootstrap-related fields may be `NA`, regardless of supplied values, to reflect the fact that they were not used.


## Genes

`Genes` is a [data.table](https://cran.r-project.org/web/packages/data.table/) listing 
results at the gene level. For your convenience, the respective aggregated transcript-level DTU calls are 
also included here (defined as at least one isoform being called DTU individually).


```r
# Genes table's fields.
print( names(mydtu$Genes) )
```

```
##  [1] "parent_id"      "elig"           "sig"            "elig_fx"       
##  [5] "quant_reprod"   "DTU"            "transc_DTU"     "known_transc"  
##  [9] "detect_transc"  "elig_transc"    "maxDprop"       "pval"          
## [13] "pval_corr"      "quant_na_freq"  "quant_dtu_freq"
```

The first few columns show the result of each decision step in a boolean manner, allowing easy filtering of the table. 
The remaining columns list the values based on which these decisions were made.
Pseudo-code formulas are shown here to help understand how the different fields interact when making decisions.

* `parent_id` - (str) Gene identifier.
* `elig` - (bool) Eligible for testing. Whether the gene met the pre-filtering criteria. `= (elig_transc >= 2)`.
* `sig` - (bool) Statistically significant. `= (pval_corr < Parameters$p_thresh)`.
* `elig_fx` - (bool) Eligible effect size. Whether at least one of the isoforms meets the effect size criterion. `= any(Transcript[parent_id, elig_fx])`.
* `quant_reprod` - (bool) Result reproducible across quantification iterations.
For positive DTU, `= (quant_dtu_freq >= Parameters$quant_reprod_thresh)`, for non-DTU `= (quant_dtu_freq <= 1 - Parameters$quant_reprod_thresh)`.
* `rep_reprod` - (bool) Result reproducible across replicates.
For positive DTU, `= (rep_dtu_freq >= Parameters$rep_reprod_thresh)`, for non-DTU `= (rep_dtu_freq <= 1 - Parameters$rep_reprod_thresh)`.
* `DTU` - (bool) Gene-level Differential Transcript Usage. `= (sig & elig_fx & quant_reprod & rep_reprod)`.
* `transc_DTU` - (bool) Aggregated from `Transcripts$DTU`: At least one isoform was individually reported as DTU. `= (Transcript[ , any(DTU), by=parent_id])`.
* `known_transc` - (int) Number of known trascripts for this gene according to the given annotation.
* `detect_trancs` - (int) Number of detected (expressed) transcripts in the given dataset.
* `elig_transc` - (int) Number of eligible transcripts, aggregated from `Transcripts$elig`.
* `pval` - (num) G test of independence p-value for the isoform ratios of the gene.
* `pval_corr` - (num) Multiple testing corrected p-value. `pval_corr= p.adjust(pval, Parameters$correction)`.
* `quant_p_median` - (num) Median of corrected p-values across quantification bootstraps (not available in `lean` mode).
* `quant_p_min` - (num) Minimum observed (corrected) p-value across quantification bootstraps (not available in `lean` mode).
* `quant_p_max` - (num) Maximum observed (corrected) p-value across quantification bootstraps (not available in `lean` mode).
* `quant_na_freq` - (num) Fraction of quantification iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
* `quant_dtu_freq` - (num) Fraction of replicate iterations that support a positive `DTU` classification.
* `rep_p_median` - (num) Median p-value across replication bootstraps (not available in `lean` mode).
* `rep_p_min` - (num) Minimum observed p-value across replication bootstraps (not available in `lean` mode).
* `rep_p_max` - (num) Maximum observed p-value across replication bootstraps (not available in `lean` mode).
* `rep_na_freq` - (num) Fraction of replication iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
* `rep_dtu_freq` - (num) Fraction of replication iterations that support a positive `DTU` classification.

**Note:** The fields reporting on the bootstraps will not be shown when bootstrapping is disabled.


## Transcripts

`Transcripts` is a [data.table](https://cran.r-project.org/web/packages/data.table/) listing 
results at the transcript level. For your convenience, the respective gene-level DTU calls are also included here.


```r
# Transcripts table's fields.
print( names(mydtu$Transcripts) )
```

```
##  [1] "target_id"      "parent_id"      "elig_xp"        "elig"          
##  [5] "sig"            "elig_fx"        "quant_reprod"   "DTU"           
##  [9] "gene_DTU"       "abundA"         "abundB"         "stdevA"        
## [13] "stdevB"         "log2FC"         "totalA"         "totalB"        
## [17] "propA"          "propB"          "Dprop"          "pval"          
## [21] "pval_corr"      "quant_na_freq"  "quant_dtu_freq"
```

The first few columns show the result of each decision step in a boolean manner, allowing easy filtering of the table. 
The remaining columns list the values based on which these decisions were made.
Pseudo-code formulas are shown here to help understand how the different fields interact when making decisions.

* `target_id` - (str) Transcript identifier.
* `parent_id` - (str) Gene identifier.
* `elig_xp` - (bool) Eligible expression level. `= (meanA >= Parameters$abund_thresh | meanB >= Parameters$abund_thresh)`.
* `elig` - (bool) Eligible for testing (meets the noise threshold and at least one other isoform is expressed). `= (elig_xp & totalA != 0 & totalB != 0 & (abundA != totalA | abundB != totalB))`.
* `sig` - (bool) Statistically significant. `= (pval_corr < Parameters$p_thresh)`.
* `elig_fx` - (bool) Eligible effect size. Proxy for biological significance. `= (Dprop > Parameters$dprop_thresh)`.
*. `quant_reprod` - (bool) Quantification reproducible.
For positive DTU, `= (quant_dtu_freq >= Parameters$quant_reprod_thresh)`, for non-DTU `= (quant_dtu_freq <= 1 - Parameters$quant_reprod_thresh)`.
* `rep_reprod` - (bool) Replication reproducible.
For positive DTU, `= (rep_dtu_freq >= Parameters$rep_reprod_thresh)`, for non-DTU `= (rep_dtu_freq <= 1 - Parameters$rep_reprod_thresh)`.
* `DTU` - (bool) The Transcript is Differentially Used. `= (sig & elig_fx & quant_reprod & rep_reprod)`.
* `gene_DTU` - (bool) Expanded from `Genes$DTU`. Indicates that the gene as a whole shows significant change in isoform ratios. `= (Genes[parent_id, DTU])`.
* `abundA` and `abundB` - (num) Either the mean or the sum of the abundances across replicates, depending on what is indicated by `Parameters$use_sums`. This is the value used for the tests.
* `stdevA` and `stdevB` - (num) The standard deviation of the abundance across the replicates.
* `log2FC` - (num) log2 of fold-change of transcript abundance: `= log2(abundB / abundA)`. `abundA` and `abundB` are not internally normalised for library size, so care must be taken when interpreting `log2FC`.
* `totalA` and `totalB` - (num) The total abundance for the gene: `totalA= sum(transcripts[parent_id, abundA])`.
* `propA` and `propB` - (num) The proportion of the gene expression owed to this transcript. `propA= abundA / totalA`.
* `Dprop` - (num) The difference in the proportion of the transcript between the two conditions (effect size). `= (probB - propA)`.
* `pval` - (num) The proportion equality test P-value for the transcript.
* `pval_corr` - (num) Multiple testing corrected p-value. `= p.adjust(pval, Parameters$correction)`.
* `quant_p_median` - (num) Median of corrected p-value across quantification bootstraps (not available in `lean` mode).
* `quant_p_min` - (num) Minimum observed (corrected) p-value across quantification bootstraps (not available in `lean` mode).
* `quant_p_max` - (num) Maximum observed (corrected) p-value across quantification bootstraps (not available in `lean` mode).
* `quant_Dprop_mean` - (num) Mean effect size across quantification bootstraps (not available in `lean` mode).
* `quant_Dprop_stdev` - (num) Standard deviation of effect size across quantification bootstraps (not available in `lean` mode).
* `quant_Dprop_min` - (num) Minimum observed effect size across quantification bootstraps (not available in `lean` mode).
* `quant_Dprop_max` - (num) Maximum observed effect size across quantification bootstraps (not available in `lean` mode).
* `quant_na_freq` - (num) Fraction of quantification iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
* `quant_dtu_freq` - (num) Fraction of quantification iterations that support a positive `DTU` classification.
* `rep_p_median` - (num) Median of corrected p-value across replication bootstraps (not available in `lean` mode).
* `rep_p_min` - (num) Minimum observed (corrected) p-value across replication bootstraps (not available in `lean` mode).
* `rep_p_max` - (num) Maximum observed (corrected) p-value across replication bootstraps (not available in `lean` mode).
* `rep_Dprop_mean` - (num) Mean effect size across replication bootstraps (not available in `lean` mode).
* `rep_Dprop_stdev` - (num) Standard deviation of effect size across replication bootstraps (not available in `lean` mode).
* `rep_Dprop_min` - (num) Minimum observed effect size across replication bootstraps (not available in `lean` mode).
* `rep_Dprop_max` - (num) Maximum observed effect size across replication bootstraps (not available in `lean` mode).
* `rep_na_freq` - (num) Fraction of replication iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
* `rep_dtu_freq` - (num) Fraction of replication iterations that support a positive `DTU` classification.

**Note:** The fields reporting on the bootstraps will not be shown when bootstrapping is disabled.


## Abundances
 
`Abundances` is a list of two `data.table`s, one for each condition. 
Each transcript is represented by a single abundance value per replicate, so if bootstrapped data was used, these values are the means across iterations. If plain abundances were provided as input, then `Abundances` essentially contains the input data. These abundances are included in the output because they are required for some of RATs' plotting options.


```r
# Elements of Abundances.
print( names(mydtu$Abundances) )
```

```
## [1] "condA" "condB"
```

1. `condA` - (num) The transcript abundances in the first condition.
2. `condB` - (num) The transcript abundances in the second condition.


```r
# Abundance table for first condition.
print( head(mydtu$Abundances[[1]]) )
```

```
##     V1  V2 target_id parent_id
## 1: 100 115    ALLA:1      ALLA
## 2:  40  30    ALLA:2      ALLA
## 3: 300 400    D2TE_a      D2TE
## 4: 400 500    D2TE_b      D2TE
## 5: 250 270     1D2TU      D2TU
## 6: 250 300     2D2TU      D2TU
```

***

# Quick results

For your convenience, RATs provides a few functions to give you a quick summary of the run. They all follow the same style.

These reports should not be seen as a substitute for the detailed RATs output.

## Summary of DTU

The `dtu_summary()` function lists the total number of genes and transcripts for each of 3 categories:

* DTU:  There is significant change in isoform ratios (in terms of both effect size and statistical significance).
* non-DTU:  No significant change.
* ineligible: Ineligible for testing. Genes/transcripts with read count below the set threshold, or where the gene has only one known transcript.


```r
# A tally of the outcome.
print( dtu_summary(mydtu) )
```

```
##                           category tally
## 1            DTU genes (gene test)     6
## 2        non-DTU genes (gene test)     3
## 3     ineligible genes (gene test)     3
## 4         DTU genes (transc. test)     6
## 5     non-DTU genes (transc. test)     3
## 6  ineligible genes (transc. test)     3
## 7           DTU genes (both tests)     6
## 8       non-DTU genes (both tests)     3
## 9    ineligible genes (both tests)     3
## 10                 DTU transcripts    14
## 11             non-DTU transcripts     7
## 12          ineligible transcripts     5
```

The `get_dtu_ids()` function lists the coresponding identifiers per category.
The ID lists obtained are ordered by effect size (`Dprop`).


```r
# Gene and transcript IDs corresponding to the tally above.
ids <- get_dtu_ids(mydtu)
print( names(ids) )
```

```
##  [1] "DTU genes (gene test)"           "non-DTU genes (gene test)"      
##  [3] "ineligible genes (gene test)"    "DTU genes (transc. test)"       
##  [5] "non-DTU genes (transc. test)"    "ineligible genes (transc. test)"
##  [7] "DTU genes (both tests)"          "non-DTU genes (both tests)"     
##  [9] "ineligible genes (both tests)"   "DTU transcripts"                
## [11] "non-DTU transcripts"             "ineligible transcripts"
```

```r
print( ids )
```

```
## $`DTU genes (gene test)`
## [1] "XSW"  "MIX"  "SW"   "LC"   "ALLA" "D2TU"
## 
## $`non-DTU genes (gene test)`
## [1] "D2TE" "SAME" "FAKE"
## 
## $`ineligible genes (gene test)`
## [1] "LONE" "SOLO" "NN"  
## 
## $`DTU genes (transc. test)`
## [1] "XSW"  "MIX"  "SW"   "LC"   "ALLA" "D2TU"
## 
## $`non-DTU genes (transc. test)`
## [1] "D2TE" "SAME" "FAKE"
## 
## $`ineligible genes (transc. test)`
## [1] "LONE" "SOLO" "NN"  
## 
## $`DTU genes (both tests)`
## [1] "XSW"  "MIX"  "SW"   "LC"   "ALLA" "D2TU"
## 
## $`non-DTU genes (both tests)`
## [1] "D2TE" "SAME" "FAKE"
## 
## $`ineligible genes (both tests)`
## [1] "LONE" "SOLO" "NN"  
## 
## $`DTU transcripts`
##  [1] "XSW:one" "XSW:two" "MIX.b"   "SW1"     "SW2"     "MIX.a"   "LC1"    
##  [8] "LC2"     "ALLA:1"  "ALLA:2"  "MIX.ab"  "2D2TU"   "1D2TU"   "MIX.l2" 
## 
## $`non-DTU transcripts`
## [1] "D2TE_b" "D2TE_a" "SAME_1" "SAME_2" "FAKE-2" "FAKE-1" "MIX.l1"
## 
## $`ineligible transcripts`
## [1] "LONE.a" "MIX.n"  "SOLO.1" "NNa"    "NNb"
```

## Summary of isoform switching

Isoform switching is a subset of DTU. Primary isoform switching is often considered the most likely 
type of DTU to have an effect. The following two functions summarise the extent of isoform switching 
in the results:


```r
# A tally of genes switching isoform ranks.
print( dtu_switch_summary(mydtu) )
```

```
##                            category genes
## 1        Primary switch (gene test)     4
## 2    Non-primary switch (gene test)     1
## 3     Primary switch (transc. test)     4
## 4 Non-primary switch (transc. test)     1
## 5       Primary switch (both tests)     4
## 6   Non-primary switch (both tests)     1
```

```r
# The gene IDs displaying isoform switching.
ids <- get_switch_ids(mydtu)
print( names(ids) )
```

```
## [1] "Primary switch (gene test)"        "Non-primary switch (gene test)"   
## [3] "Primary switch (transc. test)"     "Non-primary switch (transc. test)"
## [5] "Primary switch (both tests)"       "Non-primary switch (both tests)"
```

## Summary of DTU plurality

In case you want to know how many isoforms are affected per gene.


```r
# A tally of genes switching isoform ranks.
print( dtu_plurality_summary(mydtu) )
```

```
##   isof_affected num_of_genes
## 1             2            5
## 2             4            1
```

```r
# The gene IDs displaying isoform switching.
ids <- get_plurality_ids(mydtu)
```

These give you the number and IDs of genes in which 2, 3, etc... isoforms show DTU.


***

# Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged to help us understand the problem. Also, be sure 
to **read the vignette(s)**, and browse/search the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./figs/rats_logo.png)


