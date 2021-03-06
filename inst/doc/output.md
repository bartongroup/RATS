---
title: "RATs: Raw Output"
author: "Kimon Froussios"
date: "19 APR 2017"
output: 
  html_document: 
    keep_md: yes
    theme: readable
    toc: yes
    toc_float: TRUE
vignette: >
  %\VignetteIndexEntry{RATs 3: Raw Output}
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}

---



***

Set up an example.


```r
library(rats)

# Simulate some data.
simdat <- sim_sleuth_data(cnames = c("controls", "patients")) 
# For convenience let's assign the contents of the list to separate variables.
myslo <- simdat$slo
myannot <- simdat$annot

# Call DTU
mydtu <- call_DTU(annot = myannot, slo = myslo, name_A = "controls", name_B = "patients", 
                  varname= "condition", verbose= FALSE, dprop_thresh=0.1, qboot=FALSE, rboot=FALSE,
                  description="Comparison of two conditions using a simulated sleuth object 
                    for the purposes of the tutorial. Simulated using built-in functionality 
                    of RATs.")
```

We lowered `dprop_thresh` and bypassed the bootstrapping to increase the number of DTU positive 
classifications, to better illustrate the summary fucntions below for the given artificial data.

***


# Quick results

For your convenience, RATs provides a few functions to give you a quick summary of your results.
However, we do recommend you become familiar with the actual results structure and content, so that you
can judge the quality of the DTU calls and trace the reasons behind the classification of each gene or transcript.

## Summary of DTU

The `dtu_summary()` function lists the total number of genes and transcripts for each of 3 categories:

* DTU:  There is significant change in isoform ratios.
* non-DTU:  No significant change.
* NA:  Not applicable. Genes/transcripts with abundance below the noise threshold, or where the gene has only one transcript.


```r
# A really simple tally of the outcome.
print( dtu_summary(mydtu) )
```

```
##        DTU genes (gene test)    non-DTU genes (gene test) 
##                            2                            1 
##         NA genes (gene test)     DTU genes (transc. test) 
##                            7                            2 
## non-DTU genes (transc. test)      NA genes (transc. test) 
##                            8                            0 
##       DTU genes (both tests)   non-DTU genes (both tests) 
##                            2                            1 
##        NA genes (both tests)              DTU transcripts 
##                            0                            5 
##          non-DTU transcripts               NA transcripts 
##                            5                           11
```

Notice that **three types of results are given for genes**. That's because RATs uses *two* separate methods of identifying DTU:
The gene-level test works on the isoforms of each gene as a set, whereas the transcript-level test works on each isoform individually
and then aggregates the isoform results to form a gene-level result. Therefore, as of `v.0.4.2`, `dtu_summury()`
will display the DTU genes according to either method as well as the intersection of the two.

The `get_dtu_ids()` function lists the actual identifiers per category, instead of the numbers in each category.
As of `v0.4.2` it uses the same category names as `dtu_summury()` for consistency and clarity.


```r
# Gene and transcript IDs corresponding to the tally above.
ids <- get_dtu_ids(mydtu)
print( ids )
```

```
## $`DTU genes (gene test)`
## [1] "MIX6" "CC"  
## 
## $`non-DTU genes (gene test)`
## [1] "NN"
## 
## $`NA genes (gene test)`
## [1] "LC"   "1A1N" "1B1C" "1D1C" "ALLA" "ALLB" "NIB" 
## 
## $`DTU genes (transc. test)`
## [1] "MIX6" "CC"  
## 
## $`non-DTU genes (transc. test)`
## [1] "LC"   "NN"   "1A1N" "1B1C" "1D1C" "ALLA" "ALLB" "NIB" 
## 
## $`NA genes (transc. test)`
## character(0)
## 
## $`DTU genes (both tests)`
## [1] "MIX6" "CC"  
## 
## $`non-DTU genes (both tests)`
## [1] "NN"
## 
## $`NA genes (both tests)`
## character(0)
## 
## $`DTU transcripts`
## [1] "MIX6.c1" "MIX6.c2" "MIX6.c4" "CC_a"    "CC_b"   
## 
## $`non-DTU transcripts`
## [1] "LC2"     "MIX6.c3" "2NN"     "1NN"     "MIX6.nc"
## 
## $`NA transcripts`
##  [1] "LC1"      "1A1N-2"   "1B1C.1"   "1B1C.2"   "1D1C:one" "1D1C:two"
##  [7] "MIX6.d"   "ALLA1"    "ALLB1"    "ALLB2"    "NIB.1"
```

The ID lists obtained with `get_dtu_ids()` are *ordered by effect size* (`Dprop`).


## Summary of isoform switching

Isoform switching is a subset of DTU. Primary isoform switching is often considered the most likely 
type of DTU to have an effect. The following two functions summarise the extent of isoform switching 
in the results:


```r
# A tally of genes switching isoform ranks.
print( dtu_switch_summary(mydtu) )
```

```
##        Primary switch (gene test)    Non-primary switch (gene test) 
##                                 1                                 1 
##     Primary switch (transc. test) Non-primary switch (transc. test) 
##                                 1                                 1 
##       Primary switch (both tests)   Non-primary switch (both tests) 
##                                 1                                 1
```

```r
# The gene IDs displaying isoform switching.
ids <- get_switch_ids(mydtu)
print( ids )
```

```
## $`Primary switch (gene test)`
## [1] MIX6
## Levels: CC MIX6
## 
## $`Non-primary switch (gene test)`
## [1] MIX6
## Levels: CC MIX6
## 
## $`Primary switch (transc. test)`
## [1] MIX6
## Levels: CC MIX6
## 
## $`Non-primary switch (transc. test)`
## [1] MIX6
## Levels: CC MIX6
## 
## $`Primary switch (both tests)`
## [1] MIX6
## Levels: CC MIX6
## 
## $`Non-primary switch (both tests)`
## [1] MIX6
## Levels: CC MIX6
```

Again, three versions of the gene results are shown, covering the respective DTU methods in RATs.
`MIX6` is demonstrating a switch in both primary and tertiary isoforms.

## Summary of DTU plurality

In case you want to know how many isoforms are affected per gene.


```r
# A tally of genes switching isoform ranks.
print( dtu_plurality_summary(mydtu) )
```

```
## 2 3 
## 1 1
```

```r
# The gene IDs displaying isoform switching.
ids <- get_plurality_ids(mydtu)
print( ids )
```

```
## $`2`
## [1] CC
## Levels: 1A1N 1B1C 1D1C ALLA ALLB CC LC MIX6 NIB NN
## 
## $`3`
## [1] MIX6
## Levels: 1A1N 1B1C 1D1C ALLA ALLB CC LC MIX6 NIB NN
```

The categories are named by the number of isoforms affected. There is one gene (`CC`) that has `2` isoforms affected
and one gene (`MIX6`) that has `3` isoforms affected.

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
## [13] "tests"               "p_thresh"            "abund_thresh"       
## [16] "dprop_thresh"        "quant_reprod_thresh" "quant_boot"         
## [19] "quant_bootnum"       "rep_reprod_thresh"   "rep_boot"           
## [22] "rep_bootnum"         "rep_reprod_as_crit"
```

1. `description` - (str) Free-text description of the run. It is useful to record data sources, annotation source and version, experimental parameters...
2. `time` - (str) Date and time for the run. This does not represent the exact time the run started.
3. `rats_version` - (str) The version of RATs.
4. `R_version` - (str) The version of R (including OS architecture).
5. `var_name` - (str) The value passed to the `varname` parameter.
6. `cond_A` & `cond_B` - (str) The values passed to the `name_A` and `name_B` parameters.
7. `data_type` - (str) The format of the input data.
8. `num_replic_A` & `num_replic_B` - (int) The number of samples in each condition.
9. `num_genes` - (int) The number of genes in the provided annotation.
10. `num_transc` - (int) The number of transcripts in the provided annotation.
11. `tests` - (str) The value passed to the `testmode` parameter.
12. `p_thresh` - (num) The value passed to the `p_thresh` parameter.
13. `abund_thresh` - (num) The value passed to the `abund_thresh` parameter.
14. `dprop_thresh` - (num) The value passed to the `dprop_thresh` parameter.
15. `quant_reprod_thresh` - (num) The value passed to the `qrep_thresh` parameter.
16. `quant_boot` - (bool) The value passed to the `qboots` parameter.
17. `quant_bootnum` - (int) The value passed to the `qbootnum` parameter.
18. `rep_reprod_thresh` - (num) The value passed to the `rrep_thresh` parameter.
19. `rep_boots` - (bool) The value passed to the `rboot` parameter.
20. `rep_bootnum` - (int) The number of replicate bootstrapping iterations (usually  M*N, where M and N are the number of samples in the two conditions).
21. `rep_reprod_as_crit` - (bool) The value passed to the `rrep_as_crit` parameter.

**Note:** If bootstraps are disabled, the bootstrap-related fields may have `NA` values despite any values that were passed to their respective parameters.


## Genes

`Genes` is a [data.table](https://cran.r-project.org/web/packages/data.table/) with many fields, listing 
results at the gene level.  For your convenience, the respective aggregated transcript-level DTU calls are 
also included here (defined as at least one isoform being called DTU individually).

The maximum likelihood-based G-test of independence is used to compare the set of isoform ratios between the two conditions. 


```r
# Genes table's fields.
print( names(mydtu$Genes) )
```

```
##  [1] "parent_id"     "elig"          "sig"           "elig_fx"      
##  [5] "DTU"           "transc_DTU"    "known_transc"  "detect_transc"
##  [9] "elig_transc"   "pval"          "pval_corr"
```

The first few columns show the result of each decision step. The remaining columns list the values based on which the decisions were made.
Pseudo-code formulas are shown to help understand how the different fields interact when making decisions.

1. `parent_id` - (str) Gene identifier.
2. `elig` - (bool) Eligible for testing. Whether the gene met the pre-filtering criteria. `= (elig_transc >= 2)`.
3. `sig` - (bool) Statistically significant. `= (pvalAB_corr < Parameters$p_thresh) & (pvalBA_corr < Parameters$p_thresh)`.
4. `elig_fx` - (bool) Eligible effect size. Whether at least one of the isoforms meets the effect size criterion. `= any(Transcript[parent_id, elig_fx])`.
5. `quant_reprod` - (bool) Quantification reproducible.
For positive DTU, `= (quant_dtu_freq >= Parameters$quant_reprod_thresh)`, for non-DTU `= (quant_dtu_freq <= 1 - Parameters$quant_reprod_thresh)`.
6. `rep_reprod` - (bool) Replication reproducible.
For positive DTU, `= (rep_dtu_freq >= Parameters$rep_reprod_thresh)`, for non-DTU `= (rep_dtu_freq <= 1 - Parameters$rep_reprod_thresh)`.
7. `DTU` - (bool) Differential Transcript Usage. `= (sig & elig_fx & quant_reprod)`, or if `rrep_as_crit==TRUE` then `= (sig & elig_fx & quant_reprod & rep_reprod)`.
8. `transc_DTU` - (bool) Aggregated from `Transcripts$DTU`: At least one isoform was individually reported as DTU. `= any(Transcript[parent_id, DTU])`.
9. `known_transc` - (int) Number of known trascripts for this gene according to the given annotation.
10. `detect_trancs` - (int) Number of detected (expressed) transcripts in the given dataset.
11. `elig_transc` - (int) Number of eligible transcripts, aggregated from `Transcripts$elig`.
12. `pval` - (num) G test of independence p-value for the isoform ratios.
13. `pval_corr` - (num) Multiple testing corrected p-value. `pval_corr= p.adjust(pval, Parameters$correction)`.
14. `quant_p_mean` - (num) Mean p-value across quantification bootstraps.
15. `quant_p_stdev` - (num) Standard deviation of p-values across quantification bootstraps.
16. `quant_p_min` - (num) Minimum observed p-value across quantification bootstraps.
17. `quant_p_max` - (num) Maximum observed p-value across quantification bootstraps.
18. `quant_na_freq` - (num) Fraction of quantification iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
19. `quant_dtu_freq` - (num) Fraction of replicate iterations that support a positive `DTU` classification.
20. `rep_p_mean` - (num) Mean p-value across replication bootstraps.
21. `rep_p_stdev` - (num) Standard deviation of p-values across replication bootstraps.
22. `rep_p_min` - (num) Minimum observed p-value across replication bootstraps.
23. `rep_p_max` - (num) Maximum observed p-value across replication bootstraps.
24. `rep_na_freq` - (num) Fraction of replication iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
25. `rep_dtu_freq` - (num) Fraction of replication iterations that support a positive `DTU` classification.

**Note:** The fields reporting on the bootstraps will not be shown when bootstrapping is disabled.


## Transcripts

`Transcripts` is a [data.table](https://cran.r-project.org/web/packages/data.table/) with many fields, listing 
results at the transcript level. For your convenience, the respective gene-level DTU calls are also included here.

The maximum likelihood-based G-test of independence is used to compare each given isoform's proportions in the two conditions.


```r
# Transcripts table's fields.
print( names(mydtu$Transcripts) )
```

```
##  [1] "target_id" "parent_id" "elig_xp"   "elig"      "sig"      
##  [6] "elig_fx"   "DTU"       "gene_DTU"  "meanA"     "meanB"    
## [11] "stdevA"    "stdevB"    "sumA"      "sumB"      "log2FC"   
## [16] "totalA"    "totalB"    "propA"     "propB"     "Dprop"    
## [21] "pval"      "pval_corr"
```

1. `target_id` - (str) Transcript identifier.
2. `parent_id` - (str) Gene identifier.
3. `elig_xp` - (bool) Eligible expression level. `= (meanA >= Parameters$abund_thresh | meanB >= Parameters$abund_thresh)`.
4. `elig` - (bool) Eligible for testing (meets the noise threshold and at least one other isoform is expressed). `= (elig_xp & totalA != 0 & totalB != 0 & (sumA != totalA | sumB != totalB))`.
5. `sig` - (bool) Statistically significant. `= (pval_corr < Parameters$p_thresh)`.
6. `elig_fx` - (bool) Eligible effect size. Proxy for biological significance. `= (Dprop > Parameters$dprop_thresh)`.
7. `quant_reprod` - (bool) Quantification reproducible.
For positive DTU, `= (quant_dtu_freq >= Parameters$quant_reprod_thresh)`, for non-DTU `= (quant_dtu_freq <= 1 - Parameters$quant_reprod_thresh)`.
8. `rep_reprod` - (bool) Replication reproducible.
For positive DTU, `= (rep_dtu_freq >= Parameters$rep_reprod_thresh)`, for non-DTU `= (rep_dtu_freq <= 1 - Parameters$rep_reprod_thresh)`.
9. `DTU` - (bool) The Transcript is Differentially Used. `= (sig & elig_fx & quant_reprod)`, or if `rrep_as_crit==TRUE` then `= (sig & elig_fx & quant_reprod & rep_reprod)`.
10. `gene_DTU` - (bool) Expanded from `Genes$DTU`. Indicates that the gene as a whole shows significant change in isoform ratios. `= (Genes[parent_id, DTU])`.
11. `meanA` and `meanB` - (num) The mean abundance across replicates, for each condition.
12. `stdevA` and `stdevB` - (num) The standard deviation of the abundance across the replicates.
13. `sumA` and `sumB` - (num) The sum of the abundances across replicates. This is the value used for the tests, so that replication level informs the significance.
14. `FC` - (num) Fold-change of transcript abundance: `= (sumB / sumA)`.
15. `totalA` and `totalB` - (num) The total abundance for the gene: `totalA= sum(transcripts[parent_id, sumA])`.
16. `propA` and `propB` - (num) The proportion of the gene expression owed to this transcript. `propA= sumA / totalA`.
17. `Dprop` - (num) The difference in the proportion of the transcript between the two conditions (effect size). `= (probB - propA)`.
18. `pval` - (num) The proportion equality test P-value for the transcript.
19. `pval_corr` - (num) Multiple testing corrected p-value. `= p.adjust(pval, Parameters$correction)`.
20. `quant_p_mean` - (num) Mean p-value across quantification bootstraps.
21. `quant_p_stdev` - (num) Standard deviation of p-values across quantification bootstraps.
22. `quant_p_min` - (num) Minimum observed p-value across quantification bootstraps.
23. `quant_p_max` - (num) Maximum observed p-value across quantification bootstraps.
24. `quant_na_freq` - (num) Fraction of quantification iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
25. `quant_dtu_freq` - (num) Fraction of quantification iterations that support a positive `DTU` classification.
26. `rep_p_mean` - (num) Mean p-value across replication bootstraps.
27. `rep_p_stdev` - (num) Standard deviation of p-values across replication bootstraps.
28. `rep_p_min` - (num) Minimum observed p-value across replication bootstraps.
29. `rep_p_max` - (num) Maximum observed p-value across replication bootstraps.
30. `rep_na_freq` - (num) Fraction of replication iterations in which `DTU` could not be determined (not meeting pre-filtering criteria).
31. `rep_dtu_freq` - (num) Fraction of replication iterations that support a positive `DTU` classification.

**Note:** The fields reporting on the bootstraps will not be shown when bootstrapping is disabled.


## Abundances
 
`Abundances` is a list of two `data.table`s, one for each condition. Each transcript is represented by a single abundance value per replicate, 
so if bootstrapped data was used, these values are the means across iterations. If plain abundances were provided as input, then `Abundances` 
essentially contains the input data. These abundances are included in the output because they are required for some of RATs' plotting options.



```r
# Elements of ReplicateData
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
##          V1   V2 target_id parent_id
## 1: 19.66667 20.5    1A1N-2      1A1N
## 2:  0.00000  0.0    1B1C.1      1B1C
## 3: 52.33333 53.5    1B1C.2      1B1C
## 4:  0.00000  0.0  1D1C:one      1D1C
## 5: 76.00000 78.0  1D1C:two      1D1C
## 6: 50.00000 45.0     ALLA1      ALLA
```

***


# Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged to help us understand the problem. Also, be sure 
to **read the vignette(s)**, and browse/search the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./fig/rats.png)


