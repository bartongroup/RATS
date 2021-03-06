---
title: "RATs: Plots"
author: "Kimon Froussios"
date: "04 APR 2017"
output: 
  html_document: 
    keep_md: yes
    theme: readable
    toc: yes
    toc_float: TRUE
vignette: >
  %\VignetteIndexEntry{RATs 4: Plots}
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
                  varname= "condition", verbose= FALSE,
                  description="Comparison of two conditions using a simulated sleuth object 
                    for the purposes of the tutorial. Simulated using built-in functionality 
                    of RATs.")
```


***


# Visualisation of results


The RATs output object provides a host of information and we encourage users to familiarize themselves with it. 
But a good plot is worth a thousand numbers.


## Isoform abundances for a given gene

This function allows you to visualise what's going on in any particular gene. Both the absolute counts and the relative 
proportions are plotted for each transcript. This is a very useful function for inspecting a gene of interest. It enables 
quick visual evaluation of the dispersion of the replicate measurements, the magnitude of the proportion change, the 
presence of outliers, and the consistency among the replicates.

By default, these plots can be quite colourful (and possibly ugly) as their aim is to highlight which isoforms are predicted to be DTU. 
Options are provided to change which information layers are encoded and what colours are used.

There are several styles for this plot, depending on your preferences. The default plot format includes the most information:


```r
# Grouping by condition (DEAFULT):
#   plot_gene(mydtu, "MIX6")
plot_gene(mydtu, "MIX6", style="bycondition")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

The top two facets show the absolute abundances (counts) of the isoforms, the bottom two show the relative abundances (proportions).
The left two facets refer to one condition, the right two ones to the other condition. The boxplots describe the abundance measurements
of each isoform across replicates. As the proportions of isoforms are inter-dependent and add up to 100%, the coloured lines connect
the isoform abundances measured in each replicate. By default, the boxes' fill-colour encodes the transcript level DTU prediction for
each isoform: In this example, isoforms `.c1` and `.c2` change significantly, isoforms `.c3`, `.c4` and `.nc` do not change significantly 
and isoform `.d` could not be tested (for very narrow boxes, the fill colour is difficult to see).

This structure is great for seeing changes in the isoform distribution profile of a gene. A cleaner versions of the plot can also be obtained:


```r
# Grouping by condition (minimalist):
plot_gene(mydtu, "MIX6", style="linesonly")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

When there are many lines or many isoforms, it may be preferable to group abundances by isoform, making individual comparisons easier:


```r
# Grouping by isoform:
plot_gene(mydtu, "MIX6", style="byisoform")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)


### Customisation of the gene plot

#### Change information layers

The `fillby`, `colourby` and `shapeby` parameters can be respectively used to control which information layers
are encoded as fill, line/point colour, and point shape. Possible values are `c("isoform", "condition", "DTU", "none", "replicate")`.
Not all options are available in all styles, in which case they will be silently ignored.


```r
# Change the encoded information.
plot_gene(mydtu, "MIX6", style="bycondition", fillby="isoform")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
plot_gene(mydtu, "MIX6", style="byisoform", colourby="DTU", shapeby="replicate")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
# For a less busy look, any of the information layers can be disabled.
plot_gene(mydtu, "MIX6", style="byisoform", colourby="none", shapeby="none")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)

#### Change colour code


```r
# Colour codes can be customised by specifying new values for
# condcolvec, replcolvec, isofcolvec, dtucolvec and nonecol.
plot_gene(mydtu, "MIX6", style="bycondition", fillby="condition", condcolvec=c("magenta", "cyan"))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)


## Plots of overall run

Our simulated dataset is too small to properly demonstrate what these plots typically look like.
So each one is accompanied by a static image of the same plot created with a real and much larger dataset.

Possibly the most common plot in differential expression is the volcano plot, which plots the effect size against 
the statistical significance. As it is difficult to define a single p-value and a single effect size at the gene level,
the volcano can only be plotted at the transcript level.


```r
# Proportion change VS significance.
plot_overview(mydtu, type="volcano")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

This is what it looks like on a larger dataset:
![Dprop VS sig](./fig/volcano.png)

The next command plots the largest change in proportion seen within each gene, against the number of genes showing 
such change. This is a way to inspect what effect sizes are present in the data. As an additional layer of information,
they are colour-coded by their DTU result.


```r
# Distribution of maximum proportion change.
plot_overview(mydtu, type="maxdprop")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

This is what it looks like on a larger dataset:
![Max Dprop](./fig/maxdprop.png)


### Interactive plots

If you prefer picking points from a plot than sorting through tables, the volcano plot is also available through 
a `shiny` app, that brings up the relevant abundance changes plot for any point in the volcano plot.

1. By hovering over points on the volcano plot in the app, you can see the respective transcript identifier(s). 
2. Clicking will pull up information on the effect size, significance and confidence of the point(s), as well as 
the respective isoform abundance changes plot for the point nearest to the click.


```r
# Start the interactive volcano plot.
plot_shiny_volcano(mydtu)
```

This is what it looks like for the example data (remember that the emulated data example has very few transcripts).

![Transc Conf VS DTU](./fig/shiny_screenshot.png)

You will need to close down the app to return to your R terminal.


## Plot customisation

You can save any of the plots as a `ggplot2` object and use [ggplot2](http://ggplot2.org) manipulations on it, such as changing the axis scales.
Other `ggplot2` customisations include the axis tick marks, axis values, labels, titles, colours... Consult the [ggplot2](http://ggplot2.org)
documentation for more help on these.


```r
library(ggplot2)

myplot <- plot_overview(mydtu, "volcano")
myplot  # display
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
# Change title. 
myplot2 <- myplot + ggtitle("MY EPIC TITLE")
myplot2
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-2.png)


***


# Contact information

The `rats` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kimon Froussios, Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or **ask for assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/Rats/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged to help us understand the problem. Also, be sure 
to **read the vignette(s)**, and browse/search the support forum before posting a new issue, in case your question is already answered there.

Enjoy!

![](./fig/rats.png)


