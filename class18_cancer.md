Class 18: Cancer Mutations Mini-Project
================
Kira

## 1. Exploring a Cancer Sequencing Data Portal: Skin Cutaneous Melanoma

### Discussion \#1:

> Q1. How many cancer samples are included in the dataset?

There are 448 samples in this dataset.

> Q2. Which is the most mutated gene?

TTN is the most commonly mutated gene (80.0%).

> Q3. Which is the most common treatment undergone by patients?

Radiation is the most common treatment undergone by patients.

## 3. Generating Mutational Matrices and Visualizing Mutational Profiles:

Mutational matrices are the first step for mutational signature analysis
and correspond to a helpful data type, as it contains no protected
information.

``` r
mm_cumel = read.delim('http://www.tinyurl.com/skinmatrix')
dim(mm_cumel)
```

    [1]  96 440

``` r
# Install MutationalPatterns package
if (!require("MutationalPatterns")){
BiocManager::install("MutationalPatterns")
}
```

    Loading required package: MutationalPatterns

    Loading required package: GenomicRanges

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min

    Loading required package: S4Vectors


    Attaching package: 'S4Vectors'

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomeInfoDb

    Loading required package: NMF

    Loading required package: registry

    Loading required package: rngtools

    Loading required package: cluster

    NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 7/8

      To enable shared memory capabilities, try: install.extras('
    NMF
    ')


    Attaching package: 'NMF'

    The following object is masked from 'package:S4Vectors':

        nrun

``` r
# Generate mutational profiles (4 random samples)
library(MutationalPatterns)
set.seed(11111) # fixing the seed for random number generation

samples_to_plot = sample(1:ncol(mm_cumel),4) # selecting 4 random samples
plot_96_profile(mm_cumel[,samples_to_plot], condensed = T)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-2-1.png)

``` r
# Generate mutational profiles (top 4 mutated samples and top 4 less mutated)
mutations_in_samples = colSums(mm_cumel)
mutations_in_samples = sort(mutations_in_samples, decreasing = T)
samples_to_plot = names(mutations_in_samples)[1:4]
plot_96_profile(mm_cumel[,samples_to_plot], condensed = T)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-2-2.png)

``` r
mutations_in_samples = sort(mutations_in_samples, decreasing = F)
samples_to_plot = names(mutations_in_samples)[1:4]
plot_96_profile(mm_cumel[,samples_to_plot], condensed = T)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-2-3.png)

``` r
# Generate average mutational profiles
relative_mutational_profile = apply(mm_cumel, 2, prop.table) 
# obtained relative mutational matrix
average_mutational_profile = rowMeans(relative_mutational_profile)
average_mutational_profile = data.frame(average_mutational_profile)
plot_96_profile(average_mutational_profile, condensed = T)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-2-4.png)

There are additional ways to visualize the mutations that may be useful:

``` r
hist(log10(mutations_in_samples))
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
boxplot(mutations_in_samples)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-4-1.png)

## 5. Assigning Reference Mutational Signature:

``` r
# Mutational signature assignment
cosmic_signatures = get_known_signatures(source = 'COSMIC_v3.2')
fit_res = fit_to_signatures(mm_cumel, cosmic_signatures)

# Top contributing signatures
contributions = fit_res$contribution

top_contributing_signatures_abs = rowMeans(contributions)
top_contributing_signatures_abs = sort(top_contributing_signatures_abs,
                                       decreasing = T)[1:4]

## Top 4 contributing signatures (absolute values)
top_contributing_signatures_abs
```

        SBS7a     SBS7b     SBS38      SBS4 
    366.97614 340.91011 204.44450  99.49106 

``` r
relative_contributions = apply(contributions,2,prop.table)
top_contributing_signatures_rel = rowMeans(relative_contributions)
top_contributing_signatures_rel = sort(top_contributing_signatures_rel,
                                       decreasing = T)[1:4]

## Top 4 contributing signatures (relative values)
top_contributing_signatures_rel
```

         SBS7b      SBS7a      SBS38       SBS4 
    0.26336351 0.26019455 0.10885595 0.07240978 

``` r
# Mutational signature assignment strict
fit_res_strict = fit_to_signatures_strict(mm_cumel, cosmic_signatures)
fit_res_strict = fit_res_strict$fit_res
contributions_strict = fit_res_strict$contribution
```

## 6. Visualizing Mutational Signature Assignment Results:

To visualize the mutational signature assignment results, we will use
the default visualizations available in the `MutationalPatterns`
package.

``` r
# Visualization of signature assignment results (fit_to_signatures)
set.seed(11111)
samples_to_plot = sample(1:ncol(mm_cumel),4)

plot_contribution(contributions[,samples_to_plot], mode = "absolute")
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-8-1.png)

``` r
plot_contribution(contributions[,samples_to_plot], mode = "relative")
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-9-1.png)

``` r
plot_contribution_heatmap(contributions, cluster_samples = F)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
# Visualization of signature assignment results (strict)
plot_contribution(contributions_strict[,samples_to_plot], mode = "absolute")
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
plot_contribution(contributions_strict[,samples_to_plot], mode = "relative")
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
plot_contribution_heatmap(contributions_strict, cluster_samples = F)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-13-1.png)

To check the cosine similarity of the reconstruction for some specific
samples, we can use the following visualization from the
`MutationalPatterns` R package:

``` r
# Cosine similarity reconstruction vs. original mutational profile (fit_to_signatures)
set.seed(11111)
samples_to_plot = sample(1:ncol(mm_cumel),4)

plot_original_vs_reconstructed(mm_cumel[,samples_to_plot],
                               fit_res$reconstructed[,samples_to_plot], 
                               y_intercept = 0.90)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-14-1.png)

``` r
# Cosine similarity reconstruction vs. original mutational profile (strict)
plot_original_vs_reconstructed(mm_cumel[,samples_to_plot],
                               fit_res_strict$reconstructed[,samples_to_plot], 
                               y_intercept = 0.90)
```

![](class18_cancer_files/figure-commonmark/unnamed-chunk-15-1.png)

### Discussion \#2:

> Q4. Which is the etiology of the top absolute contributing signature
> for liver cancer?

The etiology of the top absolute contributing signature for liver cancer
is aristocholic acid exposure.

> Q5. Which is the most prominent mutational context for the top
> contributing signature in skin cancer?

The most prominent mutational context for the top contributing signature
(SBS7a) is C\>T in skin cutaneous melanoma.

> Q6. (True or False): The etiology of the top contributing signature
> for lung cancer corresponds to an endogenous cellular mechanism.

False. The top contributing signature for lung cancer corresponds to
tobacco smoking (SBS4).

> Q7. (True or False): SBS4 is one of the most common signatures found
> in lung cancer and is associated with tobacco smoking

True, SBS4 is common in lung cancer and is associated with smoking.

> Q8. (True or False): SBS7d is one of the most common signatures in
> skin cancer and is associated with UV light exposure and high numbers
> of C\>T mutations.

This is false. SBS7a and SBS7b are common signatures, but SBS7d is not
listed.
