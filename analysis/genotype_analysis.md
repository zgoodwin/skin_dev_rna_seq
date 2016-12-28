# Analysis of genotypes in Skin-Enhancer Knockout mice
Zane Goodwin  
December 25, 2016  
## Background
Enhancers are critical for increasing gene expression in a tissue-specific manner. My lab generated a mouse line where an enhancer required for epithelial development has been removed via genetic manipulation. As we weren't sure which genes would be affected by removing this enhancer, we sought to measure the expression of all genes in the the genomes of mice with the enhancer intact (wild type or "WT") and in the genmoes of mice with the enhancer partially intact (heterozygotes or "het") or completely removed (mutants or "mut"). By analyzing the diffrerences in gene expression between the "het" and "wt" and between the "mut" and "wt" mice, then we will gain a clearer picture of how this enhancer regulates gene expression in mice. Furthermore, we would like to see how gene epxression levels are changing between embryonic developmental time points (E14.5, E15.5, E16.5, E17.5 and newborn (NB)) when the skin barrier is established in mice.

**Disclaimer:** Since this is unpublished work, I am going to be vague about which enhancer was removed, how the enhancer was removed, and about which genes were affected when the enhancer was removed. As such, I have replaced the mouse gene names with ID numbers.  

## Methods
We measured gene expression in mice using a technique called RNA-Seq, where RNA molecules from mouse skin cells are sequenced using high-throughput sequencing technology. For each of the 5 time poins, We sequenced RNA from:

* 2 wild-type mice
* 2 heterozygous (het) mice
* 2 mutant (mut) mice

Overall, we sequenced rna from 30 mice, mapped the RNA sequences to the mouse reference genome using Bowtie and counted the number of reads matching to each gene using the `blah` from from R's `bioconductor` package. To control for biological variation and batch effects, we identified differentially-expressed genes using the `limma` pacakge. At the end of the day, we will end up with a list of differentially expressed genes for each of the time points.

## Setting the environment and getting the count data
First, I'm going to load the `limma` and `edgeR` packages, which contain tools necessary to perform the analysis. These packages contain the statistical methods that will normalize read counts in a way that makes it possible to identify differentially expressed genes. I'll also be using the stringR package for a bit of string manipulation.

```r
library(limma)
library(edgeR)
library(stringr)

# File locations
COUNT_FILE = "../data/all_counts.csv"
SAMPLE_KEY_FILE = "../data/sample_key_obs.csv"
DESIGN_MATRICES = "../design_matrices/"
RESULTS_DIR = "../results/"

# Read in the data
sampleKey = read.csv(SAMPLE_KEY_FILE, header=T)
counts = read.csv(COUNT_FILE, header=T, row.names=1)
```

Now, I will make a few DGE list objects for each mouse age in order to compare gene expression levels at each developmental time point.


```r
setup_dge = function(counts_table, search_string){
  counts = counts_table[,grep(search_string,colnames(counts))]
  dge = DGEList(counts = counts)
  dge = calcNormFactors(dge)
  
  # Get rid of all genes with less than 10 counts (counts per million < 2)
  keep = rowSums(cpm(dge)>2) >= 3
  dge = dge[keep, , keep.lib.sizes=FALSE]
   
  return(dge)
}

e14_5_dge = setup_dge(counts, "E14.5")
e15_5_dge = setup_dge(counts, "E15.5")
e16_5_dge = setup_dge(counts, "E16.5")
e17_5_dge = setup_dge(counts, "E17.5")
nb_dge = setup_dge(counts, "NB")
```


## Sample quality assessment

An important part of any RNA-seq experiment is to estimate the batch effects that might inadvertently cause differences in gene expression. One way to do this is to look at the MDS of RNA-seq data from each mouse. For this experiment, the best possible case is that all of the heterozyote, WT and mutant mice have similar variances in gene expression levels. If they do, then they will cluster together on the MDS plot. 

![](genotype_analysis_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The MDS analsyis seems to show that there are overall differences in gene expression between the f4 and f5 mice for the newborns (NB) and the mice at E16.5. This is a less-than-ideal scenario, and probably means that for the E14.5, E15.5 and E17.5 mice, any differences that we observe in gene expression have more to do with overall differences in gene expression between the f4 and f5 replicate mice than differences between the het, mut and wt mice.

Given this, we will want to down-weight low-quality samples, which can be done in the next section. 

## Library normalization
This part constructs the design matrices and performs the library normalization necessary to determine log-fold-change. Also, given that we observed some outlier samples in the MDS analysis, we will want to down-weight the outlier samples, in this case, X, Y and Z. 

```r
get_design_matrix = function(sample_list, search_string, reference){
  sample = sample_list[grep(search_string,sample_list$Age),]
  genotype = relevel(sample$Genotype, ref=reference)
  design = model.matrix(~genotype)
  return(design)
}

plot_mean_variance = function(voom_object, title){
  plot(voom_object$voom.xy, 
       main = paste(title,":Mean-variance trend"), 
       xlab = "log2( count size + 0.5 )", 
       ylab = "Sqrt( standard deviation )")
  lines(voom_object$voom.line$x, voom_object$voom.line$y, type = "l", col="red")
}

e14_5_dm = get_design_matrix(sampleKey, "E14.5", "wt")
e15_5_dm = get_design_matrix(sampleKey, "E15.5", "wt")
e16_5_dm = get_design_matrix(sampleKey, "E16.5", "wt")
e17_5_dm = get_design_matrix(sampleKey, "E17.5", "wt")
nb_dm = get_design_matrix(sampleKey, "NB", "wt")

# Normalize read counts
e14_5_norm = voom(e14_5_dge, e14_5_dm, plot=FALSE, save.plot = TRUE)
e15_5_norm = voom(e15_5_dge, e15_5_dm, plot=FALSE, save.plot = TRUE)
e16_5_norm = voom(e16_5_dge, e16_5_dm, plot=FALSE, save.plot = TRUE)
e17_5_norm = voom(e17_5_dge, e17_5_dm, plot=FALSE, save.plot = TRUE)
nb_norm = voom(nb_dge, nb_dm, plot=FALSE, save.plot = TRUE)

# Check the mean-variance trend to visualize the distribution of read counts
par(mfrow=c(2,3))
plot_mean_variance(e14_5_norm, "E14.5")
plot_mean_variance(e15_5_norm, "E15.5")
plot_mean_variance(e16_5_norm, "E16.5")
plot_mean_variance(e17_5_norm, "E17.5")
plot_mean_variance(nb_norm, "NB")
```

![](genotype_analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


```r
library(ggplot2)
plot_weights = function(voom_object, title){
  df = data.frame(samples = rownames(voom_object$targets), 
                  wts = voom_object$sample.weights)
  
  text_color = "black"
  p = ggplot(data=df, aes(x=samples, y=wts)) + 
      ggtitle(title) + 
      xlab("") + 
      ylab("Sample Weight") + 
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1, colour=text_color),
            axis.text.y = element_text(colour=text_color)) +
      geom_bar(stat="identity") + 
      geom_hline(yintercept=1.0, colour = "red", linetype="dashed")
  return(p)
}

# Down-weight outlier samples 
e14_5_wts = voomWithQualityWeights(e14_5_dge, e14_5_dm, normalization = "none", plot=FALSE)
e15_5_wts = voomWithQualityWeights(e15_5_dge, e15_5_dm, normalization = "none", plot=FALSE)
e16_5_wts = voomWithQualityWeights(e16_5_dge, e16_5_dm, normalization = "none", plot=FALSE)
e17_5_wts = voomWithQualityWeights(e17_5_dge, e17_5_dm, normalization = "none", plot=FALSE)
nb_wts = voomWithQualityWeights(nb_dge, nb_dm, normalization = "none", plot=FALSE)

# Visualize weights
source("./multiplot.R")
multiplot(plot_weights(e14_5_wts, "E14.5"),
         plot_weights(e15_5_wts, "E15.5"),
         plot_weights(e16_5_wts, "E16.5"),
         plot_weights(e17_5_wts, "E17.5"),
         plot_weights(nb_wts, "NB"), cols=3)
```

```
## Loading required package: grid
```

![](genotype_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Calling differentially-expressed genes
Fit linear model and perform an empirical Bayes test to model gene expression data. The purpose of this step is to fit a linear model to all the gene expression levels for each sample, then do an empirical bayes test to distingush whether or not each gene is differentially expressed.

```r
# Fit bayes models to the weighted samples
e14_5_vfit = lmFit(e14_5_wts) %>% eBayes(.)
e15_5_vfit = lmFit(e15_5_wts) %>% eBayes(.)
e16_5_vfit = lmFit(e16_5_wts) %>% eBayes(.)
e17_5_vfit = lmFit(e17_5_wts) %>% eBayes(.)
nb_vfit = lmFit(nb_wts) %>% eBayes(.)

diffexp_14_5 = topTable(e14_5_vfit, adjust="BH", number = 10, sort.by = "F")
diffexp_15_5 = topTable(e15_5_vfit, adjust="BH", number = 10, sort.by = "F")
diffexp_16_5 = topTable(e16_5_vfit, adjust="BH", number = 10, sort.by = "F")
diffexp_17_5 = topTable(e17_5_vfit, adjust="BH", number = 10, sort.by = "F")
diffexp_nb = topTable(nb_vfit, adjust="BH", number = 10, sort.by = "F")
```

## Visualization of differentially expressed genes



