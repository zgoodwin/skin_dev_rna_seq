# Analysis of gene expression levels in developing mice with a skin-specific enhancer knockout.

The code in this document reprepents an analysis that I performed in R to analyze gene expression levels in developing mouse skin. We measured gene expresion levels using [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) at the E14.5, E15.5, E16.5, E17.5 and Newborn (NB) time points. Furthermore, these mice were genetically modified via CRISPR-cas to be deficient for an [enhancer](https://en.wikipedia.org/wiki/Enhancer_(genetics)) that we believe is important for the formation of the epidermis. RNA-seq for each time point was measured for mice that were homozygous knockout (hom) and wild-type (wt) for the enhancer.

To do this, I used bioconductor's [limma voom](https://bioconductor.org/packages/release/bioc/html/limma.html) package for RNA-seq data. Limma also provides a useful set of tools for performing statistical tests on count data derived from RNA-seq experiments.

*DISCLAIMER:* This is unpublished data. As such, the gene names in the count data have been de-identified, and the original BAM files from which the gene counts were derived will not be provided.


You can run the entire analysis as follows:
```
Rscript ./analysis/genotype_analysis_script.R
```

## Plot info

The ./plots/ directory contains diagnostic plots drawn by the genotype_analysis_script program. Below are some brief descriptions on how to interpret them. More information can be found on the [limma voom](https://bioconductor.org/packages/release/bioc/html/limma.html) documentation.

* ./plots/mds_plots.pdf - Diplays the mulidimensional scaling (MDS) for each time point.
* ./plots/mean_var_plots.pdf - Plots the gene expression couns versus the variance in gene exprssion counts for each gene (the curve should look like a sideways "S")
* ./plots/weight_plots.pdf - Quality weights for each of the experiments. The lower the weight, the lower the sample quality.

To plot expression levels for your favorite gene in the hom mice \(relative to wild type\) for each of the five developmental time points, you can use the goi\_vs\_time script, which is not yet implemented.

*To do's:*

- [ ] Finish script to plot gene of interest versus developmental time point.
- [ ] Fix plot ordering
- [x] Add plots
- [x] Re-organize code
