# Analysis of expression levels in developing mice with a skin-specific enhancer knockout.

The code in this document reprepents an analysis that I performed in R to analyze gene expression levels in developing mouse skin. We measured gene expresion levels using [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) at the E14.5, E15.5, E16.5, E17.5 and Newborn (NB) time points. Furthermore, these mice were genetically modified via CRISPR-cas to be deficient for an [enhancer](https://en.wikipedia.org/wiki/Enhancer_(genetics)) that we believe is important for the formation of the epidermis. RNA-seq for each time point was measured for mice that were homozygous knockout (hom) and wild-type (wt) for the enhancer.

To do this, I used bioconductor's [limma voom](https://bioconductor.org/packages/release/bioc/html/limma.html) package for RNA-seq data. Limma also provides a useful set of tools for performing statistical tests on 

You can run the entire analysis as follows:
```
Rscript genotype_analysis_script.R
```

To plot expression levels for your favorite gene in the hom mice \(relative to wild type\) for each of the five developmental time points, you can use the following command:
```
Rscript whatever
```
