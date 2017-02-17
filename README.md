# Analysis of gene expression levels in developing mice with a skin-specific enhancer knockout.

The code in this repository is an analysis that I performed in R to analyze gene expression levels in developing mouse skin. I sought to study the function of an [enhancer](https://en.wikipedia.org/wiki/Enhancer_(genetics)) that I believe is important for the formation of the mouse epidermis. To study the function of this enhancer, my lab made mice that had intact copies of the enhancer (wt), mice where only one copy of the enhancer was deleted (het), and mice where both copies of the enhancer were deleted (mut). Since this enhancer is usually active early in mouse development, I measured genome-wide gene expression levels at the E14.5, E15.5, E16.5, E17.5 and Newborn (NB) time points. Expression levels for all genes in the mouse genome were measured using [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) for each time point. 

I determined the fold change in gene expression for hom mice relative to the wt mice and for the het mice relative to wt mice, at each developmental time point. To do this, I used bioconductor's [limma voom](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) package. 

**DISCLAIMER:** This is unpublished data. As such, the gene names in the count data have been de-identified, and the original BAM files from which the gene counts were derived will not be provided.

The **analysis** directory contains two scripts: one that performs the differential expression analysis for read count data (genotype_analysis_script.R) and one that plots gene expressio values versus time (goi_vs_time.R). Below are the help files for each script: 
```
Rscript genotype_analysis_script.R -h 
Usage: genotype_analysis_script.R [options]


Options:
	-c CHARACTER, --countFile=CHARACTER
		Read count file (see ./data/data_formats.md for format of this file)

	-k CHARACTER, --keyFile=CHARACTER
		Key File (see ./data/data_formats.md for the format of this file)

	-r CHARACTER, --results=CHARACTER
		Directory for output files

	-p CHARACTER, --plots=CHARACTER
		Directory for plots

	-h, --help
		Show this help message and exit
```

```
Rscript goi_vs_time.R -h
Plot gene expression values versus development time.
	Rscript goi_vs_time.R <differenially_expressed_genes.tsv> <gene1> ... <gene N>
```

## Input Data

* **./data/all_counts.csv -** Transcript counts for the genes in all all mice (all timepoints and genotypes)
* **./data/sample_key.csv -** Metadata for each sample, specifying the timepoint, genotype and unique ID for each mouse.

The ./data/data_formats.md file contains definitions for each column of each table.

## Output Data

* **./data/all_diffexp_genes.csv -** List of differentially expressed genes (relative to wt) for all time points.

For each of the 5 time points, `genotype_analysis_script.R` prints out a table of differential gene expression values (log2-fold change) for each EDC gene. The ./data/data_formats.md file contains definitions for each column of this table.


## Plot info

The ./plots/ directory contains diagnostic plots drawn by the genotype_analysis_script program. Below are some brief descriptions on how to interpret them. More information can be found on the [limma voom](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) documentation.

* **./plots/mds_plots.pdf -** Diplays the mulidimensional scaling (MDS) for each time point.
* **./plots/mean_var_plots.pdf -** Plots the gene expression couns versus the variance in gene exprssion counts for each gene (the curve should look like a sideways "S")
* **./plots/weight_plots.pdf -** Quality weights for each of the experiments. The lower the weight, the lower the sample quality.

To plot expression levels for your favorite gene in the hom mice \(relative to wild type\) for each of the five developmental time points, you can use the goi\_vs\_time script.

## To do's:

- [x] Finish script to plot gene of interest versus developmental time point.
- [x] Add plots
- [x] Re-organize code
- [ ] Fix plot label ordering
