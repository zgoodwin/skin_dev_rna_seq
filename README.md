# Analysis of gene expression levels in developing mice with a skin-specific enhancer knockout.

The code in this repository is an analysis that I performed in R to analyze gene expression levels in developing mouse skin. I sought to study the function of an [enhancer](https://en.wikipedia.org/wiki/Enhancer_(genetics)) that I believe is important for the formation of the mouse epidermis. To study the function of this enhancer, my lab made mice that had intact copies of the enhancer (wt), mice where only one copy of the enhancer was deleted (het), and mice where both copies of the enhancer were deleted (mut). Since this enhancer is usually active early in mouse development, I measured genome-wide gene expression levels at the E14.5, E15.5, E16.5, E17.5 and Newborn (NB) time points. Expression levels for all genes in the mouse genome were measured using [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) for each time point. 

I determined the fold change in gene expression for hom mice relative to the wt mice and for the het mice relative to wt mice, at each developmental time point. To do this, I used bioconductor's [limma voom](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) package. 

**DISCLAIMER:** This is unpublished data. As such, the gene names in the count data have been de-identified, and the original BAM files from which the gene counts were derived cannot be provided.

The **analysis** directory contains two scripts: one that performs the differential expression analysis for read count data (differentialExpression.R) and one that plots gene expression values versus time (genesVsTime.R). Below are the help files for each script: 

## differentialExpression.R

```
Rscript differentialExpression.R
Loading required package: pacman
Usage: differentialExpression.R [options]


Options:
	-c CHARACTER, --controlCountFile=CHARACTER
		Count file for control mice (should contain at least 2 replicates)

	-t CHARACTER, --testCountFile=CHARACTER
		Count file for test mice (either HET or MUT, should also contain at least 2 replicates)

	-r CHARACTER, --results=CHARACTER
		Directory for output files

	-p CHARACTER, --plots=CHARACTER
		Directory for plots

	-a EXPERIMENTNAME, --experimentName=EXPERIMENTNAME
		Name for this experiment

	-h, --help
		Show this help message and exit

```

## genesVsTime.R

```
Rscript genesVsTime.R -h
Loading required package: pacman
Usage: genesVsTime.R [options]


Options:
	-i CHARACTER, --hitTableFile=CHARACTER
		Table containing differentially expressed genes, fold changes and p-values.

	-p CHARACTER, --plots=CHARACTER
		Directory for plots

	-l CHARACTER, --geneList=CHARACTER
		Comma-separated list of genes to be plotted

	-h, --help
		Show this help message and exit


```

There is also a driver script (written in bash) called `analyzeTimePoints.sh`, which runs the differential expression analysis for each time point, combines the results, and then shows the change in gene expression over time for two genes of interest. Tables of genes containing differential expression data are stored in the `results` directory, and volcano plots, plus plots showing change in gene expression over time are in the `plots` directory.

## Input Data

* **A count file in CSV format (Control)-** Transcript counts for the genes in mice from each developmental time point, for the control genotype (wild-type). Columns correspond to biological replicates.
* **A count file in CSV format (Test)-** Transcript counts for the genes in mice from each developmental time point, for the test genotype (heterozygote or mutant). Again, columns correspond to biological replicates.

## Output Data

* **./data/all_diffexp_genes.csv -** List of differentially expressed genes (relative to wt) for all time points.

For each of the 5 time points, `differentialExpression.R` prints out a table of differential gene expression values (log2-fold change) for each EDC gene. The format of the output file is as follows:

1. **logFC -**  Average log-fold change in expression for the test group (heterozygote or mutant) versus the controls (wild type).
3. **AveExpr -**  Average expression level for the gene across all samples.
4. **F -**  F-statistic.
5. **P.Value -**  Raw P-value.
6. **adj.P.Val -**  FDR-corrected P-value.
7. **B -** B-statistic.
8. **ExperimentName -**  Name of the experiment. In this case, it's the timepoint and the genotype of the mouse.


## Plot info

The ./plots/ directory contains diagnostic plots drawn by the `differentialExpression.R` program. Below are some brief descriptionsfor each type of plot. More information can be found on the [limma voom](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) documentation.

# Data plots

* **goi_v_time plots -** Shows gene expression versus time in the developing mouse. GOI stands for "gene of interest"
* **Volcano plots -** Shows gene expression versus P-value to visualize differentially expressed genes for each experiment. Significant genes are labeled with text.

# Diagnostic plots
* **./plots/mds_plots.pdf -** Diplays the mulidimensional scaling (MDS) for a given time point, to assess inter-sample variation.
* **./plots/mean_var_plots.pdf -** Plots the gene expression couns versus the variance in gene exprssion counts for each gene (the curve should look like a sideways "S")
* **./plots/weight_plots.pdf -** Quality weights for each of the experiments. The lower the weight, the lower the sample quality.*

