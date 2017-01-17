# all_counts.csv

This file contains the read counts for each gene. The first column is a list of identifiers for genes in the mouse genome.

The remainder of the columns are identifers for each RNA-seq experiment. They are formatted as follows:

```
  [genotype]_[time point]_[litter number]
```

# all_diffexp_genes.csv

This file is the output from genotype_analysis_script.R. The column names are as follows:

1. **genotypehet -**  Average log-fold change in expression for the heterozygous mice, relative to wt (wild-type).
2. **genotypemut -**  Average log-fold change in expression for the homozygous mutant mice, relative to wt (wild-type).
3. **AveExpr -**  Average expression level for the gene across all samples.
4. **F -**  F-statistic.
5. **P.Value -**  Raw P-value.
6. **adj.P.Val -**  FDR-corrected P-value.
7. **gene -**  Gene identifier.
8. **timepoint -**  Developmental time point.

# sample_key.csv

This file contains information about each RNA-seq experiment (in long format) and it is used to construct the deisign matrices.

Column names:

1. **Genotype -**  Mouse genotype
2. **Age -**  Mouse developmental time point
3. **Group -**  Litter number

