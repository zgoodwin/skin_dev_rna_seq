###############################################################################
# Program: goi_vs_time.R
# 
# Project: Analysis of RNA-seq data from the developing mouse epidermis
#
# Description: Plots log fold change in gene expression over time for user-specified
#           genes. (plots are stored in the ../plots/ directory)
#
# Author: Zane Goodwin
# Date: 01/16/2017
#
###############################################################################

########################
#       PACKAGES       #
########################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(stringr, reshape, ggplot2)
source('multiplot.R')

#######################
#      FUNCTIONS      #
#######################

read_dge = function(file){
  # Reads in log-fold change data from a limma voom analysis 
  #   (see genotype_analysis_script.R) 
  #
  # Args:
  #   file: String of the file name for the gene expression table
  #
  # Returns: A data frame for log-fold changes from the RNA-seq experiment 
  #          specified in search_string
  # 
  df = read.csv(file, header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = c("gene_Id",
                   "het_v_wt_FC",
                   "mut_v_wt_FC",
                   "AveExpr",
                   "F", "P.Value", "adj.P.Val", "gene", "timepoint")
  return(df)
}

#######################
# EXECUTED STATEMENTS #
#######################

# dge = read_dge("../data/all_diffexp_genes.csv")
# gene_list = commandArgs(trailingOnly = TRUE)
args = commandArgs(trailingOnly = TRUE)
dgeFile = args[1]
gene_list = args[2:length(args)]

# Check command line arguments
if (is.na(dgeFile)){
  stop("The file with differential gene expression counts is not specified.",call.=FALSE)
} else if (length(gene_list) == 0){
  stop("No genes specified.", call.=FALSE)
} else if (dgeFile == "-h"){
  stop("Plot gene expression values versus development time.\n\tRscript goi_vs_time.R <differenially_expressed_genes.tsv> <gene1> ... <gene N>", call.=FALSE)
}

# Read in the list of differentially expressed genes.
dge = read_dge(dgeFile)

# Select the following fields: 
#   het_v_wt_FC (col 2)
#   mut_v_wt_FC (col 3)
#   gene        (col 8)
#   timepoint   (col 9)
dge_select = dge[which(dge$gene %in% gene_list),
                 c(2,3,8,9)]

# Melt the log-FC values into long format to make them
# compatible with ggplot
mdf = melt(dge_select, id.vars=3:4, 
           measure.vars=1:2, variable_name = "genotype")

# Set the name of the plot
plot_name = paste(c("../plots/time_v_fc",
                    paste(gene_list, collapse="_")),
                  collapse="_")

pdf(paste0(plot_name,".pdf"))
ggplot(data=mdf, aes(x=timepoint, y = value, colour=genotype)) +
       geom_point() +
       geom_line(aes(group=genotype)) +
       facet_grid(~ gene) +
       xlab("Timepoint") +
       ylab("Log Fold-Change") +
       theme(axis.text.x = element_text(colour="black", angle=45, hjust=1, vjust=1),
             strip.text.x = element_text(face="bold"))
dev.off()
################### END ###################