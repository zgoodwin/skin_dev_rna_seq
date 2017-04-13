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
# Updated: 04/12/2017
#
###############################################################################

########################
#       PACKAGES       #
########################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(stringr, reshape, ggplot2, optparse)
source('multiplot.R')

#######################
#      FUNCTIONS      #
#######################

readHitTable <- function(file){
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
  colnames(df)[1] <- "gene_Id"
  return(df)
}

#######################
# EXECUTED STATEMENTS #
#######################

option_list <- list(
  make_option(c("-i", "--hitTableFile"), type = "character", default = NULL,
              help = "Table containing differentially expressed genes, fold changes and p-values.", 
              metavar = "character"),
  make_option(c("-p", "--plots"), type = "character", default = NULL,
              help = "Directory for plots", 
              metavar = "character"),
  make_option(c("-l","--geneList"), type = "character", default = NULL,
              help = "Comma-separated list of genes to be plotted",
              metavar = "character"))


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Parse command line arguments
# args <- commandArgs(trailingOnly = TRUE)
hitTableFile <- opt$hitTableFile
geneList <- strsplit(opt$geneList, ",")[[1]]

# Check command line arguments
if ( is.null(hitTableFile) ){
  stop("The file with differential gene expression counts is not specified.",call.=FALSE)
} else if (length(geneList) == 0){
  stop("No genes specified.", call.=FALSE)
}

# Read in the list of differentially expressed genes.
hits <- readHitTable(hitTableFile)

## Select the following fields: 
#   gene name (col 1)
#   logFC (col 2)
#   ExperimentName (col 8)
index <- which(hits$gene %in% geneList)
values <- hits[index, c(1,2,8)]
values <- transform(values, RESULT = colsplit(ExperimentName, split = "_", names = c("Timepoint", "Genotype")))
values$ExperimentName = NULL
colnames(values)[c(3,4)] = c("Timepoint", "Genotype")

## Set the name of the plot
plot_name = paste(geneList, collapse="_")
pdf(paste0(opt$plots,"/goi_v_time_",plot_name,".pdf"))

## Draw the plot
ggplot(data=values, aes(x = Timepoint, y = logFC, colour = Genotype)) +
       geom_point() +
       geom_line(aes(group = Genotype)) +
       facet_grid(~ gene_Id) +
       xlab("Timepoint") +
       ylab("Log Fold-Change") +
       theme_bw() + 
       theme(axis.text.x = element_text(colour="black", angle=45, hjust=1, vjust=1),
             strip.text.x = element_text(face="bold"))
dev.off()
################### END ###################