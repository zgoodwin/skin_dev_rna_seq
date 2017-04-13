###############################################################################
# Program: differentialExpression.R
# Project: Analysis of RNA-seq data from the developing mouse epidermis
#          For mice lacking a skin-specific enhancer
# Description: Performs the limma voom RNA-seq analysis for genome-wide
#               read-count data.
#
# Author: Zane Goodwin
# Date Written: 01/11/2017
# Last Updated: 04/12/2017
# 
###############################################################################

########################
#       GLOBALS        #
########################

NUM_HITS <- Inf # Return all differentially expressed genes
MULTIPLE_TEST_CORRECTION <- "BH"


########################
#       PACKAGES       #
########################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(limma, edgeR, stringr, 
               reshape, ggplot2, optparse, ggrepel)

#######################
#      FUNCTIONS      #
#######################

setupDge <- function(countsTable){
  # Creates a dge object.
  #
  # Args:
  #   counts_table: Table of read counts
  #
  # Returns: A voom dge object for the RNA-seq experiment
  #
  thisDge <- DGEList(counts = countsTable)
  thisDge <- calcNormFactors(thisDge)

  ## Get rid of all genes with less than 10 counts (counts per million < 2)
  #  If you don't do this, the distribution of read counts will be highly skewed!
  keep <- rowSums(cpm(thisDge)>2) >= 3
  filteredDge <- thisDge[keep, , keep.lib.sizes=FALSE]

  return(filteredDge)
}

getMdsPlot <- function(dgeObject, title, sampleList){
  # Generate an MDS plot for count data wrapped in a DGE object
  #
  # Args:
  #   dgeObject: A limma dge object representing an RNA-seq experiment
  #   title: Key name to access an experiment from counts_table
  #   sampleList: Data frame of sample metadata
  #
  # Returns: an R plot object
  #

  samples = sampleList[, 1]
  plot = plotMDS(dgeObject, labels=samples,
                 main=paste(title, "MDS", sep=" "),
                 xlim=c(-1.5,1.5),
                 ylim=c(-1.5,1.5))
  return(plot)
}

plotMeanVariance <- function(voomObject, title){
  # Produces a plot that shows the read counts for each gene (x-axis)
  #   count standard deviation for each gene.
  #
  # Args:
  #   voomObject: Sample metadata table
  #   title: The title of the plot
  #
  # Returns: an R plot object
  #

  p <- plot(voomObject$voom.xy,
           main = paste(title,":Mean-variance trend"),
           xlab = "log2( count size + 0.5 )",
           ylab = "Sqrt( standard deviation )")

  lines(voomObject$voom.line$x, voomObject$voom.line$y, type = "l", col="red")
  return(p)
}

plotWeights <- function(voomObject, title){
  # Plot the weights for each RNA-seq experiment.
  #
  # Args:
  #   voomObject: Normalized limma voom object.
  #   title: Plot title
  #
  # This part puts the samples and weights into a data frame so ggplot can
  # read them.
  #
  # Returns: an R plot object

  df <- data.frame(samples = rownames(voomObject$targets),
                  wts = voomObject$sample.weights)

  x <- barplot(height = df[, 2], 
          xaxt = "n",
          main = paste0(title, " - Quality Weights"))
  abline(h=1, lty=5, col = "red")
  labs <- as.vector(df[, 1])
  text(cex = 1, x = x, y = -0.02, labs, xpd = TRUE, srt = 45, pos = 2)
  
  return()
}

plotVolcano <- function(hitTable, title){
  # Create volcano plots to show differentially expressed genes,
  #   Show only the top 10 differentially expressed genes
  #
  # Args:
  #   hitTable: Table containing log fold change and P-values
  #   title: Desired plot title
  #
  # Returns: a ggplot plot object
  
  textSize <- 15
  labelSize <- 15
  pointLabelSize <- 3
  
  ## Label the genes with the most differential expression
  top <- head(hitTable, n=10)
  geneNames <- rownames(top)
  hitTable$pointTitle <- rownames(hitTable)
  hitTable$pointTitle[which(!(rownames(hitTable) %in% geneNames))] <- ""
  hitTable$significant <- "0"
  hitTable$significant[which(hitTable$P.Value < 0.001)] <- "1"
  
  ## Draw the plot object
  p <- ggplot(data = hitTable, aes(x = logFC, y = (-1*log10(P.Value)))) + 
    scale_color_manual(values = c("black", "red")) +
    geom_point(aes(color = significant)) + 
    theme_bw() + 
    geom_label_repel(aes(label = pointTitle), size = pointLabelSize) + 
    ggtitle(title) + 
    ylab("-log10 p-value") +
    xlab("Log Fold-Change") + 
    theme(axis.text.x = element_text(size = textSize),
          axis.text.y = element_text(size = textSize),
          axis.title.x = element_text(size = labelSize),
          axis.title.y = element_text(size = labelSize))
  return(p)
}

#######################
# EXECUTED STATEMENTS #
#######################

## Parse command line arguments
option_list <- list(
  make_option(c("-c", "--controlCountFile"), type = "character", default = NULL, 
              help = "Count file for control mice (should contain at least 2 replicates)", 
              metavar = "character"),
  make_option(c("-t", "--testCountFile"), type = "character", default = NULL,
              help = "Count file for test mice (either WT or MUT)",
              metavar = "character"),
  make_option(c("-r", "--results"), type = "character", default = NULL,
              help = "Directory for output files", 
              metavar = "character"),
  make_option(c("-p", "--plots"), type = "character", default = NULL,
              help = "Directory for plots", 
              metavar = "character"),
  make_option(c("-a", "--experimentName"), type = "character", default = NULL,
              help ="Name for this experiment"))


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Check command line arguments
if ( is.null(opt$controlCountFile) ){
  print_help(opt_parser)
  stop("The control count file is not specified.", call.=FALSE)
} else if( is.null(opt$testCountFile) ){
  print_help(opt_parser)
  stop("The test count file is not specified.", call.=FALSE)
} else if( is.null(opt$result) ){
  print_help(opt_parser)
  stop("The result directory file is not specified.", call.=FALSE)
} else if( is.null(opt$plots) ){
  print_help(opt_parser)
  stop("The plot directory file is not specified.", call.=FALSE)
}

## Load the read counts for the control mice and knockout mice
controlCounts <- read.csv(opt$controlCountFile, header=T, row.names=1)
testCounts <- read.csv(opt$testCountFile, header=T, row.names=1)
allCounts <- cbind(controlCounts, testCounts)
experimentName <- opt$experimentName

numControlReplicates <- ncol(controlCounts)
numTestReplicates <- ncol(testCounts)

## Set up the design matrix
#   This is a key step that defines how the gene expression 
#   values should be compared
designMatrix <- data.frame(sampleName = c(colnames(controlCounts),
                                           colnames(testCounts)),
                           groupId = c(rep("REF", numControlReplicates),
                                        rep("TEST", numTestReplicates)))
print("Building design matrix")
design <- model.matrix(~groupId, designMatrix)

## Set up the data structure for the DGE object
print("Building DGE object")
dge <- setupDge(allCounts)

## Perform voom normalization
print("Voom normalization")
v <- voom(dge, design, plot=FALSE, save.plot = TRUE)

## Obtain quality weights
print("Weight calculation")
weights <- voomWithQualityWeights(dge, design,
                                 normalization = "none", plot=FALSE)

## Get the list of the most differentially expressed genes
print("Hits table")
vfit = lmFit(weights) %>% eBayes(.)
hits = topTable(vfit, 
                adjust = MULTIPLE_TEST_CORRECTION, 
                number = NUM_HITS, 
                sort.by = "p")
hits$ExperimentName = rep(experimentName, nrow(hits))


## Draw Quality control plots
print("Drawing plots...")
pdf(file = paste0(opt$plots, "/", experimentName, "_quality_plots.pdf"), width = 8, height = 8)
par(mfrow = c(2,2))
plotMeanVariance(v, experimentName)
plotWeights(weights, experimentName)
getMdsPlot(v, experimentName, designMatrix)
dev.off()


## Draw the volcano plot
print("Drawing the Volcano plot ...")
pdf(file = paste0(opt$plots, "/", experimentName, "_volcano.pdf"), width = 5, height = 5)
plot(plotVolcano(hits, experimentName))
dev.off()


## Write the results table
print("Writing hit file...")
write.csv(hits, file = paste0(opt$results, "/",experimentName, "_results.csv"),
          row.names = T, quote = F)

################### END ###################