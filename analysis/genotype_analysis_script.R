###############################################################################
# Program: genotype_analysis.R
# Project: Analysis of RNA-seq data from the developing mouse epidermis
# Description: Performs the limma voom RNA-seq analysis for genome-wide
#               read-count data.
#
# Author: Zane Goodwin
# Date: 01/11/2017
#
###############################################################################

########################
#       PACKAGES       #
########################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(limma, edgeR, stringr, reshape, ggplot2, optparse)

#######################
#      FUNCTIONS      #
#######################

setup_dge = function(counts_table, search_string){
  # Creates a dge object. 
  #
  # Args:
  #   counts_table: Table of read counts 
  #   search_string: Key name to access an experiment from counts_table
  #
  # Returns: A voom dge object for the RNA-seq experiment specified in
  #   search_string
  # 
  
  counts = counts_table[,grep(search_string,colnames(counts))]
  dge = DGEList(counts = counts)
  dge = calcNormFactors(dge)
  
  ## Get rid of all genes with less than 10 counts (counts per million < 2)
  #  If you don't do this, the distribution of read counts will be highly skewed!
  keep = rowSums(cpm(dge)>2) >= 3
  dge = dge[keep, , keep.lib.sizes=FALSE]
  
  return(dge)
}

get_mds_plot = function(dge_object, search_string, sample_list){
  # Generate an MDS plot for count data wrapped in a DGE object
  # 
  # Args:
  #   dge_object: A limma dge object representing an RNA-seq experiment
  #   search_string: Key name to access an experiment from counts_table
  #   sample_list: Data frame of sample metadata
  #
  # Returns: an R plot object
  # 
  
  sample = sample_list[grep(search_string,sample_list$Age),]
  plot = plotMDS(dge_object, labels=paste(sample$Genotype, sample$Group),
                 main=paste(search_string, "MDS", sep=" "),
                 xlim=c(-1.5,1.5),
                 ylim=c(-1.5,1.5))  
  return(plot)
}

get_design_matrix = function(sample_list, search_string, reference){
  # Create a design matrix for fitting linear models to the count data.
  # 
  # Args:
  #   dge_object: Sample metadata table 
  #   search_string: Key name to access an experiment from counts_table
  #   sample_list: Name of the reference genotype
  #
  # Returns: a text object representing a design matrix for a linear model
  # 
  
  sample = sample_list[grep(search_string,sample_list$Age),]
  genotype = relevel(sample$Genotype, ref=reference)
  design = model.matrix(~genotype)
  return(design)
}

plot_mean_variance = function(voom_object, title){
  # Produces a plot that shows the read counts for each gene (x-axis)
  #   count standard deviation for each gene.
  # 
  # Args:
  #   voom_object: Sample metadata table 
  #   title: The title of the plot
  #
  # Returns: an R plot object
  # 
  
  p = plot(voom_object$voom.xy, 
       main = paste(title,":Mean-variance trend"), 
       xlab = "log2( count size + 0.5 )", 
       ylab = "Sqrt( standard deviation )")
  
  lines(voom_object$voom.line$x, voom_object$voom.line$y, type = "l", col="red")
  return(p)
}

plot_weights = function(voom_object, title){
  # Plot the weights for each RNA-seq experiment.
  #
  # Args:
  #   voom_object: Normalized limma voom object.
  #   title: Plot title
  #
  # This part puts the samples and weights into a data frame so ggplot can 
  # read them.
  #
  # Returns: an R plot object 
  
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

#######################
# EXECUTED STATEMENTS #
#######################

## Parse command line arguments
option_list = list(
  make_option(c("-c", "--countFile"), type="character", default=NULL, 
              help="Read count file", metavar="character"),
  make_option(c("-k", "--keyFile"), type="character", default=NULL, 
              help="Key File", metavar="character"),
  make_option(c("-r", "--results"), type="character", default=NULL,
              help="Directory for output files", metavar="character"),
  make_option(c("-p", "--plots"), type="character", default=NULL,
              help="Directory for plots", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check command line arguments
if (is.null(opt$countFile)){
  print_help(opt_parser)
  stop("The count file is not specified.", call.=FALSE)
} else if(is.null(opt$keyFile)){
  print_help(opt_parser)
  stop("The Key file is not specified.", call.=FALSE)
} else if(is.null(opt$result)){
  print_help(opt_parser)
  stop("The result directory file is not specified.", call.=FALSE)
} else if(is.null(opt$plots)){
  print_help(opt_parser)
  stop("The plot directory file is not specified.", call.=FALSE)
}

## Load the sample key, which contains metadata about each RNA-seq run 
#  (format is specified in the data/data_format.md file)
sampleKey = read.csv(opt$keyFile, header=T)

## Load the read counts for each gene and for each experiment 
#  (format is specified in the data/data_format.md file)
counts = read.csv(opt$countFile, header=T, row.names=1)

all = data.frame()
referenceGenotype = "wt"
ages = unique(sampleKey$Age)
dgeObjects = rep(NA, length(ages))
weightObjects  = rep(NA, length(ages))

for (i in 1:length(ages)){
  currentAge = ages[i]

  designMatrix = get_design_matrix(sampleKey, currentAge, referenceGenotype)
  dge = setup_dge(counts, currentAge)

  # Need to find a better way to plot mean vs variance
  dgeObjects[i] = list(voom(dge, designMatrix, plot=FALSE, save.plot = TRUE))
  
  weights = voomWithQualityWeights(dge, designMatrix,
                                   normalization = "none", plot=FALSE)
  weightObjects[i] = list(weights)
  
  vfit = lmFit(weights) %>% eBayes(.)
  hits = topTable(vfit, adjust="BH", number = Inf, sort.by = "F")
  hits$gene = rownames(hits)
  hits$timepoint = rep(currentAge,nrow(hits))
  
  all = rbind(all, hits)
}

# Write the gene expression data for each experiment
write.csv(all,file = "../data/all_diffexp_genes.csv",quote=T)

# Draw the mds plots
# pdf("../plots/mds_plots.pdf")
pdf(paste0(opt$plots, "/mds_plots.pdf"))
par(mfrow=c(2,3))
for (i in 1:length(ages)){
  get_mds_plot(dgeObjects[[i]], ages[i], sampleKey)
}
dev.off()

# Draw the mean-variance plots
pdf(paste0(opt$plots, "/mean_var_plots.pdf"))
par(mfrow=c(2,3))
for (i in 1:length(ages)){
  plot_mean_variance(dgeObjects[[i]], ages[i])
}
dev.off()

# Draw the sample weight plots
pdf(paste0(opt$plots, "/weight_plots.pdf"))
ggplotObjects = rep(NA, length(ages))
for (i in 1:length(ages)){
  ggplotObjects[i] = list(plot_weights(weightObjects[[i]], ages[i]))
}
source("./multiplot.R")
multiplot(ggplotObjects[[1]],
          ggplotObjects[[2]],
          ggplotObjects[[3]],
          ggplotObjects[[4]],
          ggplotObjects[[5]],
          cols = 3)
dev.off()
################### END ###################