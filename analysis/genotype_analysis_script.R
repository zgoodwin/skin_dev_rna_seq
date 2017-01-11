###############################################################################
# Program: genotype_analysis.R
# Project: Analysis of 
# Description: 
#
# Author: Zane Goodwin
# Date: 01/11/2017
#
###############################################################################

########################
#       PACKAGES       #
########################

library(limma)
library(edgeR)
library(stringr)
library(reshape)
library(ggplot2)

#######################
#      FUNCTIONS      #
#######################

####
# Name: setup_dge
# Purpose: 
# Args:
#     
setup_dge = function(counts_table, search_string){
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
  sample = sample_list[grep(search_string,sample_list$Age),]
  plot = plotMDS(dge_object, labels=paste(sample$Genotype, sample$Group),
                 main=paste(search_string, "MDS", sep=" "),
                 xlim=c(-1.5,1.5),
                 ylim=c(-1.5,1.5))  
  return(plot)
}

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

#######################
# EXECUTED STATEMENTS #
#######################

## File locations
COUNT_FILE = "../data/all_counts.csv"
SAMPLE_KEY_FILE = "../data/sample_key_obs.csv"
DESIGN_MATRICES = "../design_matrices/"
RESULTS_DIR = "../results/"

## Load the sample key, which contains metadata about each RNA-seq run 
#  (format is in the README file)
sampleKey = read.csv(SAMPLE_KEY_FILE, header=T)

## Load the read counts for each gene and for each experiment 
#  (format is in the README file)
counts = read.csv(COUNT_FILE, header=T, row.names=1)


## IMPORTANT: set up design matrices for each experiment. This is needed to fit
#  linear models between the wild type and knockout mice
e14_5_dm = get_design_matrix(sampleKey, "E14.5", "wt")
e15_5_dm = get_design_matrix(sampleKey, "E15.5", "wt")
e16_5_dm = get_design_matrix(sampleKey, "E16.5", "wt")
e17_5_dm = get_design_matrix(sampleKey, "E17.5", "wt")
nb_dm = get_design_matrix(sampleKey, "NB", "wt")


## Get DGE objects for each time point
e14_5_dge = setup_dge(counts, "E14.5")
e15_5_dge = setup_dge(counts, "E15.5")
e16_5_dge = setup_dge(counts, "E16.5")
e17_5_dge = setup_dge(counts, "E17.5")
nb_dge = setup_dge(counts, "NB")


## For each experiment, draw diagnostic MDS plots, write to plots folder
pdf("mds_plots.pdf")
par(mfrow=c(2,3))
p1 = get_mds_plot(e14_5_dge, "E14.5", sampleKey)
p2 = get_mds_plot(e15_5_dge, "E15.5", sampleKey)
p3 = get_mds_plot(e16_5_dge, "E16.5", sampleKey)
p4 = get_mds_plot(e17_5_dge, "E17.5", sampleKey)
p5 = get_mds_plot(nb_dge, "NB", sampleKey)
dev.off()


## Normalize read counts by fitting read counts to a negative binomial
#  distribution.
e14_5_norm = voom(e14_5_dge, e14_5_dm, plot=FALSE, save.plot = TRUE)
e15_5_norm = voom(e15_5_dge, e15_5_dm, plot=FALSE, save.plot = TRUE)
e16_5_norm = voom(e16_5_dge, e16_5_dm, plot=FALSE, save.plot = TRUE)
e17_5_norm = voom(e17_5_dge, e17_5_dm, plot=FALSE, save.plot = TRUE)
nb_norm = voom(nb_dge, nb_dm, plot=FALSE, save.plot = TRUE)

# Visualize the distribution of read counts to check the mean-variance trend. 
pdf("mean_variants_plots.pdf")
par(mfrow=c(2,3))
plot_mean_variance(e14_5_norm, "E14.5")
plot_mean_variance(e15_5_norm, "E15.5")
plot_mean_variance(e16_5_norm, "E16.5")
plot_mean_variance(e17_5_norm, "E17.5")
plot_mean_variance(nb_norm, "NB")
dev.off()

#### SAMPLE WEIGHT NORMALIZATION ####

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

######## FIT BAYES MODELS ########

# Fit bayes models to the weighted samples
e14_5_vfit = lmFit(e14_5_wts) %>% eBayes(.)
e15_5_vfit = lmFit(e15_5_wts) %>% eBayes(.)
e16_5_vfit = lmFit(e16_5_wts) %>% eBayes(.)
e17_5_vfit = lmFit(e17_5_wts) %>% eBayes(.)
nb_vfit = lmFit(nb_wts) %>% eBayes(.)

#### Obtain list of the top differentially expressed genes ####

diffexp_14_5 = topTable(e14_5_vfit, adjust="BH", number = Inf, sort.by = "F")
diffexp_15_5 = topTable(e15_5_vfit, adjust="BH", number = Inf, sort.by = "F")
diffexp_16_5 = topTable(e16_5_vfit, adjust="BH", number = Inf, sort.by = "F")
diffexp_17_5 = topTable(e17_5_vfit, adjust="BH", number = Inf, sort.by = "F")
diffexp_nb = topTable(nb_vfit, adjust="BH", number = Inf, sort.by = "F")

# Shift gene names to their own column for each dataset
diffexp_14_5$gene = rownames(diffexp_14_5)
diffexp_15_5$gene = rownames(diffexp_15_5)
diffexp_16_5$gene = rownames(diffexp_16_5)
diffexp_17_5$gene = rownames(diffexp_17_5)
diffexp_nb$gene = rownames(diffexp_nb)

# Add a a column for the timepoint.
diffexp_14_5$timepoint = rep(14.5,nrow(diffexp_14_5))
diffexp_15_5$timepoint = rep(15.5,nrow(diffexp_15_5))
diffexp_16_5$timepoint = rep(16.5,nrow(diffexp_16_5))
diffexp_17_5$timepoint = rep(17.5,nrow(diffexp_17_5))
diffexp_nb$timepoint   = rep("NB",nrow(diffexp_nb))

# Combine all data into a single table, long format.
all_diffexp_genes = rbind(diffexp_14_5,
                          diffexp_15_5,
                          diffexp_16_5,
                          diffexp_17_5,
                          diffexp_nb)

# Write and format gene expression tables.
write.csv(all_diffexp_genes,
          file = "../data/all_diffexp_genes.csv",
          quote=T)

################### END ###################