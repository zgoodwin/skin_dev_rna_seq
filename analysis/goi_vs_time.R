###############################################################################
# Program: goi_vs_time.R
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

library(reshape)
library(ggplot2)

######## ANALYZE GENES OF INTEREST ########

# Pull out the expression levels for the genes of interest for each time point.
goi = data.frame()
genes = c("g06152", "g15049", "g15227")
goi = rbind(goi,
            diffexp_14_5[which(diffexp_14_5$gene %in% genes),],
            diffexp_15_5[which(diffexp_15_5$gene %in% genes),],
            diffexp_16_5[which(diffexp_16_5$gene %in% genes),],
            diffexp_17_5[which(diffexp_17_5$gene %in% genes),],
            diffexp_nb[which(diffexp_nb$gene %in% genes),])

# Select columns
goi_select = goi[,c(1,2,7,8)]

# Melt the data frame
mdf = melt(goi_select, id.vars=3:4, measure.vars=1:2, variable_name = "genotype")
mdf$genotype = gsub("genotype", "", mdf$genotype)
mdf$genotype = gsub("het", "Het vs. WT", mdf$genotype)
mdf$genotype = gsub("mut", "Mut vs. WT", mdf$genotype)

# Plot the fold change in expression versus time
ggplot(data=mdf, aes(x=timepoint, y = value, colour=genotype)) + 
  geom_point() + 
  geom_line(aes(group=genotype)) + 
  facet_grid(~ gene) + 
  xlab("Timepoint") + 
  ylab("Log Fold-Change")
