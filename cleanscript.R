### libraries ###
library(ggplot2)
library(DESeq2)
library(tximport)
library(tidyverse)
library(ggsignif)
library(ggrepel)
library(clusterProfiler)
library(VennDiagram)
library(enrichplot)
library(stringr)
library(forcats)
library(org.Hs.eg.db)
library(ggpubr)
library(egg)


## Yi's tx2gene ## 
# tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
# setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

### DESeq2 ###

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds") # Yi's tx2gene
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

### all MALE samples ###
###### AD vs Old ######

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf',
          './Old/13-11T-Old.fastq.gz_counts/quant.sf',
          './Old/14-7A-Old.fastq.gz_counts/quant.sf',
          './Old/15-13T-Old.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11",
             "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

# Create the PCA plot with custom colors

PCA_ADOld <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue")) +
  scale_colour_manual(values = c("darkslateblue", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADOld  

# Darisia's edit plots function
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADOld <- edit_plots(PCA_ADOld)
PCA_ADOld 

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADOld, "./PCA_ADOld_all_11sept.rds")

###### Old vs Young ######

files = c('./Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf',
          './Old/13-11T-Old.fastq.gz_counts/quant.sf',
          './Old/14-7A-Old.fastq.gz_counts/quant.sf',
          './Old/15-13T-Old.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_OldYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "darkslateblue")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_OldYoung 

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_OldYoung <- edit_plots(PCA_OldYoung)
PCA_OldYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_all_11sept.rds")

###### AD vs Young ######

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_ADYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADYoung

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADYoung <- edit_plots(PCA_ADYoung)
PCA_ADYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_all_11sept.rds")

#### 3 samples per condition ####

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds") # Yi's tx2gene
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_ADOld <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue")) +
  scale_colour_manual(values = c("darkslateblue", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADOld 

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADOld <- edit_plots(PCA_ADOld)
PCA_ADOld

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADOld, "./PCA_ADOld_3each_11sept.rds")

###### Old vs Young ######
files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Young","Young","Young")))
names <-  c("Old4","Old5","Old6", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_OldYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "darkslateblue")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_OldYoung 


edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_OldYoung <- edit_plots(PCA_OldYoung)
PCA_OldYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_3each_11sept.rds")

###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)


sampleTable <- data.frame(condition = factor(c("AD","AD","AD","Young","Young","Young")))
names <-  c("AD3","AD6", "AD9", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_ADYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADYoung

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADYoung <- edit_plots(PCA_ADYoung)
PCA_ADYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_3each_11sept.rds")

## DESeq2 DE analysis results ##

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)
res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")
summary(res_ADvsOld)
# sum(res_ADvsOld$padj < 0.05, na.rm = TRUE)

## saving ##
setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_ADvsOld, file = "./DESEQres_pairwise_ADvsOld(3each-bulk).rds")

###### Old vs Young ######

files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Young","Young","Young")))
names <-  c("Old4","Old5","Old6", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)
res_OldvsYoung <- lfcShrink(dds_OldvsYoung, coef="condition_Old_vs_Young", type="apeglm")
summary(res_OldvsYoung)
# sum(res_OldvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##
setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_OldvsYoung, file = "./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")

###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD","AD","AD","Young","Young","Young")))
names <-  c("AD3","AD6", "AD9", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)
res_ADvsYoung <- lfcShrink(dds_ADvsYoung, coef="condition_AD_vs_Young", type="apeglm")
summary(res_ADvsYoung)
# sum(res_ADvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##
setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_ADvsYoung, file = "./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")


###### Looking at DESeq2 results ######
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

## AD vs Old ## 

res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

base_mean_values = res_ADvsOld@listData[["baseMean"]]
hist(log10(base_mean_values), ## checking base mean distribution for base mean cutoff
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution ADvsOld",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- filtered_ADvsOld_basemean[filtered_ADvsOld_basemean$padj < 0.05, ]
filtered_ADvsOld_LFC1 <- filtered_ADvsOld_basemean[filtered_ADvsOld_basemean$LFC > 1 | filtered_ADvsOld_basemean$LFC < -1, ]
filtered_ADvsOld_LFC1 <- na.omit(filtered_ADvsOld_LFC1)
saveRDS(filtered_ADvsOld_LFC1, "./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")

## Old vs Young ##

res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])

base_mean_values = res_OldvsYoung@listData[["baseMean"]]
hist(log10(base_mean_values), 
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution OldvsYoung",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")

filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- filtered_OldvsYoung_basemean[filtered_OldvsYoung_basemean$padj < 0.05, ]
filtered_OldvsYoung_LFC1 <- filtered_OldvsYoung_basemean[filtered_OldvsYoung_basemean$LFC > 1 | filtered_OldvsYoung_basemean$LFC < -1, ]
filtered_OldvsYoung_LFC1 <- na.omit(filtered_OldvsYoung_LFC1)
saveRDS(filtered_OldvsYoung_LFC1, "./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs YOUNG

res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])

base_mean_values = res_ADvsYoung@listData[["baseMean"]]
hist(log10(base_mean_values), 
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution ADvsYoung",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")

filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- filtered_ADvsYoung_basemean[filtered_ADvsYoung_basemean$padj < 0.05, ]
filtered_ADvsYoung_LFC1 <- filtered_ADvsYoung_basemean[filtered_ADvsYoung_basemean$LFC > 1 | filtered_ADvsYoung_basemean$LFC < -1, ]
filtered_ADvsYoung_LFC1 <- na.omit(filtered_ADvsYoung_LFC1)
saveRDS(filtered_ADvsYoung_LFC1, "./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## Venn diagrams ##
### filtered DESeq2 genes ###

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]

## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]

## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]


ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene

####  Venn UPREGULATED ####

gene_lists_up <- list(
  ADvsOld = ADOld_genes_up,
  OldvsYoung = OldYoung_genes_up,
  ADvsYoung = ADYoung_genes_up)

venn_result_up <- venn.diagram(
  x = gene_lists_up,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  # Place labels at default position (inside the sets)
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5, # transparency of colour circles
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold")

grid.newpage()
grid.draw(venn_result_up)

saveRDS(venn_result_up, "D://ABBY.windows_surface_pro/diffexp_results/venn_up_11sept.rds")

#### Venn DOWNREGULATED ####

gene_lists_down <- list(
  ADvsOld = ADOld_genes_down,
  OldvsYoung = OldYoung_genes_down,
  ADvsYoung = ADYoung_genes_down)

venn_result_down <- venn.diagram(
  x = gene_lists_down,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5,
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold") # Make category names bold

grid.newpage()
grid.draw(venn_result_down)

saveRDS(venn_result_down, "D://ABBY.windows_surface_pro/diffexp_results/venn_down_11sept.rds")


###### Looking at/plotting (volcano) DESeq2 results ######
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

### AD vs Old ### 

res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- na.omit(filtered_ADvsOld_basemean)
filtered_ADvsOld_basemean <- mutate(filtered_ADvsOld_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_ADvsOld_basemean, "./ADvsOld_dataframe_basemeancutoff_direction.rds")

## volcano plot no LFC cutoff ##
volcano <- ggplot(filtered_ADvsOld_basemean, aes(x = LFC, y = -log10(padj), colour = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "mediumpurple", "DOWN" = "mediumpurple", "notDE" = "grey82")) + # Color for non-significant and significant genes
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = " ") +
  theme(legend.position = "none")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

volcano <- edit_plots(volcano) 
volcano

## to get numbers of up- and down-regulated genes for labels on plot ##
num_upregulated <- sum(filtered_ADvsOld_basemean$Direction == "UP")
num_downregulated <- sum(filtered_ADvsOld_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano <- volcano +
  geom_text(aes(label = "650", #num_upregulated
                x = 7,
                y = 10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = "601", #num_downregulated
                x = -6,
                y=10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)

volcano
saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_ADOld.rds")

### Old vs Young ## #

res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- na.omit(filtered_OldvsYoung_basemean)
filtered_OldvsYoung_basemean <- mutate(filtered_OldvsYoung_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_OldvsYoung_basemean, "./OldvsYoung_dataframe_basemeancutoff_direction.rds")
## volcano plot no LFC cutoff
volcano <- ggplot(filtered_OldvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "lightgreen", "DOWN" = "lightgreen", "notDE" = "grey82")) + # Color for non-significant and significant genes
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = ""
  ) +
  theme(legend.position = "none")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano) 
volcano

num_upregulated <- sum(filtered_OldvsYoung_basemean$Direction == "UP")
num_downregulated <- sum(filtered_OldvsYoung_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano <- volcano +
  geom_text(aes(label = "1023",
                x = 7,
                y = 10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = "1406",
                x = -6,
                y=10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)
volcano
saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_OldYoung.rds")

### AD vs Young ###
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- na.omit(filtered_ADvsYoung_basemean)
filtered_ADvsYoung_basemean <- mutate(filtered_ADvsYoung_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_ADvsYoung_basemean, "./ADvsYoung_dataframe_basemeancutoff_direction.rds")

## volcano plot no LFC cutoff
volcano <- ggplot(filtered_ADvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "orange", "DOWN" = "orange", "notDE" = "grey82")) + # Color for non-significant and significant genes
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "") +
  theme(legend.position = "none")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano) 
volcano

num_upregulated <- sum(filtered_ADvsYoung_basemean$Direction == "UP")
num_downregulated <- sum(filtered_ADvsYoung_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano <- volcano +
  geom_text(aes(label = "2243",
                x = 7,
                y = 12),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = "2716",
                x = -6,
                y=12),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)
volcano
saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_ADYoung.rds")



##### CLUSTERPROFILER - (CONFIDENCE that genes are AD-REALTED) #####
## GENES i put in = all filters applied ##
## BACKGROUND/UNIVERSE = only base mean filter (what IS expressed in cells) ##

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

### UNIVERSES ###
## AD vs OLD
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])
filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]

## OLD vs YOUNG
# res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
# OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
#                             LFC = res_OldvsYoung@listData[["log2FoldChange"]],
#                             pvalue = res_OldvsYoung@listData[["pvalue"]],
#                             padj = res_OldvsYoung@listData[["padj"]],
#                             basemean = res_OldvsYoung@listData[["baseMean"]])
# filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
# 
# ## AD vs YOUNG
# res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
# ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
#                            LFC = res_ADvsYoung@listData[["log2FoldChange"]],
#                            pvalue = res_ADvsYoung@listData[["pvalue"]],
#                            padj = res_ADvsYoung@listData[["padj"]],
#                            basemean = res_ADvsYoung@listData[["baseMean"]])
# filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]

## INPUT GENES ##
### filtered DESeq2 genes ###
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

### REMOVING VERSION NUMBERS OF GENES
ADvsOld$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsOld$Gene)
filtered_ADvsOld_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsOld_basemean$Gene)
ADvsOld_genes <- ADvsOld$Gene

OldvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",OldvsYoung$Gene)
#filtered_OldvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_OldvsYoung_basemean$Gene)
OldvsYoung_genes <- OldvsYoung$Gene

ADvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsYoung$Gene)
#filtered_ADvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsYoung_basemean$Gene)
ADvsYoung_genes <- ADvsYoung$Gene

### to get BPs for AD Old up/down separately
## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]

## genes
ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene

## AD-related genes = (common in ADOld & ADYoung) - (common in OldYoung & ADYoung)
unique_ADOld_up <- intersect(ADOld_genes_up, ADYoung_genes_up)
unique_ADOld_up <- setdiff(unique_ADOld_up, OldYoung_genes_up)

unique_ADOld_down <- intersect(ADOld_genes_down, ADYoung_genes_down)
unique_ADOld_down <- setdiff(unique_ADOld_down, OldYoung_genes_down)

AD_genes <- c(unique_ADOld_up, unique_ADOld_down)
AD_df <- ADvsOld[ADvsOld$Gene %in% AD_genes, ]

### AD vs Old ###

## upregulated BPs
## no background
GOenrich_ADOld <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld2 <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld, "./GOenrich_diseasegenes_up_nobackground.rds")
saveRDS(GOenrich_ADOld2, "./GOenrich_diseasegenes_up_+background.rds")

## downregulated BPs
## no background
GOenrich_ADOld3 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld4 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld3, "./GOenrich_diseasegenes_down_nobackground.rds")
saveRDS(GOenrich_ADOld4, "./GOenrich_diseasegenes_down_+background.rds")

### PLOTTING CLUSTERPROFILER ###
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

enrich_ADOld_up <- readRDS("./GOenrich_diseasegenes_up_+background.rds")
enrich_ADOld_down <- readRDS("./GOenrich_diseasegenes_down_+background.rds")
res_enrich_ADOld_up <- as.data.frame(enrich_ADOld_up@result)
res_enrich_ADOld_down <- as.data.frame(enrich_ADOld_down@result)

## filtering top 20 BPs
### ordering by decreasing count ###

column_to_order <- "Count" # Column to order by

# Order the data frame based on the specified column
ordered_res_enrich_ADOld_up <- res_enrich_ADOld_up[order(res_enrich_ADOld_up[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADOld_up <- ordered_res_enrich_ADOld_up[1:20, ]
df_enrich_ADOld_up <- data.frame(BP = ordered_res_enrich_ADOld_up$Description, padj = ordered_res_enrich_ADOld_up$p.adjust, counts = ordered_res_enrich_ADOld_up$Count)

ordered_res_enrich_ADOld_down <- res_enrich_ADOld_down[order(res_enrich_ADOld_down[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADOld_down <- ordered_res_enrich_ADOld_down[1:20, ]
df_enrich_ADOld_down <- data.frame(BP = ordered_res_enrich_ADOld_down$Description, padj = ordered_res_enrich_ADOld_down$p.adjust, counts = ordered_res_enrich_ADOld_down$Count)

disease_up <- readRDS("./GOenrich_diseasegenes_up_+background.rds")
disease_down <- readRDS("./GOenrich_diseasegenes_down_+background.rds")

res_disease_up <- as.data.frame(disease_up@result)
res_disease_down <- as.data.frame(disease_down@result)

# Column to order by
column_to_order <- "Count"
# Order the data frame based on the specified column
ordered_disease_up <- res_disease_up[order(res_disease_up[[column_to_order]], decreasing = TRUE), ]
ordered_disease_down <- res_disease_down[order(res_disease_down[[column_to_order]], decreasing = TRUE), ]

saveRDS(ordered_disease_down, "./ordered_diseasegenes_down.rds")
saveRDS(ordered_disease_up, "./ordered_diseasegenes_up.rds")

up <- readRDS("./ordered_diseasegenes_up.rds")
down <- readRDS("./ordered_diseasegenes_down.rds")

up_sig <- up[up$p.adjust<0.05, ]
down_sig <- down[down$p.adjust <0.05, ]

column_to_order <- "Count"

# Order the data frame based on the specified column
ordered_up <- up_sig[order(up_sig[[column_to_order]], decreasing = TRUE), ]
df_up <- data.frame(BP = ordered_up$Description, padj = ordered_up$p.adjust, counts = ordered_up$Count)

ordered_down <- down_sig[order(down_sig[[column_to_order]], decreasing = TRUE), ]
ordered_down <- ordered_down[1:20, ]
df_down <- data.frame(BP = ordered_down$Description, padj = ordered_down$p.adjust, counts = ordered_down$Count)

## up
plot_enrich_up <- ggplot(df_up, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "springgreen") +  
  labs(title = "AD-specific Significant Enriched Terms (UP)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0)) 


#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

plot_enrich_up<- edit_plots(plot_enrich_up)
plot_enrich_up

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(plot_enrich_up, "./diseasegenes_sigBPs_up.rds")


## down
plot_enrich_down <- ggplot(df_down, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "red") +  
  labs(title = "AD-specific Significant Enriched Terms (DOWN)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0))

plot_enrich_down <- edit_plots(plot_enrich_down)
plot_enrich_down

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(plot_enrich_up, "./diseasegenes_sigBPs_down.rds")


## up (all processes - including not sig) ## 
plot_enrich_ADOld_up <- ggplot(df_enrich_ADOld_up, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "springgreen") +  
  labs(title = "AD-specific Significant Enriched Terms (UP)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0)) 

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

plot_enrich_ADOld_up<- edit_plots(plot_enrich_ADOld_up)
plot_enrich_ADOld_up

saveRDS(plot_enrich_ADOld_up, "D://ABBY.windows_surface_pro/diffexp_results/clusterprof_disease_up.rds")
saveRDS(plot_enrich_ADOld_up, "./clusterprof_UP_09-09.rds")

## down
plot_enrich_ADOld_down <- ggplot(df_enrich_ADOld_down, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "red") +  
  labs(title = "AD-specific Significant Enriched Terms (DOWN)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0))

plot_enrich_ADOld_down <- edit_plots(plot_enrich_ADOld_down)
plot_enrich_ADOld_down

saveRDS(plot_enrich_ADOld_down, "D://ABBY.windows_surface_pro/diffexp_results/clusterprof_disease_down.rds")
saveRDS(plot_enrich_ADOld_down, "./clusterprof_DOWN_09-09.rds")




####### QC PLOTS #######
### fixing multiqc ###
### downloaded data from multiQc html report ###

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

## per sequence quality score (mean) ##
seq_quality <- read_csv("./fastqc_per_sequence_quality_scores_plot.csv")
seq_quality <- as.data.frame(seq_quality)
seq_quality_long <- seq_quality %>%
  pivot_longer(cols = -`Mean Sequence Quality (Phred Score)`, 
               names_to = "Sample", 
               values_to = "Count")

seq_quality_long <- mutate(seq_quality_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

seq_quality_plot <- ggplot(seq_quality_long, aes(x = `Mean Sequence Quality (Phred Score)`, y = Count, colour = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Mean Sequence Quality (Phred Score)", y = "Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_x_continuous(limits = c(0, 36)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

seq_quality_plot <- edit_plots(seq_quality_plot) 
seq_quality_plot 
saveRDS(seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/seqqual_plot_all_11sept.rds")

## per seq qual score for selected samples ## 

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", 
"13-11T-Old", "14-7A-Old", "15-13T-Old", 
"2-12A-Young", "4-13A-Young", "6-14A-Young")

# Filter the dataframe to include only rows with a specific value in a specific column
selected_seq_quality_long <- seq_quality_long %>%
  filter(Sample %in% samples_keep)

selected_seq_quality_plot <- ggplot(selected_seq_quality_long, aes(x = `Mean Sequence Quality (Phred Score)`, y = Count, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Mean Sequence Quality (Phred Score)", y = "Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_x_continuous(limits = c(0, 36)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

selected_seq_quality_plot <- edit_plots(selected_seq_quality_plot)
selected_seq_quality_plot

saveRDS(selected_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/seqqual_plot_3each_11sept.rds")

### per base quality score ###

pb_seq_quality <- read_csv("./fastqc_per_base_sequence_quality_plot.csv")
pb_seq_quality <- as.data.frame(pb_seq_quality)
pb_seq_quality_long <- pb_seq_quality %>%
  pivot_longer(cols = -`Position (bp)`, 
               names_to = "Sample", 
               values_to = "Phred Score")

pb_seq_quality_long <- mutate(pb_seq_quality_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

pb_seq_quality_plot <- ggplot(pb_seq_quality_long, aes(x = `Position (bp)`, y = `Phred Score`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Phred Score", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_seq_quality_plot <- edit_plots(pb_seq_quality_plot)
pb_seq_quality_plot 

saveRDS(pb_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_seqqual_plot_all_11sept.rds")

## per base quality score for selected samples 

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", 
"13-11T-Old", "14-7A-Old", "15-13T-Old", 
"2-12A-Young", "4-13A-Young", "6-14A-Young")

# Filter the dataframe to include only rows with a specific value in a specific column
pb_selected_seq_quality_long <- pb_seq_quality_long %>%
  filter(Sample %in% samples_keep)

pb_selected_seq_quality_plot <- ggplot(pb_selected_seq_quality_long, aes(x = `Position (bp)`, y = `Phred Score`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Phred Score", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_selected_seq_quality_plot <- edit_plots(pb_selected_seq_quality_plot)
pb_selected_seq_quality_plot

saveRDS(pb_selected_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_seqqual_plot_3each_11sept.rds")

### per base N content ##

pb_ncont <- read_csv("./fastqc_per_base_n_content_plot.csv")
pb_ncont <- as.data.frame(pb_ncont)
pb_ncont_long <- pb_ncont %>%
  pivot_longer(cols = -`Position in Read (bp)`, 
               names_to = "Sample", 
               values_to = "Percentage N-Count")

pb_ncont_long <- mutate(pb_ncont_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

pb_ncont_plot <- ggplot(pb_ncont_long, aes(x = `Position in Read (bp)`, y = `Percentage N-Count`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Percentage N-Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 6)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_ncont_plot <- edit_plots(pb_ncont_plot)
pb_ncont_plot

saveRDS(pb_ncont_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_Ncontent_plot_all_11sept.rds")

## per base N content for selected samples ##

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", 
"13-11T-Old", "14-7A-Old", "15-13T-Old", 
"2-12A-Young", "4-13A-Young", "6-14A-Young")

# Filter the dataframe to include only rows with a specific value in a specific column
selected_pb_ncont_long <- pb_ncont_long %>%
  filter(Sample %in% samples_keep)

selected_pb_ncont_long_plot <- ggplot(selected_pb_ncont_long, aes(x = `Position in Read (bp)`, y = `Percentage N-Count`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Percentage N-Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 6)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

selected_pb_ncont_long_plot <- edit_plots(selected_pb_ncont_long_plot)
selected_pb_ncont_long_plot

saveRDS(selected_pb_ncont_long_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_Ncontent_plot_3each_9sept.rds")

###### PLOT CUMULATIVE OTAL GENE COUNTS PER SAMPLE (proxy for sequencing depth) ######
## all samples ## 

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/27-5A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', 
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', 
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', 
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', 
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', 
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
            "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
            "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

# Create a bar plot using ggplot2
plot_count <- ggplot(total_counts_df, aes(x = SampleLabel, y = TotalCount, fill = condition, colour = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Cumulative Gene Count", fill = "", colour = "") +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) # Rotate x-axis labels 

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

plot_count <- edit_plots(plot_count)
plot_count

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_bar_all_11sept.png", plot_count, width = 15, height = 10)
saveRDS(plot_count, "./cumulativegenecount_bar_all_11sept.rds")

### box plot total gene counts all samples ###

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 27900000, fontface="bold")
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

boxplot_counts

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_box_all_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts, "./cumulativegenecount_box_all.rds")

###### PLOT TOTAL GENE COUNTS PER SAMPLE (selected samples) ######
#### boxlot ####

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))
# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

boxplot_counts

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_box_3each_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts, "./cumulativegenecount_box_3each_11.rds")


### plot salmon mapping rate for all samples ###

sample_names <- c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
                  "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
                  "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")

map_perc <- c(45.8, 44.9, 52.8, 52.4, 49.4, 51.8, 53, 48, 49.3, 59.1, 49.1, 44.4,
              55.7, 50.1, 54.8, 54.1, 49.4, 49.7, 53.3, 57.3, 56.8, 43.5,
              68.3, 44.4, 65.5, 56.3, 53.3, 52.5, 59.3, 46.2)

condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

# Bar plot
plot_count <- ggplot(data, aes(x = SampleLabel, y = MappedPercentage, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Mapped Percentage (%)", fill = "", colour = "") +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) # Rotate x-axis labels 

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
plot_count <- edit_plots(plot_count)
plot_count

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_bar_all_11sept.png", plot_count, width = 15, height = 10)
saveRDS(plot_count,"./maprate_bar_all_11sept.rds")

## boxplot ## 

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, ## 'ns' determined by kruskal-wallis test ##
              textsize = 4, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

boxplot_counts
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_box_all_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts,"./maprate_box_all_11sept.rds")

### plot salmon mapping rate selected samples ###

sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition<- c("AD", "AD", "AD", "Old", "Old", "Old", "Young", "Young", "Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

### boxplot mapping rate ### 

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

boxplot_counts
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_box_3each_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts,"./maprate_box_3each_11sept.rds")

#### boxplot for ages (selected samples) ####

samples <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
ages <- c(71, 61, 79, 68, 70, 77, 42, 57, 59)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
donorinfo <- data.frame(Sample = samples, Age = ages, Condition = condition)

boxplot_ages <- ggplot(donorinfo, aes(x = Condition, y = Age)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Age",
       colour = "") +
  ylim(30, 80)   # Set y-axis limits

boxplot_ages <- edit_plots(boxplot_ages)
boxplot_ages

boxplot_ages <- boxplot_ages +
  geom_point(aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= donorinfo$Sample, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

boxplot_ages
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./ages_box_3each_11sept.png", boxplot_ages, width = 15, height = 10)
saveRDS(boxplot_ages, "./ages_box_3each.rds")



### CNET PLOTS ###
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

enrich_AD_up <- readRDS("./GOenrich_diseasegenes_up_+background.rds") # cnet requires enrichGO result
enrich_AD_down <- readRDS("./GOenrich_diseasegenes_down_+background.rds")


cnet_up <- cnetplot(enrich_AD_up, showCategory=c("gliogenesis", "glial cell differentiation", "ensheathment of neurons", "axon ensheathment", "myelination"), 
                    circular = TRUE, colorEdge = TRUE, node_label = "gene", 
                    color.params = list(category = c("#00B0F6", "#00BF7D", "#A3A500", "#F8766D", "#E76BF3"), gene = "grey60")) +
  theme(legend.position = "none")

cnet_down <- cnetplot(enrich_AD_down, showCategory=c("modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling", "regulation of monoatomic ion transport", "regulation of monoatomic ion transmembrane transport", "calcium ion transport"), 
                      circular=TRUE, colorEdge = TRUE, node_label = "gene", 
                      color.params = list(category = c("#A3A500", "#E76BF3", "#00B0F6", "#00BF7D",  "#F8766D"), gene = "grey60")) +
  theme(legend.position = "none")

setwd("D://ABBY.windows_surface_pro/diffexp_results/")

saveRDS(cnet_up, "./cnet_up_final.rds")
saveRDS(cnet_down, "./cnet_down_final.rds")





######### STATS ######### 
library(dplyr)     # for data manipulation (optional)

# Create a data frame with your data 
tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

##### CUMULATIVE GENE COUNTS PER SAMPLE (all samples) ######

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/27-5A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', 
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', 
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', 
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', 
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', 
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
            "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
            "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(TotalCount = totals, condition)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(TotalCount ~ condition, data = total_counts_df)
kruskal_result


###### CUMULATIVE GENE COUNTS PER SAMPLE (selected samples) ######

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(TotalCount = totals, condition)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(TotalCount ~ condition, data = total_counts_df)
kruskal_result

### salmon mapping rate for all samples

map_perc <- c(45.8, 44.9, 52.8, 52.4, 49.4, 51.8, 53, 48, 49.3, 59.1, 49.1, 44.4,
              55.7, 50.1, 54.8, 54.1, 49.4, 49.7, 53.3, 57.3, 56.8, 43.5,
              68.3, 44.4, 65.5, 56.3, 53.3, 52.5, 59.3, 46.2)

condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Create a data frame
data <- data.frame(MappedPercentage = map_perc, condition)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(MappedPercentage ~ condition, data = data)
kruskal_result

### salmon mapping rate for selected samples

condition <- c("AD", "AD", "AD", "Old", "Old", "Old", "Young", "Young", "Young")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)

# Create a data frame
data <- data.frame(MappedPercentage = map_perc, condition)

kruskal_result <- kruskal.test(MappedPercentage ~ condition, data = data)
kruskal_result









#### GGARRANGE MULTI-PANEL FIGURES ####


setwd("D://ABBY.windows_surface_pro/diffexp_results/")

## MultiQC ALL

plot1 <- readRDS("./seqqual_plot_all_11sept.rds")
plot2 <- readRDS("./pb_seqqual_plot_all_11sept.rds")
plot3 <- readRDS("./pb_Ncontent_plot_all_11sept.rds")

seqqual <- plot1 +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
edit_plots <- function( ## change axis text size and axis label size
  ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.40, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

seqqual <- edit_plots(seqqual)
seqqual

# saved pdf 10 x 10
saveRDS(seqqual, "./seqqual_newest_all.rds")

pbseqqual <- plot2 +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
pbseqqual <- edit_plots(pbseqqual)
pbseqqual

# saved pdf 10 x 10
saveRDS(pbseqqual, "./pbseqqual_newest_all.rds")

ncont <- plot3 +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

ncont <- edit_plots(ncont)
ncont

# saved pdf 10 x 10
saveRDS(ncont, "./ncont_newest_all.rds")

multi_panel_figure_all <- ggarrange(seqqual, pbseqqual, ncont, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.17,
                                    align = "hv",
                                    vjust = 1.8)
multi_panel_figure_all


## MultiQC plots

plot1_all <- readRDS("./seqqual_biglegend_all.rds")
plot2_all <- readRDS("./pbseqqual_biglegend_all.rds")
plot3_all <- readRDS("./ncont_biglegend_all.rds")

multi_panel_figure_all <- ggarrange(plot1_all, plot2_all, plot3_all, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.17,
                                    align = "hv",
                                    vjust = 1.8)
multi_panel_figure_all

## MultiQC 3 each

plot1_all <- readRDS("./seqqual_plot_3each_11sept.rds")
plot2_all <- readRDS("./pb_seqqual_plot_3each_11sept.rds")
plot3_all <- readRDS("./pb_Ncontent_plot_3each_9sept.rds")

seqqual <- plot1_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
edit_plots <- function( ## change axis text size and axis label size
  ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.40, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

seqqual <- edit_plots(seqqual)
seqqual

# saved pdf 10 x 10
saveRDS(seqqual, "./seqqual_newest_3each.rds")

pbseqqual <- plot2_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
pbseqqual <- edit_plots(pbseqqual)
pbseqqual


# saved pdf 10 x 10
saveRDS(pbseqqual, "./pbseqqual_newest_3each.rds")

ncont <- plot3_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

ncont <- edit_plots(ncont)
ncont

# saved pdf 10 x 10
saveRDS(ncont, "./ncont_newest_3each.rds")

multi_panel_figure_3each <- ggarrange(seqqual, pbseqqual, ncont, ncol = 3,
                                      common.legend = TRUE,   # Add a single legend for all plots
                                      legend = "top", # Adjust the legend position as needed
                                      labels = "AUTO",
                                      label.x = 0.17,
                                      align = "hv",
                                      vjust = 1.8)
multi_panel_figure_3each



## MultiQC plots

plot1 <- readRDS("./seqqual_biglegend_3each.rds")
plot2 <- readRDS("./pbseqqual_biglegend_3each.rds")
plot3 <- readRDS("./ncont_biglegend_3each.rds")

multi_panel_figure_3each <- ggarrange(plot1, plot2, plot3, ncol = 3,
                                      common.legend = TRUE,   # Add a single legend for all plots
                                      legend = "top", # Adjust the legend position as needed
                                      labels = "AUTO",
                                      label.x = 0.17,
                                      align = "hv",
                                      vjust = 1.8)
multi_panel_figure_3each


## mapping rate

### plot salmon mapping rate for all samples

sample_names <- c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
                  "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
                  "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")

map_perc <- c(45.8, 44.9, 52.8, 52.4, 49.4, 51.8, 53, 48, 49.3, 59.1, 49.1, 44.4,
              55.7, 50.1, 54.8, 54.1, 49.4, 49.7, 53.3, 57.3, 56.8, 43.5,
              68.3, 44.4, 65.5, 56.3, 53.3, 52.5, 59.3, 46.2)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "", colour=" ") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./maprate_box_all_new.rds")


### plot salmon mapping rate selected samples ###

sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition<- c("AD", "AD", "AD", "Old", "Old", "Old", "Young", "Young", "Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

### boxplot mapping rate ### 

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits


boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./maprate_box_3each_new.rds")




#### cumulative gene counts
## all samples

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/27-5A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', 
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', 
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', 
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', 
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', 
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
            "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
            "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./cumulativegenecount_box_all_new.rds")


## selected samples

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))
# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()
boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(boxplot_counts, "./cumulativegenecount_box_3each_new.rds")





#### PCA
### all 

setwd("D://ABBY.windows_surface_pro/diffexp_results/")

plot1_all <- readRDS("./PCA_ADOld_all_11sept.rds")
plot2_all <- readRDS("./PCA_OldYoung_all_11sept.rds")
plot3_all <- readRDS("./PCA_ADYoung_all_11sept.rds")

ADOld <- plot1_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADOld
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADOld <- edit_plots(ADOld)
ADOld
# saved pdf 9 x 9
saveRDS(ADOld, "./ADOld_PCA_biglegend_all.rds")

OldYoung <- plot2_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
OldYoung

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
OldYoung <- edit_plots(OldYoung)
OldYoung

# saved pdf 9 x 9
saveRDS(OldYoung, "./OldYoung_PCA_biglegend_all.rds")

ADYoung <- plot3_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADYoung
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADYoung <- edit_plots(ADYoung)
ADYoung

# saved pdf 9 x 9
saveRDS(ADYoung, "./ADYoung_PCA_biglegend_all.rds")


## MultiQC plots

plot1 <- readRDS("./ADOld_PCA_biglegend_all.rds")
plot2 <- readRDS("./OldYoung_PCA_biglegend_all.rds")
plot3 <- readRDS("./ADYoung_PCA_biglegend_all.rds")

multi_panel_figure_all <- ggarrange(plot1, plot2, plot3, ncol = 3,
                                    common.legend = FALSE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.17,
                                    label.y = 0.7,
                                    align = "hv",
                                    vjust = 1.8)

multi_panel_figure_all
saveRDS(multi_panel_figure_all, "./multi_PCA_all.rds")

### 3 each

setwd("D://ABBY.windows_surface_pro/diffexp_results/")

plot1_all <- readRDS("./PCA_ADOld_3each_11sept.rds")
plot2_all <- readRDS("./PCA_OldYoung_3each_11sept.rds")
plot3_all <- readRDS("./PCA_ADYoung_3each_11sept.rds")

ADOld <- plot1_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADOld
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADOld <- edit_plots(ADOld)
ADOld
# saved pdf 9 x 9
saveRDS(ADOld, "./ADOld_PCA_biglegend_3each.rds")

OldYoung <- plot2_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
OldYoung

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
OldYoung <- edit_plots(OldYoung)
OldYoung

# saved pdf 9 x 9
saveRDS(OldYoung, "./OldYoung_PCA_biglegend_3each.rds")

ADYoung <- plot3_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADYoung
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADYoung <- edit_plots(ADYoung)
ADYoung

# saved pdf 9 x 9
saveRDS(ADYoung, "./ADYoung_PCA_biglegend_3each.rds")


## MultiQC plots

plot1 <- readRDS("./ADOld_PCA_biglegend_3each.rds")
plot2 <- readRDS("./OldYoung_PCA_biglegend_3each.rds")
plot3 <- readRDS("./ADYoung_PCA_biglegend_3each.rds")

multi_panel_figure_all <- ggarrange(plot1, plot2, plot3, ncol = 3,
                                    common.legend = FALSE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = c("D", "E", "F"),
                                    label.x = 0.17,
                                    label.y = 0.75,
                                    align = "hv",
                                    vjust = 1.8)

multi_panel_figure_all
saveRDS(multi_panel_figure_all, "./multi_PCA_3each_DEF.rds")


## ages

#### boxplot for ages ####

samples <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
ages <- c(71, 61, 79, 68, 70, 77, 42, 57, 59)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
donorinfo <- data.frame(Sample = samples, Age = ages, Condition = condition)
boxplot_ages <- ggplot(donorinfo, aes(x = Condition, y = Age)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Age",
       colour = "") +
  ylim(30, 80)   # Set y-axis limits

boxplot_ages <- boxplot_ages +
  geom_point(aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= donorinfo$Sample, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_ages
boxplot_ages <- boxplot_ages +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
boxplot_ages

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
boxplot_ages <- edit_plots(boxplot_ages)
boxplot_ages

saveRDS(boxplot_ages, "./boxplot_ages_3each_newest.rds")



### venn diagrams

### filtered DESeq2 genes ###
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")
ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]

## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]

## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]


ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene


## Venn diagrams ##
#### UPREGULATED ####

gene_lists_up <- list(
  ADvsOld = ADOld_genes_up,
  OldvsYoung = OldYoung_genes_up,
  ADvsYoung = ADYoung_genes_up)

venn_up <- venn.diagram(
  x = gene_lists_up,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  # Place labels at default position (inside the sets)
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5,
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold"
)

grid.newpage()
grid.draw(venn_up)

saveRDS(venn_up, "D://ABBY.windows_surface_pro/diffexp_results/venn_up_new.rds")

#### DOWNREGULATED ####
gene_lists_down <- list(
  ADvsOld = ADOld_genes_down,
  OldvsYoung = OldYoung_genes_down,
  ADvsYoung = ADYoung_genes_down)

venn_down <- venn.diagram(
  x = gene_lists_down,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5,
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold" # Make category names bold
  #fontfamily = "sans"  # Font family (adjust as needed)
)

grid.newpage()
grid.draw(venn_down)

saveRDS(venn_down, "D://ABBY.windows_surface_pro/diffexp_results/venn_down_new.rds")

setwd("D://ABBY.windows_surface_pro/diffexp_results")
venn_up <- readRDS("./venn_up_new.rds")
venn_down <- readRDS("./venn_down_new.rds")

multi_venn <- ggarrange(venn_up, venn_down, ncol = 2,
                        common.legend = FALSE,   # Add a single legend for all plots
                        legend = "top", # Adjust the legend position as needed
                        labels = "AUTO",
                        label.x = 0.17,
                        label.y = 0.75,
                        align = "hv",
                        vjust = 1.8)
grid.newpage()
grid.draw(multi_venn)
saveRDS(multi_venn, "./multi_venn.rds")

venn_down

grid.newpage()
grid.draw(venn_down)




### volcanos
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ADOld <- readRDS("./volcano_ADOld.rds")
OldYoung <- readRDS("./volcano_OldYoung.rds")
ADYoung <- readRDS("./volcano_ADYoung.rds")

ADOld <- ADOld +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADOld
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADOld <- edit_plots(ADOld)
ADOld

OldYoung <- OldYoung +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
OldYoung
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
OldYoung <- edit_plots(OldYoung)
OldYoung

ADYoung <- ADYoung +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ADYoung
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
ADYoung <- edit_plots(ADYoung)
ADYoung

multi_volcano <- ggarrange(ADOld, OldYoung, ADYoung, ncol = 3,
                           common.legend = TRUE,   # Add a single legend for all plots
                           legend = "none", # Adjust the legend position as needed
                           labels = "AUTO",
                           label.x = 0.17,
                           label.y = 1,
                           align = "hv",
                           vjust = 1.8)

multi_volcano


## clusterProfiler

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
UP <- readRDS("./clusterprof_disease_up.rds")
DOWN <- readRDS("./clusterprof_disease_down.rds")

multi <- ggarrange(UP, DOWN, ncol = 2,
                   common.legend = TRUE,   # Add a single legend for all plots
                   legend = "none", # Adjust the legend position as needed
                   labels = "AUTO",
                   label.x = 0,
                   label.y = 1,
                   align = "hv",
                   vjust = 1.8)

multi
