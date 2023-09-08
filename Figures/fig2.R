library(ggplot2)
library(plyr)

source('../R_functions/processing_prior_pm_clustering.R')


# Colorblind-friendly palette
cbfly <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

set.seed(1) #To ensure legends don't overlap data points

# Plot PCA TM Lung:
mydata <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_saver.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_labels.csv')
var_cutoff=2

rownames(mydata) = mydata[,1]
mydata <- mydata[,-1]
true_labels <- true_labels[,-1]

pcadata <- processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)
pcadata <- as.data.frame(t(pcadata))
pcadata$true_labels <- as.factor(true_labels)
pcadata$true_labels <- revalue(pcadata$true_labels,c('1'='lung', '2'='alveolar' ,
                                                   '3'= 'mast & im.',
                                                   '4'= 'circ. mono.' ,
                                                   '5'='inv. mono.',
                                                   '6'='multicil.',
                                                   '7'='dendritic'))

p1 <- ggplot(pcadata, aes(x = PC1, y = PC2)) + geom_point(aes( color = true_labels))
p1 <- p1 + theme(legend.position = c(0.88, 0.22), aspect.ratio = 1)
p1 <- p1 + scale_colour_manual(values=cbfly)
p1

pdf('../Figures/Fig_2_PCA_TM_Lung.pdf')
p1
dev.off()

# Plot UMAP TM Lung:
linnormdata_full <- Linnorm::Linnorm(mydata)
linnormdata <- processing_prior_pm_clustering(dataset=linnormdata_full,
                                              LogNormalization=FALSE,
                                              var_cutoff=var_cutoff,
                                              LocalAvgAllPts = FALSE) #PCA of linnorm

umap_embedding_full <- umap::umap(t(linnormdata), n_components = 2)
umap_embedding <- as.data.frame(umap_embedding_full$layout)
umap_embedding$true_labels <- pcadata$true_labels


p2<-ggplot(umap_embedding,aes(x=V1, y=V2)) + geom_point(aes( color=true_labels)) 
p2<-p2+labs(x ='UMAP1',y='UMAP2') + theme(legend.position = c(0.12, 0.22),
                                          aspect.ratio=1)
p2 <- p2 + scale_colour_manual(values=cbfly)
p2

pdf('../Figures/Fig_2_UMAP_TM_Lung.pdf')
p2
dev.off()

# Plot PCA TM Panc:

mydata <- read.csv(file = '../Data_after_Imputation/pancreas_tabula_muris_saver.csv')
true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_tabula_muris_labels.csv')
var_cutoff = 1

rownames(mydata) = mydata[,1]
mydata <- mydata[,-1]
true_labels <- true_labels[,-1]

pcadata <- processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,
                                          var_cutoff=var_cutoff)
pcadata <- as.data.frame(t(pcadata))
pcadata$true_labels <- as.factor(true_labels)
pcadata$true_labels <- revalue(pcadata$true_labels,c('1'='acinar','2'='beta',
                                                     '3'='stellate',
                                                     '4'='pancreatic A',
                                                     '5'='ductal',
                                                     '6'='pancreatic D',
                                                     '7'='pancreatic PP'))


p3 <- ggplot(pcadata,aes(x=PC1, y=PC2)) +geom_point(aes( color=true_labels))
p3 <- p3 + theme(legend.position = c(0.14, 0.22),aspect.ratio=1)
p3 <- p3 + scale_colour_manual(values=cbfly)
p3

pdf('../Figures/Fig_2_PCA_TM_Panc.pdf')
p3
dev.off()

# Plot UMAP TM Panc:

linnormdata_full <- Linnorm::Linnorm(mydata)
linnormdata <- processing_prior_pm_clustering(dataset=linnormdata_full, 
                                              LogNormalization=FALSE, 
                                              var_cutoff=var_cutoff,
                                              LocalAvgAllPts = FALSE) #PCA of linnorm

umap_embedding_full <- umap::umap(t(linnormdata),n_components=2)
umap_embedding <- as.data.frame(umap_embedding_full$layout)
umap_embedding$true_labels <- pcadata$true_labels
#axis.title.x='UMAP1',axis.title.y='UMAP2'

p4 <- ggplot(umap_embedding,aes(x=V1, y=V2)) + geom_point(aes( color=true_labels))
p4 <- p4 + labs(x = 'UMAP1', y = 'UMAP2') + theme(legend.position = c(0.14, 0.78), aspect.ratio = 1)
p4 <- p4 + scale_colour_manual(values = cbfly)
p4

pdf('../Figures/Fig_2_UMAP_TM_Panc.pdf')
p4
dev.off()