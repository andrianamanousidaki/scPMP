library(ggplot2)
library(plyr)

source('../R_functions/processing_prior_pm_clustering.R')
source('../R_functions/Path_Metrics_Clustering_R.R')
source('../R_functions/geometric_perturbation.R')

set.seed(1) #To ensure legends don't overlap data points

# Colorblind-friendly palette

cbfly <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load and revalue labels
true_labels_num <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')
true_labels_num <- true_labels_num[,-1]
true_labels <- as.factor(true_labels_num)
true_labels <- revalue(true_labels,c('1'='HCC827', '2'='H1975' ,'3'= 'H838', '4'='H2228' , '5'= 'A549'))

# Load and process cell mix SCT:
mydata <- read.csv(file = '../Data_after_Imputation/Cellmix_sng_SCT.csv')
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
sctdata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,
                                        var_cutoff=NULL) #performs PCA and denoising

# Cellmix PCA plot:
pcadata = as.data.frame(t(sctdata))
pcadata$true_labels <- true_labels

p1 <- ggplot(pcadata,aes(x=PC1, y=PC2)) +geom_point(aes( color=true_labels))
p1 <- p1 + theme(legend.position = c(0.11, 0.17),aspect.ratio=1)
p1 <- p1 + scale_colour_manual(values=cbfly)
p1

pdf('../Figures/Fig_4_PCA_Cellmix.pdf')
p1
dev.off()

gp_pca<-geometric_perturbation(X=sctdata,U=sctdata,Labels=true_labels_num) #Use all 40 dimensions for PCA tree
hclust_avg <- hclust(as.dist(gp_pca$Cluster_Centroid_Distances), method = 'average')
plot(hclust_avg, labels = c('HCC827', 'H1975' , 'H838', 'H2228' , 'A549'),main="", xlab="", ylab="", sub="")

# Cellmix PM plot:

pm_output <- Path_Metrics_Clustering_R(dataset=sctdata,
                                       n_clust=length(unique(unlist(true_labels))),
                                       p=2,
                                       Fast_mds=TRUE,
                                       silhouette=FALSE)

pmdata <- as.data.frame(t(pm_output$U1[1:2,]))
pmdata$true_labels <- true_labels

p2 <- ggplot(pmdata,aes(x=V1, y=V2)) +geom_point(aes( color=true_labels))
p2 <- p2 + labs(x ='PM1',y='PM2')+theme(legend.position = c(0.11, 0.17),
                                        aspect.ratio=1)
p2 <- p2 + scale_colour_manual(values=cbfly)

pdf('../Figures/Fig_4_PM_Cellmix.pdf')
p2
dev.off()

gp_pm<-geometric_perturbation(X=sctdata,U=pm_output$U1,Labels=true_labels_num)
hclust_avg <- hclust(as.dist(gp_pm$Cluster_Centroid_Distances), method = 'average')
plot(hclust_avg, labels = c('HCC827', 'H1975' , 'H838', 'H2228' , 'A549'),main="", xlab="", ylab="", sub="")

# Load and process (linnorm) cell mix count data:

mydata <- read.csv(file = '../Count_data/cellmix_sng.csv')
rownames(mydata) = mydata[,1]
mydata <- mydata[,-1]
mydata <- Linnorm::Linnorm(mydata)
lindata <- processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,
                                          var_cutoff=NULL,LocalAvgAllPts = FALSE) #performs PCA

# Cell mix umap

umap_embedding2d <- umap::umap(t(lindata),n_components=2) #Use same dimension as for PM
umapdata2d = as.data.frame(umap_embedding2d$layout)
umapdata2d$true_labels <- true_labels

p3 <- ggplot(umapdata2d, aes(x=V1, y=V2)) + geom_point(aes( color=true_labels))
p3 <-p3 + labs(x ='UMAP1',y='UMAP2') + theme(legend.position = c(0.89, 0.83), aspect.ratio=1)
p3 <-p3 + scale_colour_manual(values=cbfly)

pdf('../Figures/Fig_4_UMAP_Cellmix.pdf')
p3
dev.off()

umap_embedding<-umap::umap(t(lindata),n_components=pm_output$r) #Use same dimension as for PM
umapdata = as.data.frame(umap_embedding$layout)

gp_umap<-geometric_perturbation(X=lindata,U=t(umap_embedding$layout),Labels=true_labels_num)
hclust_avg <- hclust(as.dist(gp_umap$Cluster_Centroid_Distances), method = 'average')
plot(hclust_avg, labels = c('HCC827', 'H1975' , 'H838', 'H2228' , 'A549'),main="", xlab="", ylab="", sub="")

# Cell mix t-sne

tsne_embedding2d<-Rtsne::Rtsne(t(lindata),dims=2) #Use same dimension as for PM
tsnedata2d = as.data.frame(tsne_embedding2d$Y)
tsnedata2d$true_labels <- true_labels

p4 <- ggplot(tsnedata2d,aes(x=V1, y=V2)) + geom_point(aes( color=true_labels))
p4 <- p4 + labs(x ='TSNE1',y='TSNE2')+theme(legend.position = c(0.11, 0.17),aspect.ratio=1)
p4 <- p4 + scale_colour_manual(values=cbfly)

pdf('../Figures/Fig_4_Rtsne_Cellmix.pdf')
p4
dev.off()

set.seed(1) 

#tsne_embedding <- tsne::tsne(t(lindata),k=pm_output$r) #Use same dimension as for PM

tsne_embedding <- Rtsne::Rtsne(t(lindata),dims=3) 
tsnedata = as.data.frame(tsne_embedding$Y)

#gp_tsne<-geometric_perturbation(X=lindata,U=t(tsne_embedding),Labels=true_labels_num)

gp_tsne<-geometric_perturbation(X=lindata,U=t(tsne_embedding),Labels=true_labels_num)

hclust_avg <- hclust(as.dist(gp_tsne$Cluster_Centroid_Distances), method = 'average')

plot(hclust_avg, labels = c('HCC827', 'H1975' , 'H838', 'H2228' , 'A549'),
         main="", xlab="", ylab="", sub="")

pdf('../Figures/Fig_4_hierarchical_plot_Rtsne_Cellmix.pdf')

plot(hclust_avg, labels = c('HCC827', 'H1975' , 'H838', 'H2228' , 'A549'),
     main="", xlab="", ylab="", sub="")

dev.off()