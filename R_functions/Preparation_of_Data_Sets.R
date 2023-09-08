# This script presents how we downloaded, filtered and imputed the data sets
# used for comparisons. The resulted data sets are saved in the folders
# Count_data and Data_after_Imputation.


# Libraries -------------------------------------------------------------------


library(SAVER)
library(dplyr)
library(Seurat)
library(patchwork)
source("../R_functions/seurat_analysis.R")


# Create directory of pre-processed data sets ----------------------------------

folder1 <- c('../Data_after_Imputation')
folder2 <- c( '../Count_data')

#Make  folder of imputed data

if (file.exists(folder1)) {
  
  cat("The folder already exists")
  
}else{
  
  dir.create(folder1)
  
}

#Make  folder of cleaned count data

if (file.exists(folder2)) {
  
  cat("The folder already exists")
  
}else{
  
  dir.create(folder2)
  
}

################################################################################
# RNAmix1 and RNAmix2---------------------------------------------------------
################################################################################

#Download data from 'https://github.com/LuyiTian/sc_mixology/tree/master/data'

url='https://github.com/LuyiTian/sc_mixology/raw/master/data/mRNAmix_qc.RData'
download.file(url,'mRNAmix_qc.RData')
load('mRNAmix_qc.RData')

# The data are already cleaned

RNAmix1<-sce8_qc
l1<-paste(sce8_qc@colData$H2228_prop,sce8_qc@colData$H1975_prop,sce8_qc@colData$HCC827_prop,sep="_")
c1<-sce8_qc@colData$cell_name
l1<-cbind(c1,l1)

RNAmix2<-sce2_qc
l2<-paste(sce2_qc@colData$H2228_prop,sce2_qc@colData$H1975_prop,sce2_qc@colData$HCC827_prop,sep="_")
c2<-sce2_qc@colData$cell_name
l2<-cbind(c2,l2)

rnamix1<-as.Seurat(RNAmix1,data=NULL)
saveRDS(GetAssayData(rnamix1),"../Count_data/rnamix1_filtered.rds")  

rnamix2<-as.Seurat(RNAmix2,data=NULL)
saveRDS(GetAssayData(rnamix2),"../Count_data/rnamix2_filtered.rds")  

# Imputation

rnamix1_saver<-saver(GetAssayData(rnamix1),ncores=10,size.factor=1)
write.csv(rnamix1_saver$estimate,"../Data_after_Imputation/RNAmix1_original.csv")

rnamix2_saver<-saver(GetAssayData(rnamix2),ncores=10,size.factor=1)
write.csv(rnamix2_saver$estimate,"../Data_after_Imputation/RNAmix2_original.csv")

library(plyr)
lr1<-revalue(as.factor(l1[,2]), c( "0.68_0.16_0.16"="1",
                                   "0.33_0.33_0.33"="2",
                                   "0.16_0.68_0.16"="3", 
                                   "1_0_0"="4",
                                   "0_1_0"="5",   
                                   "0_0_1" ="6", 
                                   "0.16_0.16_0.68"="7" ))

write.csv(lr1, "../Data_after_Imputation/rnamix1_original_labels.csv")

lr2<-revalue(as.factor(l2[,2]), c( "0.68_0.16_0.16"="1",
                                   "0.33_0.33_0.33"="2",
                                   "0.16_0.68_0.16"="3",
                                   "1_0_0"="4",
                                   "0_1_0"="5",   
                                   "0_0_1" ="6", 
                                   "0.16_0.16_0.68"="7" ))
write.csv(lr2, "../Data_after_Imputation/rnamix2_original_labels.csv")

##################################################################################
#Lung Tabula Muris - no quality control filtering is needed -------------------
################################################################################

# Download data file "https://ndownloader.figshare.com/files/10700143"

url<-"https://ndownloader.figshare.com/files/10700143/FACS.zip"
download.file(url,'FACS.zip')
unzip('FACS.zip')

download.file(url = 'https://figshare.com/ndownloader/files/13088129/annotations_facs.csv','annotations_facs.csv')
annotation<-read.csv('annotations_facs.csv')

tm_lung <-read.csv('FACS/Lung-counts.csv')
rownames(tm_lung)<-tm_lung[,1]
tm_lung<-tm_lung[,-1]

#SAVER is computationally intensive for large datasets, we recommend 
#running it in parallel on a cluster either by creating your own parallel
#environment or by specifying `ncores`.

tm_lung_saver<-saver(as.matrix(tm_lung),ncores=8,size.factor=1)
tm_lung_saver<-tm_lung_saver$estimate
saveRDS(tm_lung_saver, '../Data_after_Imputation/lung_saver.rds')


counts<-readRDS("lung_saver.rds")
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]

write.csv(counts,'../Data_after_Imputation/lung_tabula_muris_saver.csv')


labels[,2]<-revalue(labels[,2],c("alveolar epithelial type 1 cells, alveolar epithelial type 2 cells, club cells, and basal cells"=
                                   "alveolar ep 1&2,club & basal",
                                 "dendritic cells, alveolar macrophages, and interstital macrophages"=
                                   "dendritic,alveolar macro,interstital macro"))
num_labels<-revalue(labels[,2],c("lung neuroendocrine cells and unknown cells"="1",
                                 "alveolar ep 1&2,club & basal"  ="2"             ,
                                 "mast cells and unknown immune cells"="3"        ,
                                 "circulating monocytes" = "4"                      ,
                                 "invading monocytes"  = "5"                       ,
                                 "multiciliated cells"  = "6"                      ,
                                 "dendritic,alveolar macro,interstital macro" ="7"))

write.csv(num_labels,"../Data_after_Imputation/lung_tabula_muris_labels.csv")

cores<-cbind(1:7,c("lung neuroendocrine cells and unknown cells",
                   "alveolar ep 1&2,club & basal"          ,
                   "mast cells and unknown immune cells"      ,
                   "circulating monocytes"                       ,
                   "invading monocytes"                     ,
                   "multiciliated cells"                       ,
                   "dendritic,alveolar macro,interstital macro" ))
write.csv(cores,'../Data_after_Imputation/num_to_word_tabula_muris_lung.csv')

# Prep for RACEID3 ---------------------------------------------------------------

mydata<-read.csv('FACS/Lung-counts.csv')
rownames(mydata)<-mydata[,1]
mydata<-mydata[,-1]
labels<-read.csv('../Data_after_Imputation/lung_tabula_muris_labels.csv')

counts<-mydata
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]
#load('facs_Pancreas_seurat_tiss.Robj',verbose = TRUE)
write.csv(counts,'../Count_data/lung_counts_updated.csv')

types<-unique(labels[,2])

num_labels<-revalue(labels[,2],c("lung neuroendocrine cells and unknown cells"="1" ,                                                   
                                 "alveolar epithelial type 1 cells, alveolar epithelial type 2 cells, club cells, and basal cells"="2",
                                 "mast cells and unknown immune cells"="3",                                                            
                                 "circulating monocytes"="4",                                                                          
                                 "invading monocytes" ="5",                                                                            
                                 "multiciliated cells"="6" ,                                                                           
                                 "dendritic cells, alveolar macrophages, and interstital macrophages" ="7"))
labels[,2]<-num_labels
rownames(lables)<-labels[,1]
labels<-labels[,-1]

write.csv(labels,"../Count_data/lung_counts_updated_labels.csv")
################################################################################
#Beta Simulated Dataset --------------------------------------------------------
################################################################################


print_qc_plots = FALSE
source('../R_functions/simulation_of_beta_data_set.R')


################################################################################
#Pancreas Tabula Muris - no quality control filtering is needed---------------
################################################################################

# Download data file "https://ndownloader.figshare.com/files/10700143"
url<-"https://ndownloader.figshare.com/files/10700143/FACS.zip"
download.file(url,'FACS.zip')
unzip('FACS.zip')
tm_panc <-read.csv('FACS/Pancreas-counts.csv')
rownames(tm_panc)<-tm_panc[,1]
tm_panc<-tm_panc[,-1]
tm_panc_saver<-saver(tm_panc,ncores=20,size.factor=1)
tm_panc_saver<-tm_panc_saver$estimate
saveRDS(tm_panc_saver, 'pancreas_saver.rds')

# labels
download.file(url = 'https://figshare.com/ndownloader/files/13088129/annotations_facs.csv','annotations_facs.csv')
annotation<-read.csv('annotations_facs.csv')

counts<-readRDS('pancreas_saver.rds')
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]
write.csv(counts,'../Data_after_Imputation/pancreas_tabula_muris_saver.csv')

labels[,2]<-revalue(labels[,2],c("acinar cell"="1",
                                 "beta cell"="2",
                                 "stellate cell"="3",
                                 "pancreatic A cell"="4",
                                 "ductal cell"="5",
                                 "pancreatic D cell"="6",
                                 "pancreatic PP cell"="7"))


write.csv(num_labels,"../Data_after_Imputation/pancreatic_tabula_muris_labels.csv")

# Prep for RACEID3 ------------------------------------------------------------

tm_panc <-read.csv('FACS/Pancreas-counts.csv')
rownames(tm_panc)<-tm_panc[,1]
tm_panc<-tm_panc[,-1]

##labels

counts<-tm_panc
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]
write.csv(counts,'../Count_data/pancreas_tm_counts_updated.csv')

library(plyr)
labels[,2]<-revalue(labels[,2],c("acinar cell"="1",
                                 "beta cell"="2",
                                 "stellate cell"="3",
                                 "pancreatic A cell"="4",
                                 "ductal cell"="5",
                                 "pancreatic D cell"="6",
                                 "pancreatic PP cell"="7"))


write.csv(labels,"../Count_data/pancreatic_tm_updated_labels.csv")





################################################################################
#Baron's Pancreatic data--------------------------------------------------------
################################################################################
#Download data from 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230757'

url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2230757&format=file&file=GSM2230757%5Fhuman1%5Fumifm%5Fcounts%2Ecsv%2Egz/GSM2230757_human1_umifm_counts.csv.gz'
download.file(url = url, 'GSM2230757_human1_umifm_counts.csv.gz')

x <- read.csv("GSM2230757_human1_umifm_counts.csv.gz", header = TRUE, check.names = FALSE)
cell_clusters<-x[,2:3]
rownames(cell_clusters)<-make.names(x[,2], unique=TRUE)
cell_clusters[,1]<-rownames(cell_clusters)

x.dat <- t(as.matrix(x[, 4:ncol(x)]))
colnames(x.dat) <- make.names(x[,2], unique=TRUE)
x <- x.dat

saveRDS(x,'baron_pancreatic_counts_before_filtering.rds')

filtered_data<-seurat_analysis(counts=x,processing =TRUE,project="Baron's Pancreatic data",
                               normalization = FALSE,clustering =FALSE,
                               signature_genes = FALSE)

data_for_saver<-GetAssayData(filtered_data)
saveRDS(data_for_saver, 'pancreatic_for_saver.rds')

panc_saver<-saver(data_for_saver,ncores=20,size.factor=1)
wrte.csv(panc_saver,'pancreatic_saver_no_normalization.csv')

# SCT normalization
Data <- CreateSeuratObject(counts = panc_saver$estimate,min.cell=0,min.feat=0)

#For tabula muris data sets: pattern="^mt-"
Data <- PercentageFeatureSet(Data, pattern = "^MT-", col.name = "percent.mt")
Data <- SCTransform(Data, vars.to.regress = "percent.mt")
tosave<-GetAssayData(Data,assay="SCT")
genes<-Data@assays$SCT@var.features

#Choose how to save output data set

write.csv(tosave[genes,],'Baron_Pancreatic_SCT.csv')

# Remove tiny clusters

mydata <- read.csv(file = 'Baron_Pancreatic_SCT.csv')
true_labels <- read.csv(file = 'pancreatic_num_labels.csv')
rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels[,1]<-colnames(mydata)

freq<-data.frame(table(true_labels[,2]))
ord<-order(freq[,2],decreasing = TRUE)
ord<-ord[1:8]
labels<-freq[ord,1]
true_labels<-true_labels[which(true_labels[,2] %in% labels),]

mydata<-mydata[,true_labels[,1]]
write.csv(mydata,'../Data_after_Imputation/BaronPancSCT_filtered.csv')

rownames(true_labels) = true_labels[,1]
true_labels<-true_labels[,-1]
write.csv(true_labels,"../Data_after_Imputation/pancreatic_num_labels_filtered.csv")

# Prep for RACEID3 ------------------------------------------------------------
original<-readRDS('baron_pancreatic_counts_before_filtering.rds')

pancreatic<-read.csv('../Data_after_Imputation/BaronPancSCT_filtered.csv')
rownames(pancreatic)<-pancreatic[,1]
pancreatic<-pancreatic[,-1]

original<-original[,colnames(pancreatic)]
saveRDS(original, '../Count_data/barons_pancreatic_for_Raceid3.rds')
################################################################################
# PBMC4k-------------------------------------------------------------------------
################################################################################
#Download data from 'https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k'
#and SingleR annotation from 'https://github.com/dviraran/SingleR/blob/master/manuscript_figures/FiguresData/SingleR.PBMC.4K.RData'

url='https://github.com/dviraran/SingleR/raw/master/manuscript_figures/FiguresData/SingleR.PBMC.4K.RData'
download.file(url,'SingleR.PBMC.4K.RData')

url='https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz'
download.file(url,'pbmc4k_filtered_gene_bc_matrices.tar.gz')
untar("pbmc4k_filtered_gene_bc_matrices.tar.gz")

## Find cell names annotated by SingleR

load('SingleR.PBMC.4K.RData')
singler2 = singler$singler[[2]]
sr_cells<-singler2$SingleR.single$cell.names

#Prepare data set for SAVER
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/GRCh38/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc4k", min.cells = 3, min.features = 200)
cell <- colnames(pbmc)
cell<-substr(cell,1,nchar(cell)-2)
id<-which(cell%in%sr_cells)
pbmc <- subset(pbmc, cells=colnames(pbmc)[id])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = percent.mt < 20)
saveRDS(GetAssayData(pbmc),'pbmc4k_for_saver.rds')

#Imputation
c<-readRDS("pbmc4k_for_saver.rds")
set.seed(1)
saverc<- saver(c,ncores=20)
saveRDS(saverc$estimate,'pbmc4k_qc_saver.rds')

data<-readRDS('pbmc4k_qc_saver.rds')
write.csv(data,'pbmc4k_qc_saver.csv')

#Labels

pbmc_cells<-colnames(data)
pbmc_cells<-substr(pbmc_cells,1,nchar(pbmc_cells)-2)

main_labels<-data.frame(cells=singler2$SingleR.single.main$cell.names,labels=singler2$SingleR.single.main$labels)
id_cells<-c()

for( i in 1:length(pbmc_cells)){
  id_cells[i]<-which(main_labels[,1]== pbmc_cells[i])
}

main_labels<-main_labels[id_cells,]
library(plyr)

num_main_labels<-revalue(main_labels[,2],c("Monocytes"="1","CD8+ T-cells"="2","CD4+ T-cells"="3",
                                           "B-cells"="4","HSC"="5","NK cells"="6","DC"="7"))

main_labels$num_labels=num_main_labels

#Check that colnames in data set are in correct order
identical(as.character(main_labels$cells),as.character(pbmc_cells))
identical(as.character(main_labels$cells),as.character(rownames(main_labels)))
main_labels<-main_labels[pbmc_cells,]

#Word to numerical labels
d<-data.frame(levels=unique(main_labels[,3]),num=1:7)
write.csv(d,row.names = FALSE,'pbmc4k_num_to_main_labels.csv')
write.csv(num_main_labels,row.names = FALSE,'pbmc4k_main_labels.csv')

#Remove cluster 5,7 
l<-read.csv('pbmc4k_main_labels.csv')
identical(as.character(l$x),as.character(main_labels$num_labels))
rownames(l)<-colnames(data)

ce<-which(l$x!=5 & l$x!=7)
data<-data[,rownames(l)[ce]]
colnames(data)<-substr(colnames(data),1,nchar(colnames(data))-2)
write.csv(data,'../Data_after_Imputation/subset_pbmc4k.csv')

#Merge the two T cell cell types
l$x<-revalue(as.character(l$x),c("3"="2"))

write.csv(l[ce,1],'../Data_after_Imputation/subset_pbmc4k_labels.csv')

# Prep for Raceid3 -------------------------------------------------------------

load('SingleR.PBMC.4K.RData')
singler2 = singler$singler[[2]]
sr_cells<-singler2$SingleR.single$cell.names
pbmc <- Read10X(data.dir = "filtered_gene_bc_matrices/GRCh38/")

cell <- colnames(pbmc)
cell<-substr(cell,1,nchar(cell)-2)
id<-which(cell%in%sr_cells)
pbmc <- pbmc[,id]

saveRDS(pbmc,'pbmc4k_before_mt.rds')


pbmc_seurat <- CreateSeuratObject(counts = pbmc, project = "pbmc4k", min.cells = 3, min.features = 200)
cell <- colnames(pbmc_seurat)
cell<-substr(cell,1,nchar(cell)-2)
id<-which(cell%in%sr_cells)
pbmc_seurat <- subset(pbmc_seurat, cells=colnames(pbmc_seurat)[id])
pbmc_seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc_seurat, pattern = "^MT-")
pbmc_seurat <- subset(pbmc_seurat, subset = percent.mt < 20)

mt_affected_cells<-colnames(pbmc_seurat)


pbmc<-readRDS('pbmc4k_before_mt.rds')
pbmc<-pbmc[,mt_affected_cells]

saveRDS(pbmc,'pbmc4k_only_mt.rds')
pbmc<-readRDS('pbmc4k_only_mt.rds')
#mainlabels
pbmc_cells<-colnames(pbmc)
pbmc_cells<-substr(pbmc_cells,1,nchar(pbmc_cells)-2)

main_labels<-singler2$SingleR.single.main$labels

id_cells<-c()

for( i in 1:length(pbmc_cells)){
  id_cells[i]<-which(rownames(main_labels)== pbmc_cells[i])
}

main_labels<-main_labels[id_cells]
library(plyr)

#merge the labels
num_main_labels<-revalue(main_labels,c("Monocytes"="1","CD8+ T-cells"="2","CD4+ T-cells"="2",
                                       "B-cells"="4","HSC"="5","NK cells"="6","DC"="7"))
saveRDS(num_main_labels,'pbmc4k_only_mt_7_MAIN_LABELS.rds')
l<-as.data.frame(num_main_labels)

rownames(l)<-colnames(pbmc)
#l<-l[,-1]
ce<-which(l!=5 & l!=7)
pbmc<-pbmc[,rownames(l)[ce]]
write.csv(pbmc,'../Count_data/subset_pbmc4k_only_mt.csv')

l2<-as.data.frame(l[rownames(l)[ce],])
write.csv(l[rownames(l)[ce],1],'../Count_data/pbmc4k_only_mt_5_MAIN_LABELS.csv')


################################################################################
# CellmixSNG -------------------------------------------------------------------
################################################################################
##No imputation was applied to CELLMIXSNG data set 

#Download data from  'https://github.com/LuyiTian/sc_mixology/tree/master/data'
url='https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class_5cl.RData'
download.file(url,'sincell_with_class_5cl.RData')
load("sincell_with_class_5cl.RData")

url='https://github.com/LuyiTian/sc_mixology/raw/master/data/csv/sc_10x_5cl.metadata.csv.gz'
download.file(url,'sc_10x_5cl.metadata.csv.gz')

## Find singletons of CELLMIX

data<-sce_sc_10x_5cl_qc
metadata<-read.csv('sc_10x_5cl.metadata.csv.gz')
counts<-data@assays$data$counts
sng_index<-data.frame(rownames(metadata),metadata$demuxlet_cls)
rownames(sng_index)<-sng_index[,1]

counts<-counts[,which(sng_index[,2]=='SNG')]
dim(counts)
write.csv(counts,"../Count_data/cellmix_sng.csv")


## Filtered labels

info<-sce_sc_10x_5cl_qc@colData@listData
cell_labels<-info$cell_line_demuxlet
c_labels<-cbind(colnames(sce_sc_10x_5cl_qc), cell_labels)
rownames(c_labels)<-c_labels[,1]
filtered_labels<-c_labels[which(sng_index[,2]=='SNG'),]
library(plyr)
revallab<-revalue(filtered_labels[,2], c("HCC827"=1, "H1975"=2 , "H838"=3 ,  "H2228"=4 , "A549"=5))
write.csv(revallab,"../Data_after_Imputation/cellmix_sng_labels.csv")

#   SCT normalization

Data <- CreateSeuratObject(counts = counts,min.cell=0,min.feat=0)

#For tabula muris data sets: pattern="^mt-"
Data <- PercentageFeatureSet(Data, pattern = "^MT-", col.name = "percent.mt")
Data <- SCTransform(Data, vars.to.regress = "percent.mt")
tosave<-GetAssayData(Data,assay="SCT")
genes<-Data@assays$SCT@var.features

#Choose how to save output data set

write.csv(tosave[genes,],'../Data_after_Imputation/Cellmix_sng_SCT.csv')

