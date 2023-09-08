Seurat_clustering_PM_processing<-function(dataset=c('RNAMix1','RNAMix2','TMLung','Beta2',
                                                  'TMPanc','BaronPanc','PBMC4k','CellMixSng',
                                                  'Balls','Enlongated_With_Bridge',
                                                  'Swiss','SO(3)')){
  
  #Required Libraies ----------------------------------------------------------
  s1<-proc.time()
  
  require(Seurat)
  
  source('../R_functions/processing_prior_pm_clustering.R')

  #Input dataset ---------------------------------------------------------------
  
  
  if( dataset =='RNAMix1'){
    mydata <- read.csv(file = '../Data_after_Imputation/RNAmix1_original.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/rnamix1_original_labels.csv')
    var_cutoff=0.3
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  
  if( dataset =='RNAMix2'){
    mydata <- read.csv(file = '../Data_after_Imputation/RNAmix2_original.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/rnamix2_original_labels.csv')
    var_cutoff=0.2
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  
  if( dataset =='TMLung'){
    mydata <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_saver.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/lung_tabula_muris_labels.csv')
    var_cutoff=2
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  if( dataset =='Beta2'){
    mydata <- read.csv(file = '../Data_after_Imputation/beta_3_4_10_filtered_saver.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/beta_cell_groups_3_4_10_after_filtering.csv')
    true_labels<-true_labels[,-c(1,2)]
    var_cutoff=0.4
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  if(dataset =='TMPanc'){
    mydata <- read.csv(file = '../Data_after_Imputation/pancreas_tabula_muris_saver.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_tabula_muris_labels.csv')
    var_cutoff=1
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  if(dataset =='BaronPanc'){
    mydata <- read.csv(file = '../Data_after_Imputation/BaronPancSCT_filtered.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/pancreatic_num_labels_filtered.csv')
    var_cutoff=NULL
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
   }
  
  if(dataset =='PBMC4k'){
    mydata <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k.csv')
    true_labels <- read.csv(file = '../Data_after_Imputation/subset_pbmc4k_labels.csv')
    var_cutoff=2
    
    
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    }
  
  if(dataset =='CellMixSng'){
    mydata <- read.csv(file = '../Data_after_Imputation/Cellmix_sng_SCT.csv')
    true_labels <- read.csv('../Data_after_Imputation/cellmix_sng_labels.csv')  
 
    rownames(mydata) = mydata[,1]
    mydata<-mydata[,-1]
    true_labels<-true_labels[,-1]
    var_cutoff=NULL
    }
  # Manifold datasets
  
  if(dataset =='SO(3)'){
    mydata <- read.csv(file = '../Count_data/GLmanifold_d9_N3000_k3_sig0075.csv',header = FALSE)
    true_labels <- read.csv('../Count_data/GLmanifold_d9_N3000_k3_sig0075_Labels.csv',header = FALSE)
    mydata <-t(mydata)
    var_cutoff=NULL}
  
  if(dataset =='Balls'){
    mydata <- read.csv(file = '../Count_data/Balls_n400_r53.csv',header = FALSE)
    true_labels <- read.csv('../Count_data/Balls_n400_r53_Labels.csv',header = FALSE)
    mydata <-t(mydata)
    var_cutoff=NULL}
  
  if(dataset =='Enlongated_With_Bridge'){
    mydata <- read.csv(file = '../Count_data/ElongatedGaussiansWithBridge3.csv',header = FALSE)
    true_labels <- read.csv('../Count_data/ElongatedGaussiansWithBridge3_Labels.csv',header = FALSE)
    mydata <-t(mydata)
    var_cutoff=NULL}
  
  if(dataset =='Swiss'){
    mydata <- read.csv(file = '../Count_data/SwissRoll1.csv',header = FALSE)
    true_labels <- read.csv('../Count_data/SwissRoll1_Labels.csv',header = FALSE)
    mydata <-t(mydata)
    var_cutoff=NULL}
  

  
  # Process data set : 
  set.seed(20)
  if(dataset=='BaronPanc' || dataset=='CellMixSng'){ 
    
    mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,var_cutoff=var_cutoff)}else{
   
      if(dataset!='Balls'&& dataset!='Enlongated_With_Bridge'&&dataset!='Swiss' && dataset!='SO(3)'){
        mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=TRUE,var_cutoff=var_cutoff)}else{
          mydata<-processing_prior_pm_clustering(dataset=mydata,LogNormalization=FALSE,EliminateGenes =FALSE,var_cutoff=var_cutoff)}
          
     }
  
  # Create Seurat Object ---------------------------------------------------------
  rownames(mydata) = 1:nrow(mydata)
  colnames(mydata) =paste0("Cell_", 1:ncol(mydata)) 
  Seurat_mydata <- CreateSeuratObject(counts = mydata, project = "mydata", min.cells = 0, min.features = 0)
  
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(mydata))
  Seurat_mydata@assays$RNA@scale.data <- as.matrix(mydata)
  
  rownames(mydata) <- paste0("PC_", 1:min(nrow(mydata),40))
  #Seurat_mydata <- RunPCA(Seurat_mydata, features = rownames(Seurat_mydata),return.only.var.genes = FALSE, npcs = 2)
  Seurat_mydata[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(t(mydata)), key = "PC_",assay='RNA')
  
  #Find Neighbors ---------------------------------------------------------------
  s2<-proc.time()
  
  Seurat_mydata <- FindNeighbors(Seurat_mydata, dims = 1:min(nrow(mydata),40))
  
  # Add labels in Seurat Object -------------------------------------------------
  
  colnames(true_labels)=NULL
  #rownames(true_labels) <- colnames(x =Seurat_mydata)
  Seurat_mydata<-AddMetaData(Seurat_mydata, true_labels, col.name = 'true_labels')
  
  # Clustering ---------------------------------------------------------------
  
  
  if( dataset =='RNAMix1'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.11)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.11}
  
  if( dataset =='RNAMix2'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.15)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.15}
  
  if( dataset =='TMLung'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.2)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.2}

  if( dataset =='Beta2'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.1)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.1}
  
  if(dataset =='TMPanc'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.02)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.02}
  
  if(dataset =='BaronPanc'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.03)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.03}

  if(dataset =='PBMC4k'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.003)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0.003}
  
  if(dataset =='CellMixSng'){
  Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0)
  Seurat_labels <- Seurat_mydata$RNA_snn_res.0}
  
  if(dataset =='SO(3)'){
    Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
    Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05}
  
  if(dataset =='Balls'){
    Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
    Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05}
  
  if(dataset =='Enlongated_With_Bridge'){
    Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.03)
    Seurat_labels <- Seurat_mydata$RNA_snn_res.0.03}
  
  if(dataset =='Swiss'){
    Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
    Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05}
  
  end<-proc.time()
  runtime<-data.frame(complete_rt=(end-s1)[3]/60,clustering_rt=(end-s2)[3]/60)
  
  
  # Evaluation ------------------------------------------------------
  eval<-clustering_evaluation(Seurat_labels,t(true_labels))
  eval= t(data.frame(unlist(eval)))
  rownames(eval)<-dataset
  accuracy<-cbind(eval,runtime)
  return(accuracy)
}
