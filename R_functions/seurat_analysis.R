#This function processes and clusters an input data set with Seurat

#Parameter explanation
#counts: data path of 10x or counts matrix
#project :project name
#processing,normalization,clustering,signature_genes= if TRUE the corresponding will be done  
# min.cell,min.feat: for Create seurat object
#species= species the data come from
#mt= threshold for mt genes per cell
#n_hvg = number of highly variable genes for vst
#dim= number of pcs for jack straw
#res= resolution for findClusters
#
#The defaults
#processing=TRUE
#normalization=TRUE
#clustering=FALSE
#signature_genes=FALSE
#min.cell=3
#min.feat=200 
#species='human'
#mt= 20
#n_hvg =2000
#dim= 100
#res= 0.5
  

  
seurat_analysis<-function(counts,
                          project ,
                          processing=TRUE,
                          normalization=TRUE,
                          clustering=FALSE,
                          signature_genes=FALSE,
                          min.cell=3,
                          min.feat=200,
                          species='human',
                          mt=20,
                          n_hvg=2000,
                          dim=100,
                          res=0.5)
{
  
 # Input data
 
  
  if(is.character(counts)){
    
    x <- Read10X(data.dir = counts)
    seurat<-CreateSeuratObject( x,min.cells = min.cell,min.features =min.feat,
                                project=project )
    
  }else if(is.matrix(counts) |  is.data.frame(counts)){
    
    seurat<-CreateSeuratObject(counts=counts,min.cells = min.cell,
                               min.features =min.feat,project =project )
    
  }else{
    seurat<-counts
  }
  
  
  
  #Preprocessing
  
  if(processing){  
    
    if(species=='human'){
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    }else{
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")  
    }
    
    seurat <- subset(seurat, subset = percent.mt <mt)
  }
  
  #Normalization and dim reduction
  if(normalization){  
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", 
                                   nfeatures = n_hvg,verbose = FALSE)
   }
    
  #Clustering
  if(clustering){
    seurat <- ScaleData(seurat, features = rownames(seurat))
    seurat<- RunPCA(seurat, features = VariableFeatures(object = seurat))
    
    seurat <- JackStraw(seurat, num.replicate = 100, dims = dim)
    seurat <- ScoreJackStraw(seurat, dims = 1:dim)
   
    #find important pcs
    pcid<-which(seurat@reductions$pca@jackstraw@overall.p.values[,2]<0.001)
    seurat<-RunTSNE(seurat, reduction = "pca", dims=pcid)
    seurat <- FindNeighbors(seurat, dims =pcid,do.plot = FALSE )
    seurat <- FindClusters(seurat, resolution = res,plot.SNN = FALSE,  
                           save.SNN = TRUE)
    seurat <- RunUMAP(seurat, dims =pcid)
    
    }


  # Signature genes
  if(signature_genes){
  seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE)
  top.markers<-seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  seurat<-AddMetaData(seurat,metadata =top.markers,col.name = 'top.markers')
  }
  
  #Assign cell type to clusters
  #if(cell_annotation){
  
  #B = GetAssayData(seurat)
  
  #common <- intersect(rownames(B), rownames(ref))
  #trained <- trainSingleR(ref[common,], labels=ref$label.main)
  #pred.ref <- classifySingleR(B[common,], trained)
  #predicted_Cell_types<-table(pred.ref$labels)
  
  #celltype_per_cell<-pred.ref@listData$labels
  
  #seurat<-AddMetaData(seurat,metadata = celltype_per_cell,col.name = 'cell_types')
  
 # }

    return(seurat)  }

