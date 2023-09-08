#Processing of Data sets for Path Metrics clustering ---------------------------

processing_prior_pm_clustering<-function(dataset, # rows with features, columns with points
                                         LogNormalization= TRUE,
                                         var_cutoff=NULL,
                                         EliminateGenes = TRUE,
                                         NumberOfGenes=2000,
                                         LocalAvgAllPts = TRUE,
                                         LocalAvgNbhdSize = 12){
  
# Libraries ---------------------------------------------------------------------
require(igraph)
require(RANN)

# Input, Parmeters, Processing & Denoising Options (Defaults used in paper) ----

Data<-t(dataset)
LogNormalization= LogNormalization

EliminateGenes = EliminateGenes
NumberOfGenes = NumberOfGenes # If EliminateGenes = TRUE, number of genes to keep

LocalAvgAllPts = LocalAvgAllPts
LocalAvgNbhdSize = LocalAvgNbhdSize

NumberOfPCs = 40              #Number of principal components to use

if(LogNormalization){RescaleHighVarGenes = TRUE}else{
  RescaleHighVarGenes = FALSE}
var_cutoff<-var_cutoff

# Normalization ----------------------------------------------------------------

if(LogNormalization){
  Data = diag(1/apply(Data,1,sum))%*%Data*10000
  Data = log(Data+1)}


#Optional: restrict to high variance genes; rescale super high variance genes --

if(EliminateGenes){
  gene_var<-apply(Data,2,function(x) var(x))
  sorted_var = sort(gene_var, decreasing=TRUE)
  cutoff = sorted_var[NumberOfGenes+1]
  high_var_genes = which(gene_var>cutoff)
  X = Data[,high_var_genes] }else{
  X = Data}

# Optional: rescale super high variance genes ----------------------------------

if(RescaleHighVarGenes){
  col.var_X<-apply(X,2,function(x) var(x))
  #See which genes are outliers in terms of very high variance:
  if(is.null(var_cutoff)){
    var_cutoff = quantile(col.var_X, probs = .75)+1.5*quantile(col.var_X, probs =.75)-quantile(col.var_X, probs =.25)}
  
  super_high_var_genes = which(col.var_X>var_cutoff)
  
  var_HG<-apply(X[,super_high_var_genes],2,function(x) var(x))
  X[,super_high_var_genes] = sqrt(var_cutoff)*X[,super_high_var_genes]%*%diag(1/sqrt(var_HG))
}


# Use PCA to reduce dimension --------------------------------------------------

pca_res <- prcomp(X, center = TRUE,scale. =FALSE, rank.=min(NumberOfPCs,ncol(X)))
X = pca_res$x

# Denoise Data (optional) -----------------------------------------------------

if(LocalAvgAllPts){
  knn=nn2(data= X,query=X, k = LocalAvgNbhdSize)
  D_KNN=as.data.frame(knn$nn.dists)
  IDX=knn$nn.idx
  
  LocalAverages =  matrix(0, nrow=nrow(X),ncol=ncol(X))
  
  for(i in 1:nrow(X)){
   LocalAverages[i,] = apply(X[IDX[i,],],2,function(x) mean(x))}
  
  X_before_local_avg = X
  X = LocalAverages
  #Redo PCA --------------------------------------------------
  
  pca_res <- prcomp(X, center = TRUE,scale. =FALSE, rank.=ncol(X))
  X = pca_res$x
  
  
}

return(t(X))
}
