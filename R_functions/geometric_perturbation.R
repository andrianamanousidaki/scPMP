geometric_perturbation<-function(X,U,Labels){

 # input: data set X
 # input: data set U, which is a low dimensional embedding of X pointsXfeatures
 # input: Labels, predetermined groups in X
 # output: measurement of global geometric perturbation of the groups when
 # data forced into lower dimensions
X=t(X)
U=t(U)
Labels=unlist(Labels)
k=length(unique(Labels))
uniqueGroups = sort(unique(Labels))
EmbeddingCentroids = matrix(0,nrow=k,ncol=dim(U)[2] )
ClusterCentroidDistances = matrix(0,nrow=k,ncol=k)
ClusterCentroidDistancesFull = matrix(0,nrow=dim(U)[1],ncol=dim(U)[1])
PMClusterCentroidDistances = matrix(0,nrow=k,ncol=k)
PMClusterCentroidDistancesFull = matrix(0,nrow=dim(U)[1],ncol=dim(U)[1])



for( i in 1:(k-1)){
  ind_i = which(Labels==uniqueGroups[i])
  centroid_i = colMeans(X[ind_i,])
  EmbeddingCentroids[i,] = colMeans(U[ind_i,])
  
  for( j in (i+1):k){
    ind_j = which(Labels==uniqueGroups[j])
    centroid_j = colMeans(X[ind_j,])
    pm_centroid_j = colMeans(U[ind_j,])
    ClusterCentroidDistances[i,j] = norm(centroid_i - centroid_j, type = "2") 
    ClusterCentroidDistances[j,i] = ClusterCentroidDistances[i,j]
    ClusterCentroidDistancesFull[ind_i,ind_j]=ClusterCentroidDistances[i,j]
    ClusterCentroidDistancesFull[ind_j,ind_i]=ClusterCentroidDistances[i,j]
    PMClusterCentroidDistances[i,j] = norm(EmbeddingCentroids[i,] - pm_centroid_j, type = "2") 
    PMClusterCentroidDistances[j,i] = PMClusterCentroidDistances[i,j]
    PMClusterCentroidDistancesFull[ind_i,ind_j]=PMClusterCentroidDistances[i,j]
    PMClusterCentroidDistancesFull[ind_j,ind_i]=PMClusterCentroidDistances[i,j]
  }
}

c = sum(rowSums( ClusterCentroidDistancesFull*PMClusterCentroidDistancesFull ))/ sum(rowSums( PMClusterCentroidDistancesFull^2 ))
Geometric_Perturbation = sum(rowSums( (ClusterCentroidDistancesFull - c*PMClusterCentroidDistancesFull)^2 ))/sum(rowSums( ClusterCentroidDistancesFull^2))

output<-list(Geometric_Perturbation =Geometric_Perturbation , Cluster_Centroid_Distances=PMClusterCentroidDistances)

return(output)
}