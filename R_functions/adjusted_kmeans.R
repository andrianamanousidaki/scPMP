adjusted_kmeans<-function(U,k){
  #this function needs library(pracma)
  km<-kmeans(U,k,nstart=20)
  pm_labels<-km$cluster
  H=table(pm_labels)
  
  if(min(H)< sqrt(nrow(U)/2)){
    num_tiny_clusters=length(which(H<sqrt(nrow(U)/2)))
    km2<- kmeans(U,k+num_tiny_clusters,nstart=20)
    pm_labels_original<-pm_labels
    pm_labels<-km2$cluster
    C<-km2$centers
    H<-table(pm_labels)
    h2<-sort.int(H,index.return=TRUE)
    H2<-h2$x
    idxH<-h2$ix
    tiny_cluster_idx = idxH[1:num_tiny_clusters]
    CD =squareform(as.vector(dist(C, method = "euclidean")))# P CENTROID DISTANCE MATRIX 
    options = CD[,tiny_cluster_idx]
    
    if(length(tiny_cluster_idx)==1){
      v=options
      merge_idx = which.min(v[v>0])
      pm_labels[pm_labels==tiny_cluster_idx] = merge_idx
    }else{
      for (i in 1:length(tiny_cluster_idx)){
      v=options[,i]
      merge_idx = which( options[,i] == min(v[v>0]))
      pm_labels[pm_labels==tiny_cluster_idx[i]] = merge_idx}
      }
    
    
    remaining_labels = unique(pm_labels)
    remaining_labels = remaining_labels[order(remaining_labels)]
     for (i in 1:length(remaining_labels)){
    pm_labels[pm_labels==remaining_labels[i]]=i
     }
    
  }
   
  klabels=pm_labels
  return(klabels)
}