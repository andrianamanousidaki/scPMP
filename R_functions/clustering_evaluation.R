## Clustering evaluation function
## Given predicted and true cluster labels as input
## we get ECA,ECP,ARI, number of clusters as output

clustering_evaluation<-function(predicted_clusters,true_clusters){
  library(plyr)
  library(cluster)
  library(mclust)
  entropy=function(x){
    freqs <- table(x)/length(x)
    freqs = freqs[freqs>0]
    return(-sum(freqs * log(freqs)))
  }
  
  ecp=mean(unlist(lapply(unique(true_clusters),function(x){entropy(predicted_clusters[true_clusters==x])})))
  
  eca=mean(unlist(lapply(unique(predicted_clusters),function(x){entropy(true_clusters[predicted_clusters==x])})))
  
  evaluation <- list(
    ARI=adjustedRandIndex(true_clusters,predicted_clusters),
    ECP=ecp,
    ECA=eca,
    number_of_clusters=length(table(predicted_clusters))
  )
  
  return(evaluation)
}