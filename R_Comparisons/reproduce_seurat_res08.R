source('../R_functions/clustering.R')

dataset=c('RNAMix1','RNAMix2','TMLung','Beta2','TMPanc','BaronPanc',
          'PBMC4k','CellMixSng')
method=c('RaceID_clustering','SC3_clustering','Seurat_clustering_default_res_0.8',
         'Seurat_Clustering_PM_processing','Seurat_def_clustering',
         'SIMLR_clustering','Path_metrics_clustering','kmeans',
         'dbscan',	'umap+dbscan',	'tsne+kmeans')


method<-method[3]

acc_res<-matrix(0,nrow=1,ncol=8)
colnames(acc_res)<-c( "ARI","ECP" , "ECA","number_of_clusters","complete_rt","clustering_rt", "method","dataset_name")
for(d in dataset){
  for(m in method){
    print(d)
    res<-clustering(dataset=d,method=m)
    accuracy<-res$accuracy$accuracy
    runtime<-t(as.data.frame(res$accuracy$runtime))
   
    dataset_name<-d
    result<-cbind(accuracy,runtime,method,dataset_name)
    acc_res<-rbind(acc_res,result)
  }
  
}

saveRDS(acc_res,paste0('../Results/reproduce_',method,'.rds'))

library(writexl)
write_xlsx(acc_res,paste0('../Results/reproduce_',method,'.xlsx'))

