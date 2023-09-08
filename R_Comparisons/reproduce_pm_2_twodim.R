source('../R_functions/clustering.R')

dataset=c('RNAMix1','RNAMix2','TMLung','Beta2','TMPanc','BaronPanc',
          'PBMC4k','CellMixSng','Balls','Enlongated_With_Bridge',
          'Swiss','SO(3)')
method=c('RaceID_clustering','SC3_clustering','Seurat_clustering_default_res_0.8',
         'Seurat_Clustering_PM_processing','Seurat_def_clustering',
         'SIMLR_clustering','Path_metrics_clustering','kmeans',
         'dbscan',	'umap+dbscan',	'tsne+kmeans')


method<-method[7]

acc_res<-matrix(0,nrow=1,ncol=9)
colnames(acc_res)<-c("ARI", "ECP", "ECA",  "method",  "complete_rt","clustering_rt",
                     "gp", "d", "nclust" )

for(d in dataset){
  for(m in method){
    print(d)
    res<-clustering(dataset=d,method=m,p=2,geometric_perturbation=TRUE,
                    two_dimensions=TRUE, predict_pm_clusters=FALSE)
    accuracy<-as.data.frame(res$accuracy)
    runtime<-as.data.frame(res[6:7])
    gp<-res$Geometric_Perturbation
    nclust<-length(unique(res$pm_labels))
    result<-cbind(accuracy,runtime,gp,d,nclust)
    acc_res<-rbind(acc_res,result)
  }
  
}

saveRDS(acc_res,paste0('../Results/reproduce_',method,'_2_2d.rds'))

library(writexl)
write_xlsx(acc_res,paste0('../Results/reproduce_',method,'_2_2d.xlsx'))

