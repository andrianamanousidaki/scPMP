# a = list.files(path ='Results/',pattern = '.rds')
# files = file.path('Results', a)
# all(file.exists(files))
# 
# 
# for(i in 14:16){
#   res<-readRDS(files[i])
#   res<-cbind(res,res[,4])
#   res<-res[,-4]
#   colnames(res)[7:8]<-c("d","nclust")
#   saveRDS(res,files[i])
# }

a = list.files(path ='../Results/',pattern = '.rds')
files = file.path('../Results', a)
all(file.exists(files))

id1 = which(a == "downsampling_res.rds" )
id2 = which( a == "all_runtime_results.rds" )
a = a[-c(id1,id2)] 
files = files[-c(id1,id2)]

merged_time<-as.data.frame(matrix(0,ncol=4,nrow=1))
colnames(merged_time)<-c("method","complete_rt","clustering_rt","dataset")

for(i in 1:(length(a))){
  print(a[i])
  print(files[i])
  res<-as.data.frame(readRDS(files[i]))
  name<-substr(a[i],11,nchar(a[i])-4)
  res$method<-rep(name,times=nrow(res))
  subset_res<-cbind(res$method,res$complete_rt,res$clustering_rt,res$d)
  colnames(subset_res)<-c("method","complete_rt","clustering_rt","dataset")
  merged_time<-rbind(merged_time,subset_res)
}
scanpy<-readRDS('../Results/scanpy_ecp_eca_ari_RESULTS.rds')
method<-rep('scanpy',times=8)
complete_rt<-scanpy$complete_rt
dataset<-scanpy$dataset
dataset[6]<-'BaronPanc'
clustering_rt<-scanpy$clustering_rt

scanpy_data<-data.frame(method=method,complete_rt=complete_rt,clustering_rt=clustering_rt,dataset=dataset)
merged_time<-rbind(merged_time,scanpy_data)
saveRDS(merged_time,'../Results/all_runtime_results.rds')

merged_time<-readRDS('../Results/all_runtime_results.rds')

library(writexl)
write_xlsx(merged_time,paste0('../Results/runtime_results.xlsx'))


two_sets_time<-merged_time[which(merged_time$dataset=='PBMC4k' | merged_time$dataset=='BaronPanc'), ]
rownames(two_sets_time)<-1:nrow(two_sets_time)
time_to_plot<-two_sets_time[-c(5,6,9,10,11,12,15,16,17,18,21,22,27,28,35,36,39,40, 43),]
time_to_plot[,2]<-round(as.numeric(time_to_plot[,2]),digits = 3)
time_to_plot[,3]<-round(as.numeric(time_to_plot[,3]),digits = 3)

library(plyr)
time_to_plot$method<-revalue(time_to_plot$method,c("dbscanmdims" = "dbscan", 
                               "kmeansmdims" ="k-means",                     
                               "Path_metrics_clustering_1.5_mdims"="PM_1.5",
                               "Path_metrics_clustering_2_mdims"="PM_2" ,
                               "Path_metrics_clustering_4_mdims"="PM_4" ,
                               "RaceID_clustering"='RaceID3' ,
                               "SC3_clustering"='SC3',                   
                               "Seurat_clustering_default_res_0.8" ='Seurat_def_res_0.8',
                               "Seurat_Clustering_PM_processing"='Seurat',
                               "Seurat_def_clustering" ='Seurat(def.)',
                               "SIMLR_clustering"='SIMLR' ,                
                               "tsne+kmeansmdims"="t-SNE + k-means",
                               "umap+dbscanmdims"= "UMAP+db",
                               "scanpy" = "Scanpy") ) 
time_to_plot$dataset<-factor(time_to_plot$dataset, levels= c( "PBMC4k","BaronPanc" ))


library(ggplot2)
library(ggthemes)
library(scales)

#Complete Runtime 
cbfly <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p_complete<- ggplot(time_to_plot, aes(x=method, y=complete_rt, fill=dataset))
p_complete<-p_complete + geom_bar(stat = "identity", position = 'dodge')+ theme_classic()
p_complete<-p_complete+ labs(title = "Complete Runtime", x = "method", y = "complete_rt(min)")
p_complete<-p_complete+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=15),
                     axis.text.y = element_text(size=12),legend.position='top', 
                     legend.justification='left',
                     legend.direction='horizontal')
p_complete<-p_complete+scale_fill_manual(values=cbfly)+scale_y_continuous(breaks =seq(0,100,5))
p_complete

pdf('../Figures/Fig_5_complete_runtime.pdf')
print(p_complete)
dev.off()
              
#Clustering Runtime
p_clustering<- ggplot(time_to_plot, aes(x=method, y=clustering_rt, fill=dataset))
p_clustering<-p_complete + geom_bar(stat = "identity", position = 'dodge')+ theme_classic()
p_clustering<-p_complete+ labs(title = "Clustering Runtime", x = "method", y = "clustering_rt(min)")
p2<-p_clustering+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+scale_fill_manual(values=cbfly)

pdf('../Figures/Fig_5_clustering_runtime.pdf')
print(p1)
dev.off()