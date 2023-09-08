Path_Metrics_Clustering_R<-function(dataset, # feautures X points
                                    n_clust = NULL,# number of clusters
                                    p=2,# parameter of path metric
                                    k=NULL,# number of nearest neighbors for the kNN graph
                                    SpectralOpts.EmbeddingDimension=40,# Dimensions of PM embedding
                                    SpectralOpts.LearnMDSEmbeddingDimension= TRUE,# Learn the embedding dimesion TRUE or FALSE
                                    run.adjusted.kmeans=TRUE,# Use adjusted kmeans rather than regular kmeans
                                    silhouette=FALSE,# Predict the number of cluster using the silhouetter criterion, TRUE or FALSE
                                    klist=20, # Silhouette will find the best number of clusters among 2:klist
                                    num_landmarks =NULL,
                                    Fast_mds=FALSE # Suggested for data sets with over 2000 points
                                    ){

 require(igraph)
 require(RANN)
 require(Matrix)
 require(expm)
 require(cluster)
 require(cccd)
 require(pracma)
 require(ClusterR)


# Required functions #######################################################################
 
   # Adjusted kmeans clustering ---------------------------------------------------------
  
  adjusted_kmeans<-function(U,k){
    #this function needs library(pracma)
    km<-KMeans_rcpp(U,clusters=k,num_init = 50)
    pm_labels<-km$clusters
    H=table(pm_labels)
    
    if(min(H)< sqrt(nrow(U)/2)){
      num_tiny_clusters=length(which(H<sqrt(nrow(U)/2)))
      km2<- KMeans_rcpp(U,clusters=k+num_tiny_clusters,num_init =50)
      pm_labels_original<-pm_labels
      pm_labels<-km2$clusters
      C<-km2$centroids
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
  
  
 # Connect Graph function -------------------------------------------------------------- 
 
  ConnectGraph<-function(Ap, X, p, random.sample.size.threshold=20){ 
    
    CC=components(graph_from_adjacency_matrix(Ap,weighted = TRUE))
    NumberCC = CC$no
    if(NumberCC>1){
      
      CCindex=list()
      for (i in 1:NumberCC) {
        CCindex[[i]] = which(CC$membership==i)
      }
      
      big_C<-which(CC$csize > random.sample.size.threshold)
      if(length(big_C > 0) ){
        for(i in big_C){
          set.seed(i+53)
          CCindex[[i]]<-sample(CCindex[[i]], size=random.sample.size.threshold, replace = FALSE)
        }
      }   
      
      for (i in 1:(NumberCC-1)) {
        for(j in (i+1):NumberCC){
          I=CCindex[[i]]
          J=CCindex[[j]]
          SampleDistances = as.data.frame(matrix(0,length(I),length(J)))
          for (s in 1:length(I)) {
            for (q in 1:length(J)) {
              #rownames(SampleDistances)[s]<- I[s]
              #colnames(SampleDistances)[q]<-J[q]
              SampleDistances[s,q] = norm(X[I[s],]-X[J[q],],type='2')
            }
          }
          
          MinDistance = min(min(SampleDistances))
          idx = which(SampleDistances == MinDistance, TRUE)
          Ap[I[idx[1,1]],J[idx[1,2]]]=MinDistance^p
          Ap[J[idx[1,2]],I[idx[1,1]]]=Ap[I[idx[1,1]],J[idx[1,2]]]
          
          # s0=as.numeric(rownames(SampleDistances)[idx[1,1]])
          # q0=as.numeric(colnames(SampleDistances)[idx[1,2]])
          # 
          # Ap[s0,q0]=MinDistance^p
          # Ap[q0,s0]=Ap[s0,q0]
        }
      }
    } else {
      Ap=Ap
    }
    return(Ap)
  }
######################################################################################

#Input and Parameters of PM clustering -------------------------------------------------------------------
X=t(dataset)

n = dim(X)[1]

n_clust=n_clust

p=p

if(is.null(k)){k = min(dim(X)[1],500)}else{k=k}


SpectralOpts.EmbeddingDimension = SpectralOpts.EmbeddingDimension #use this as 
#final number of dimensions used for clustering

SpectralOpts.LearnMDSEmbeddingDimension = SpectralOpts.LearnMDSEmbeddingDimension

run.adjusted.kmeans=run.adjusted.kmeans

silhouette = silhouette

klist = klist

# Compute shortest paths --------------------------------------------------

# knn graph              --------------------------------------------------
knn=nn2(data= X,query=X, k = k)
D_knn=as.data.frame(t(knn$nn.dists))
IDX=as.data.frame(t(knn$nn.idx))


Base= matrix(1,n,min(k,dim(X)[1]))
for( i in 1:n){
  Base[i,]=i*Base[i,]
}

Base=t(Base)

zero<- which(D_knn == 0)

A=sparseMatrix(i=unlist(Base)[-zero],j=unlist(IDX)[-zero],x=unlist(D_knn)[-zero], dims = c(n,n))
max_idx<-which(A<t(A),TRUE)
A[max_idx]<-t(A)[max_idx]

Ap = A^p #weight edges by ||x_i-x_j||^p

# If kNN graph is disconnected, add small edges to obtain a connected graph ----

Ap = ConnectGraph(Ap, X, p)
remove(X)
remove(A)

# Compute l_p distances from kNN graph ------------------------------------

PD_Dis = distances( graph_from_adjacency_matrix(Ap,weighted = TRUE),algorithm = "johnson")
PD_Dis = PD_Dis^(1/p)  

# Compute path metric MDS embedding: --------------------------------------
if(Fast_mds){
  # Landmark-based PD-MDS computation:
  print("Fast MDS starts")
  if(is.null(num_landmarks )){
    num_landmarks = ceil(50*log(n))}else{
      num_landmarks =num_landmarks }
  
  m = SpectralOpts.EmbeddingDimension
  landmarks <- sort(sample.int(n, num_landmarks, replace = FALSE))
  is_landmark <- logical(length=n)
  is_landmark[landmarks] <- TRUE
  non_landmarks <- which(is_landmark==FALSE)
  Landmark_Dis <- distances( graph_from_adjacency_matrix(Ap,weighted = TRUE), v = landmarks, algorithm = "johnson")
  Landmark_Dis = Landmark_Dis^(1/p)
  Rs <- Landmark_Dis[ ,landmarks]^2 # Extract small square matrix of squared landmark distances
  R <- t(cbind(Rs, Landmark_Dis[ ,non_landmarks]^2)) 
  ev <- eigen(Rs,symmetric=TRUE)
  o <- order(abs(Re(ev$values)), decreasing = TRUE)
  ev$values <- ev$values[o]
  ev$vectors <- ev$vectors[ ,o]
  V <- ev$vectors[ ,1:ceil(num_landmarks/2)]
  Lam_inv <- diag(1/ev$values[1:ceil(num_landmarks/2)])
  T <- V %*% Lam_inv %*% (t(V))
  QR <- qr( R - ( ones(n,1) %*% (ones(1,n) %*% R) )/n )
  Rqr <- qr.R(QR)
  Q <- qr.Q(QR)
  ev2 <- eigen(-(1/2)* Rqr %*% T %*% t(Rqr),symmetric=TRUE)
  V2 <- ev2$vectors[ ,1:m] 
  Lam2root <- diag(sqrt(ev2$values[1:m]))
  Z <- Q %*% V2 %*% Lam2root
  # Re-order indices to match point order
  U1 <- zeros(n,m)
  U1[landmarks, ] <- Z[1:num_landmarks, ]
  U1[non_landmarks, ] <- Z[(num_landmarks+1):n, ]
  Eigvals <- ev2$values[1:m]
  all_eigenval<-ev2$values
}else{
  print("Classic MDS starts")
  MDS<-cmdscale(PD_Dis, k = 40, eig = TRUE)
  U1<-MDS$points
  Eigvals<-MDS$eig
  Eigvals<-Eigvals[1:40]
  all_eigenval<-MDS$eig
  }

if(SpectralOpts.LearnMDSEmbeddingDimension){
  ratios=Eigvals[1:39]/Eigvals[2:length(Eigvals)]
  max_dim = max(which(Eigvals/max(Eigvals) > .01))-1 #want BOTH eigenvals to not be too small
  ratio_dim_est = which(max(ratios[3:max_dim])==ratios[3:max_dim])+2 #at least 3 dimensions
  #ratio_dim_est = find(max(ratios(2:max_dim))==ratios(2:max_dim))+1 #at least 2 dimensions
  U = U1[,1:ratio_dim_est]
  U1=U
  r=ratio_dim_est
}else{
  U1=U1[,1:SpectralOpts.EmbeddingDimension]
  r=SpectralOpts.EmbeddingDimension
}

# Run clustering: -----------------------------------------------------------

uniqueGroups = n_clust

# Path metric clustering  -----------------------------------

if(silhouette){
  sil_scores<-matrix(0, 1, (klist-1))
  matrix_clustering_labels= matrix(0, n, (klist-1))
  if(run.adjusted.kmeans){
    
    for(k in 2:klist){
      km<-adjusted_kmeans(U1,k)
      matrix_clustering_labels[,(k-1)]=km
      ss <- silhouette(km, dist(U1, method = "euclidean"))
      ss.width<-ss[,3]
      sil_scores[(k-1)]<-mean(ss.width)
    }
    
  }else{
    
    for(k in 2:klist){
      km<-KMeans_rcpp(U1,clusters=k,num_init = 50)
      matrix_clustering_labels[,(k-1)]=km$clusters
      ss <- silhouette(km$clusters, dist(U1, method = "euclidean"))
      ss.width<-ss[,3]
      sil_scores[(k-1)]<-mean(ss.width)}
    
  }
  pm_labels<-matrix_clustering_labels[,which(sil_scores == max(sil_scores))]
}else{
  
  if(run.adjusted.kmeans){
    pm_labels = adjusted_kmeans(U1,n_clust)}else{
      
      pm_labels = KMeans_rcpp(U1,clusters=n_clust,num_init = 50)
      pm_labels<-pm_labels$clusters}
  
}
output_list<-list(pm_labels =pm_labels ,
                  U1=t(U1),
                  r=r,
                  Eigenvalues=all_eigenval)
return(output_list)
}