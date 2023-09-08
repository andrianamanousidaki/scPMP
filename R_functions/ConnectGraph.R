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