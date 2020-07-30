CompECC<-
  function(PPI){
    
    N<-vector();
    for(i in 1:ncol(PPI)){N[i]<-sum(PPI[i,]);};
    ECC<-matrix(0,nrow=nrow(PPI),ncol = ncol(PPI));
    Conode<-PPI%*%PPI;
    for(i in 1:nrow(PPI)){
      for(j in 1:nrow(PPI)){
        if(PPI[i,j]==0||N[i]==1||N[j]==1){
          ECC[i,j]=0;
        }
        else{
          ECC[i,j]<-Conode[i,j]/min(N[i]-1,N[j]-1);}
      }
    }
    rownames(ECC)<-rownames(PPI);
    colnames(ECC)<-colnames(PPI);
    
    return(ECC);
  }