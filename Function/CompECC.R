##############################################################################
##
##  **description**: For an adjacency matrix of protein-protein interaction network, computes the Edge Clustering Coefficient matrix, whose size is the same as PPI network matrix.

##  **usage**: CompECC(PPI)

##  **arguments**: 

##  PPI: The adjacency matrix of a user-given PPI network with rownames and
##  colnames labeling a gene ID.

##  **value**:

##  ECC:  The Edge Clustering Coefficient matrix computed by PPI network.
##
##############################################################################

CompECC <-
  function(PPI){
    
    N <- vector();
    for(i in 1:ncol(PPI)){
      N[i] <- sum(PPI[i,]);
    }
    ECC <- matrix(0,nrow=nrow(PPI),ncol=ncol(PPI));
    Conode <- PPI%*%PPI;
    for(i in 1:nrow(PPI)){
      for(j in 1:nrow(PPI)){
        if(PPI[i,j]==0||N[i]==1||N[j]==1){
          ECC[i,j] <- 0;
        }
        else{
          ECC[i,j] <- Conode[i,j]/min(N[i]-1,N[j]-1);
        }
      }
    }
    rownames(ECC) <- rownames(PPI);
    colnames(ECC) <- colnames(PPI);
    
    return(ECC);
  }
