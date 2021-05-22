##############################################################################
##
##  description: This function finds the common genes between the scRNA-Seq data matrix, ECC matrix and gene-to-gene functional similarity matrix, and computes the NCG value for each cell.

##  usage: CompNCG(ECC, exp, km)

##  arguments: 

##  ECC: The output value of function CompECC.
##  exp: The output value of function Processdata.
##  km: The pre-compiled pairwise Kappa similarity matrix on Gene Ontology of human genes.

##  value:

##  NCG_nor: The normalized single-cell potency measure computing by Edge Clustering Coefficient, scRNA-seq and Gene Ontology similarity scores.
##
##############################################################################

CompNCG<-
  function(ECC,exp,km){
    
    require(igraph);
    
    commonEID1.v <- intersect(rownames(km),rownames(ECC));
    commonEID.v <- intersect(commonEID1.v,rownames(exp));
    match(commonEID.v,rownames(ECC)) -> map1.idx;
    ECCint <- ECC[map1.idx,map1.idx];
    match(commonEID.v,rownames(km)) -> map2.idx;
    kmint <- km[map2.idx,map2.idx];
    match(commonEID.v,rownames(exp)) -> map3.idx;
    expint <- exp[map3.idx,];
    
    EG <- ECCint*kmint;
    ECG <- vector();
    for(i in 1:ncol(EG)){
      ECG[i] <- sum(EG[i,]);
    }
    NCG <- as.vector(t(expint)%*%ECG);
    NCG_nor <- NCG/max(NCG);  #normalization
    
    return(NCG_nor);
  }
