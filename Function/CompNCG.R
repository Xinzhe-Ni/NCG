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
    
    ECCGO<-ECCint*kmint
    NC<- vector()
    for(i in 1:ncol(ECCGO)){
      NC[i]<-sum(ECCGO[i,])}
    NCG<-as.vector(t(expint)%*%NC)
    NCG_nor<-NCG/max(NCG)
    
    return(NCG_nor)
  }