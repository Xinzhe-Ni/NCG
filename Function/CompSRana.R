CompSRana <-
function(exp.v,adj.m,local=TRUE,maxSR=NULL){

  ### compute outgoing flux around each node
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
   SR <- SR/maxSR;
  }
  if(local){
   NS.v <- apply(p.m,1,CompNS);
  }
  else {
   NS.v <- NULL;
  }
  return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}
