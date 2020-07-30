CompNS <-
function(p.v){
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }
  else { ### one degree nodes have zero entropy, avoid singularity.
    S <- 0;
  }
  return(S);
}
