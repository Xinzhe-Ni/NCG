CompS <-
function(p.v){
  tmp.idx <- which(p.v>0);
  S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(S);
}
