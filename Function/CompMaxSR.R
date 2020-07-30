CompMaxSR <-
function(adj.m){

require(igraph);
# find right eigenvector of adjacency matrix
fa <- function(x,extra=NULL) {
   as.vector(adj.m %*% x)
}
ap.o <- arpack(fa, options=list(n=nrow(adj.m),nev=1,which="LM"),sym=TRUE);
#计算特征向量
v <- ap.o$vectors;
lambda <- ap.o$values;
maxSR <- log(lambda); ### maximum entropy
return(maxSR);

}
