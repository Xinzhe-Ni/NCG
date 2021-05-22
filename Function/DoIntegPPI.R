##Teschendorff A E, Enver T. Single-cell entropy for accurate estimation of differentiation potency from a cell’s transcriptome[J]. Nature communications, 2017, 8(1): 1-15.

##############################################################################
##
##  description: This function finds the common genes between the scRNA-Seq data matrix and the genes present in the PPI network, and constructs the maximally connected subnetwork and reduced expression matrix for the computation of signaling entropy.

##  usage: DoIntegPPI(exp.m, ppiA.m)

##  arguments: 

##  exp.m: The scRNA-Seq data matrix with rows labeling genes and columns labeling single cells.
##  ppiA.m: The adjacency matrix of a user-given PPI network with rownames and colnames labeling a gene ID.

##  value:

##  expMC: Reduced expression matrix with genes in the maximally connected subnetwork.
##  adjMC: Adjacency matrix of the maximally connected subnetwork.
##
##############################################################################

DoIntegPPI <-
function(exp.m,ppiA.m){
    
require(igraph);

commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
match(commonEID.v,rownames(exp.m)) -> map1.idx;
expPIN.m <- exp.m[map1.idx,];

match(commonEID.v,rownames(ppiA.m)) -> map2.idx;
adj.m <- ppiA.m[map2.idx,map2.idx];

gr.o <- graph.adjacency(adj.m,mode="undirected");
comp.l <- clusters(gr.o);
cd.v <- summary(factor(comp.l$member));
mcID <- as.numeric(names(cd.v)[which.max(cd.v)]);
maxc.idx <- which(comp.l$member==mcID);
adjMC.m <- adj.m[maxc.idx,maxc.idx];
expMC.m <- expPIN.m[maxc.idx,];

return(list(expMC=expMC.m,adjMC=adjMC.m));
}
