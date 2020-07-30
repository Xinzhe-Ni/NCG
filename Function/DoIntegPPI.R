DoIntegPPI <-
function(exp.m,ppiA.m){
    
require(igraph);

commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
#取交集
match(commonEID.v,rownames(exp.m)) -> map1.idx;
expPIN.m <- exp.m[map1.idx,];

match(commonEID.v,rownames(ppiA.m)) -> map2.idx;
adj.m <- ppiA.m[map2.idx,map2.idx];

gr.o <- graph.adjacency(adj.m,mode="undirected");#作图
comp.l <- clusters(gr.o)#计算最大连接
cd.v <- summary(factor(comp.l$member));
mcID <- as.numeric(names(cd.v)[which.max(cd.v)]);
maxc.idx <- which(comp.l$member==mcID);
adjMC.m <- adj.m[maxc.idx,maxc.idx];
expMC.m <- expPIN.m[maxc.idx,];

return(list(expMC=expMC.m,adjMC=adjMC.m));
}
