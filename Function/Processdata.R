Processdata<-
  function(counts,PPI){
    
    require("preprocessCore");
    require("AnnotationDbi");
    require("org.Hs.eg.db");
    
    counts <- counts[!duplicated(counts[,1]),];
    rownames(counts) <- counts[,1];
    counts <- counts[,-1];
    ncounts <- normalize.quantiles(as.matrix(counts),copy=FALSE);
    ncounts[ncounts < 1] <- 1;
    ncounts <- log2(ncounts+0.1);
    rownames(ncounts) <- rownames(counts);
    colnames(ncounts) <- colnames(counts);
    
    #Gene ID Conversion
    r <- rownames(PPI);
    geneIDselect <- select(org.Hs.eg.db, keys=r,columns="SYMBOL", keytype="ENTREZID");
    rownames(PPI) <- geneIDselect[,2];
    colnames(PPI) <- geneIDselect[,2];
    
    return(list(exp=ncounts,adj=PPI));
  }
