##Teschendorff A E, Enver T. Single-cell entropy for accurate estimation of differentiation potency from a cell’s transcriptome[J]. Nature communications, 2017, 8(1): 1-15.

DoSCENT <-
function(exp.m,sr.v,pheno.v=NULL,mixmod=NULL,maxPS=5,pctG=0.01,kmax=9,pctLM=0.05,pcorTH=0.1){

    require(mclust);
    require(igraph);
    require(isva);
    require(cluster);
    require(corpcor);
    
    ntop <- floor(pctG*nrow(exp.m)); 

    print("Fit Gaussian Mixture Model to Signaling Entropies");
    logitSR.v <- log2(sr.v/(1-sr.v));
    if(is.null(mixmod)){ ## default assumes different variance for clusters
      mcl.o <- Mclust(logitSR.v,G=1:maxPS); 
    }
    else {
      mcl.o <- Mclust(logitSR.v,G=1:maxPS,modelNames=c("E")); 
    }
    potS.v <- mcl.o$class; 
    nPS <- length(levels(as.factor(potS.v)));
    print(paste("Identified ",nPS," potency states",sep=""));
    names(potS.v) <- paste("PS",1:nPS,sep=""); 
    mu.v <- mcl.o$param$mean;
    sd.v <- sqrt(mcl.o$param$variance$sigmasq); 
    avSRps.v <- (2^mu.v)/(1+2^mu.v);

    savSRps.s <- sort(avSRps.v,decreasing=TRUE,index.return=TRUE);
    spsSid.v <- savSRps.s$ix;
    ordpotS.v <- match(potS.v,spsSid.v); 
    if(!is.null(pheno.v)){
       nPH <- length(levels(as.factor(pheno.v)));
       distPSph.m <- table(pheno.v,ordpotS.v)
       print("Compute Shannon (Heterogeneity) Index for each Phenotype class");
       probPSph.m <- distPSph.m/apply(distPSph.m,1,sum);
       hetPS.v <- vector();
       for(ph in 1:nPH){
           prob.v <- probPSph.m[ph,];
           sel.idx <- which(prob.v >0);
           hetPS.v[ph] <-  - sum(prob.v[sel.idx]*log(prob.v[sel.idx]))/log(nPS);
       }
       names(hetPS.v) <- rownames(probPSph.m);
       print("Done");
     }
     else {
       distPSph.m=NULL; probPSph.m=NULL; hetPS.v=NULL;
     }
  
     ### now cluster cells independently of SR
     ### select genes over which to cluster
     print("Using RMT to estimate number of significant components of variation in scRNA-Seq data");
     tmp.m <- exp.m - rowMeans(exp.m);
     rmt.o <- EstDimRMT(tmp.m);  svd.o <- svd(tmp.m);
     tmpG2.v <- vector();
     print(paste("Number of significant components=",rmt.o$dim,sep=""));
     for(cp in 1:rmt.o$dim){
       tmp.s <- sort(abs(svd.o$u[,cp]),decreasing=TRUE,index.return=TRUE);
       tmpG2.v <- union(tmpG2.v,rownames(exp.m)[tmp.s$ix[1:ntop]]);
     }
     selGcl.v <- tmpG2.v;
     
     ### now perform clustering of all cells over the selected genes
     print("Identifying co-expression clusters");
     map.idx <- match(selGcl.v,rownames(exp.m));
     distP.o <- as.dist( 0.5*(1-cor(exp.m[map.idx,])) );
     asw.v <- vector();
     for(k in 2:kmax){
        pam.o <- pam(distP.o,k,stand=FALSE);
        asw.v[k-1] <- pam.o$silinfo$avg.width;
     }
     k.opt <- which.max(asw.v)+1;
     pam.o <- pam(distP.o,k=k.opt,stand=FALSE);
     clust.idx <- pam.o$cluster;
     print(paste("Inferred ",k.opt," clusters",sep=""));
     psclID.v <- paste("PS",ordpotS.v,"-CL",clust.idx,sep="");
     
     ### identify landmark clusters
     print("Now identifying landmarks (potency co-expression clusters)");
     distPSCL.m <- table(paste("CL",clust.idx,sep=""),paste("PS",ordpotS.v,sep=""));
     sizePSCL.v <- as.vector(distPSCL.m);
     namePSCL.v <- vector();
     ci <- 1;
      for(ps in 1:nPS){
       for(cl in 1:k.opt){
          namePSCL.v[ci] <- paste("PS",ps,"-CL",cl,sep="");
          ci <- ci+1;
       }
      }
      names(sizePSCL.v) <- namePSCL.v;
      ldmkCL.idx <- which(sizePSCL.v > pctLM*ncol(exp.m));
      print(paste("Identified ",length(ldmkCL.idx)," Landmarks",sep=""));
    
      ### distribution of phenotypes among LMs
      if(!is.null(pheno.v)){
        tab.m <- table(pheno.v,psclID.v);
        tmp.idx <- match(names(sizePSCL.v)[ldmkCL.idx],colnames(tab.m));
        distPHLM.m <- tab.m[,tmp.idx];
      }
      else {
        distPHLM.m <- NULL;
      }
    
      ### medoids
      print("Constructing expression medoids of landmarks");
      med.m <- matrix(0,nrow=length(selGcl.v),ncol=nPS*k.opt);
      srPSCL.v <- vector();
      ci <- 1;
      for(ps in 1:nPS){
       for(cl in 1:k.opt){
          tmpS.idx <- intersect(which(ordpotS.v==ps),which(clust.idx==cl));
          m<-matrix(exp.m[map.idx,tmpS.idx]);
          e<-unlist(m);
          med.m[,ci] <- apply(matrix(e,nrow=length(map.idx)),1,median);
          srPSCL.v[ci] <- mean(sr.v[tmpS.idx]);
          ci <- ci+1;
        }
       }
       names(srPSCL.v) <- namePSCL.v;
       srLM.v <- srPSCL.v[ldmkCL.idx];
       medLM.m <- med.m[,ldmkCL.idx];
       colnames(medLM.m) <- namePSCL.v[ldmkCL.idx];
       rownames(medLM.m) <- selGcl.v;
       
       ### now project each cell onto two nearest landmarks
       print("Inferring dependencies/trajectories/transitions between landmarks");
       cellLM2.v <- vector(); cellLM.v <- vector();
       for(c in 1:ncol(exp.m)){
        distCellLM.v <- 0.5*(1-as.vector(cor(exp.m[map.idx,c],medLM.m)));
        tmp.s <- sort(distCellLM.v,decreasing=FALSE,index.return=TRUE);
        cellLM2.v[c] <- paste("LM",tmp.s$ix[1],"-LM",tmp.s$ix[2],sep="");
        cellLM.v[c] <- colnames(medLM.m)[tmp.s$ix[1]];
       }
    
       adjLM.m <- matrix(0,nrow=ncol(medLM.m),ncol=ncol(medLM.m));
       rownames(adjLM.m) <- colnames(medLM.m);
       colnames(adjLM.m) <- colnames(medLM.m);
       for(lm1 in 1:ncol(medLM.m)){
        for(lm2 in 1:ncol(medLM.m)){
         adjLM.m[lm1,lm2] <- length(which(cellLM2.v==paste("LM",lm1,"-LM",lm2,sep="")));
        }
       }
       sadjLM.m <- adjLM.m + t(adjLM.m);
       corLM.m <- cor(medLM.m);
       pcorLM.m <- cor2pcor(corLM.m);
       rownames(pcorLM.m) <- rownames(corLM.m);
       colnames(pcorLM.m) <- rownames(corLM.m);    
       netLM.m <- pcorLM.m;diag(netLM.m) <- 0;
       netLM.m[pcorLM.m < pcorTH] <- 0;
       netLM.m[pcorLM.m > pcorTH] <- 1;    
    
       return(list(potS=ordpotS.v,distPSPH=distPSph.m,prob=probPSph.m,hetPS=hetPS.v,cl=clust.idx,pscl=psclID.v,distPSCL=distPSCL.m,medLM=medLM.m,srPSCL=srPSCL.v,srLM=srLM.v,distPHLM=distPHLM.m,cellLM=cellLM.v,cellLM2=cellLM2.v,adj=sadjLM.m,pcorLM=pcorLM.m,netLM=netLM.m));
}
