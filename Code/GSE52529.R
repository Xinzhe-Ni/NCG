##Data processing
source('D:/南开/国创/终极版/需要的函数/Processdata.R');
counts<-read.table("D:/南开/国创/需处理的数据/GSE52529/GSE52529_fpkm_matrix.txt",sep="\t",header=T);
load("D:/南开/国创/代码/scent/hprdAsigH-13Jun12.Rd");
data<-Processdata(counts,hprdAsigH.m);
r<-as.character(rownames(data$exp));
geneIDselect<-mapIds(org.Hs.eg.db,keys=r,column="SYMBOL",keytype="ENSEMBL",multiVals="first");
rownames(data$exp)<-geneIDselect;
data$exp=data$exp[!is.na(rownames(data$exp)),];

##Calculation of SR
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/CompMaxSR.R');
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/CompNS.R');
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/CompS.R');
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/CompSRana.R');
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoIntegPPI.R');
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoSCENT.R');
int.o <- DoIntegPPI(data$exp,data$adj);
maxSR <- CompMaxSR(int.o$adjMC);
SR.v <- vector(); 
for(s in 1:ncol(data$exp)){ SR.v[s] <-CompSRana(int.o$expMC[,s],int.o$adjMC,maxSR=maxSR)$sr;}
save(SR.v,file="D:/南开/国创/数据/SR.v/GSE52529.SR.Rdata");
#load("D:/南开/国创/数据/SR.v/GSE52529.SR.Rdata");

##Calculation of NCG
source('D:/南开/国创/终极版/需要的函数/CompECC.R');
source('D:/南开/国创/终极版/需要的函数/CompNCG.R');
ECC<-CompECC(int.o$adjMC);
save(ECC,file="D:/南开/国创/数据/ECC.m/GSE52529.RData");
#load("D:/南开/国创/数据/ECC.m/GSE52529.RData");
load("D:/南开/国创/代码/slice/data/hs_km.Rda");
NCG<-CompNCG(ECC,int.o$expMC,km);
save(NCG,file="D:/南开/国创/数据/NCG.v/GSE52529_NCG.RData");
#load("D:/南开/国创/数据/NCG.v/GSE52529_NCG.RData");

##Construction of lineage trajectory
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoSCENT.R');
source('D:/南开/国创/终极版/需要的函数/DoSCENTalt.R');
pheno.v<-c(rep(1,96),rep(2,96),rep(3,96),rep(4,84));
#0h,24h,48h,72h
scent.o <- DoSCENT(data$exp,sr.v=SR.v,pheno.v);
scent.o <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred<-prediction(c(NCG[1:96],NCG[289:372]),c(rep(1,96),rep(0,84)));
perf<-performance(pred,"tpr","fpr");
auc<-performance(pred,'auc')@y.values

##P-value evaluation
PI<-wilcox.test(NCG[1:96],NCG[289:372],alternative = "greater")
View(PI)
