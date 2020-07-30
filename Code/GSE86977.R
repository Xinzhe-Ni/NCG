##Data processing
require("preprocessCore")
counts<-read.csv("D:/南开/国创/需处理的数据/GSE86977/GSE86977_UMI_20K.2684.csv",header=T)
counts=counts[!duplicated(counts[,1]),]
row.names(counts)<-counts[,1]
counts<-counts[,-1]
counts[counts>16]<-16
counts1<-as.matrix(counts)
ncounts<-normalize.quantiles(counts1,copy=FALSE)
ncounts[ncounts<1]<-0.1
require("AnnotationDbi")
require("org.Hs.eg.db")
require("igraph")
load("D:/南开/国创/代码/scent/hprdAsigH-13Jun12.Rd");
r<-rownames(hprdAsigH.m);
geneIDselect <-select(org.Hs.eg.db, keys=r,columns="SYMBOL", keytype="ENTREZID");
rownames(hprdAsigH.m)<-geneIDselect[,2];
colnames(hprdAsigH.m)<-geneIDselect[,2];
data<-list(exp=ncounts,adj=hprdAsigH.m)

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
save(SR.v,file="D:/南开/国创/数据/SR.v/GSE86977.SR.Rdata");
#load("D:/南开/国创/数据/SR.v/GSE86977.SR.Rdata");

##Calculation of NCG
source('D:/南开/国创/终极版/需要的函数/CompECC.R');
source('D:/南开/国创/终极版/需要的函数/CompNCG.R');
ECC<-CompECC(int.o$adjMC);
save(ECC,file="D:/南开/国创/数据/ECC.m/GSE86977_ECC.RData");
#load("D:/南开/国创/数据/ECC.m/GSE86977_ECC.RData");
load("D:/南开/国创/代码/slice/data/hs_km.Rda");
NCG<-CompNCG(ECC,int.o$expMC,km);
save(NCG,file="D:/南开/国创/数据/NCG.v/GSE86977_NCG.RData");
#load("D:/南开/国创/数据/NCG.v/GSE86977_NCG.RData");

##Construction of lineage trajectory
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoSCENT.R');
source('D:/南开/国创/终极版/需要的函数/DoSCENTalt.R');
pheno.v<-c(rep(1,40),rep(2,504),rep(3,278),rep(4,595),rep(5,502),rep(6,765));
#day0,day12,day19,day26,day40,day54
scent.o <- DoSCENT(data$exp,sr.v=SR.v,pheno.v);
scent.o <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR")
pred1<-prediction(c(NCG[1:40],NCG[823:1417]),c(rep(1,40),rep(0,595)))
perf1<-performance(pred1,"tpr","fpr")
auc1<-performance(pred1,'auc')@y.values
pred2<-prediction(c(NCG[1:40],NCG[1920:2684]),c(rep(1,40),rep(0,765)))
perf2<-performance(pred2,"tpr","fpr")
auc2<-performance(pred2,'auc')@y.values
pred3<-prediction(c(NCG[823:1417],NCG[1920:2684]),c(rep(1,595),rep(0,765)))
perf3<-performance(pred3,"tpr","fpr")
auc3<-performance(pred3,'auc')@y.values

##P-value evaluation
PI1<-wilcox.test(SR.v[1:40],SR.v[823:1417],alternative = "greater")
View(PI1)
PI2<-wilcox.test(SR.v[1:40],SR.v[1920:2684],alternative = "greater")
View(PI2)
PI3<-wilcox.test(SR.v[823:1417],SR.v[1920:2684],alternative = "greater")
View(PI3)
PI<-wilcox.test(NCG[823:1417],NCG[1920:2684],alternative = "greater")
View(PI)
