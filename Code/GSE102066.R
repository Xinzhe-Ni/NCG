##Data processing
source('D:/南开/国创/终极版/需要的函数/Processdata.R');
counts<-read.table("D:/南开/国创/需处理的数据/GSE102066/GSE102066_normalized_counts.experiment1.DataMatrix.txt",sep="\t",header=T)
load("D:/南开/国创/代码/scent/hprdAsigH-13Jun12.Rd");
data<-Processdata(counts,hprdAsigH.m);

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
save(SR.v,file="D:/南开/国创/数据/SR.v/GSE102066.SR.Rdata");
#load("D:/南开/国创/数据/SR.v/GSE102066.SR.Rdata");

##Calculation of NCG
source('D:/南开/国创/终极版/需要的函数/CompECC.R');
source('D:/南开/国创/终极版/需要的函数/CompNCG.R');
ECC<-CompECC(int.o$adjMC);
save(ECC,file="D:/南开/国创/数据/ECC.m/GSE102066_ECC.RData");
#load("D:/南开/国创/数据/ECC.m/GSE102066_ECC.RData");
load("D:/南开/国创/代码/slice/data/hs_km.Rda");
NCG<-CompNCG(ECC,int.o$expMC,km);
save(NCG,file="D:/南开/国创/数据/NCG.v/GSE102066.RData");
#load("D:/南开/国创/数据/NCG.v/GSE102066_NCG.RData");

##Construction of lineage trajectory
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoSCENT.R');
source('D:/南开/国创/终极版/需要的函数/DoSCENTalt.R');
pheno.v<-c(rep(1,80),rep(2,78),rep(3,85),rep(4,80),rep(5,79),rep(6,81));
#0d,1d,5d,7d,10d,30d
scent.o <- DoSCENT(data$exp,sr.v=SR.v,pheno.v);
scent.o1 <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred<-prediction(c(NCG[1:158],NCG[403:483]),c(rep(1,158),rep(0,81)));
perf<-performance(pred,"tpr","fpr");
auc<-performance(pred,'auc')@y.values

##P-value evaluation
PI<-wilcox.test(perNCG[1:158],perNCG[403:483],alternative = "greater")
View(PI)
