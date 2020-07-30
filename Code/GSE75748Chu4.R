##Data processing
source('D:/南开/国创/终极版/需要的函数/Processdata.R');
counts<-read.table("D:/南开/国创/需处理的数据/GSE75748Chu4/GSE75748_bulk_time_course_ec.csv",sep=",",header=T)
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
save(SR.v,file="D:/南开/国创/数据/SR.v/GSE75748Chu4.SR.Rdata");
#load("D:/南开/国创/数据/SR.v/GSE75748Chu4.SR.Rdata");

##Calculation of NCG
source('D:/南开/国创/终极版/需要的函数/CompECC.R');
source('D:/南开/国创/终极版/需要的函数/CompNCG.R');
source('D:/南开/国创/终极版/需要的函数/CompNC.R');
ECC<-CompECC(int.o$adjMC);
save(ECC,file="D:/南开/国创/数据/ECC.m/GSE75748Chu4_ECC.RData");
#load("D:/南开/国创/数据/ECC.m/GSE75748Chu4_ECC.RData");
load("D:/南开/国创/代码/slice/data/hs_km.Rda");
NCG<-CompNCG(ECC,int.o$expMC,km);
NC<-CompNC(ECC,int.o$expMC,km);
save(NC,file="D:/南开/国创/数据/NCG.v/GSE75748Chu4_NC.RData");
save(NCG,file="D:/南开/国创/数据/NCG.v/GSE75748Chu4_NCG.RData");
#load("D:/南开/国创/数据/NCG.v/GSE75748Chu4_NCG.RData");

##Construction of lineage trajectory
source('D:/南开/国创/代码/scent/scent_1.0/scent/R/DoSCENT.R');
source('D:/南开/国创/终极版/需要的函数/DoSCENTalt.R');
pheno.v<-c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3));
#12h,24h,36h,72h,96h
scent.o <- DoSCENT(data$exp,sr.v=SR.v,pheno.v);
scent.o1 <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred<-prediction(c(SR.v[1:6],SR.v[10:15]),c(rep(1,6),rep(0,6)));
perf<-performance(pred,"tpr","fpr");
auc<-performance(pred,'auc')@y.values

##P-value evaluation
PI<-wilcox.test(SR.v[1:6],SR.v[10:15],alternative = "greater")
View(PI)
