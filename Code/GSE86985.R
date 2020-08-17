##Data processing
source('Processdata.R');
load("hprdAsigH-13Jun12.Rd");
counts <- read.table("GSE86985_trueseq.tpm.csv",sep=",",header=T);
data <- Processdata(counts,hprdAsigH.m);

##Calculation of NCG
source('DoIntegPPI.R');
source('CompECC.R');
source('CompNCG.R');
load("hs_km.Rda");
int.o <- DoIntegPPI(data$exp,data$adj);
ECC <- CompECC(int.o$adjMC);
NCG <- CompNCG(ECC,int.o$expMC,km);

##Construction of lineage trajectory
source('DoSCENTalt.R');
pheno.v <- c(rep(1,5),rep(2,3),rep(3,4),rep(4,6),rep(5,5),rep(6,6),rep(7,4));
#day0,day6,day9,day12,day19,day26,day54
scent.o_NCG <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred1 <- prediction(c(NCG[c(1:7,18:22)],NCG[c(13:15,29:31)]),c(rep(1,12),rep(0,6)));
perf1 <- performance(pred1,"tpr","fpr");
auc1 <- performance(pred1,'auc')@y.values;
pred2 <- prediction(c(NCG[c(1:7,18:22)],NCG[c(16,17,32,33)]),c(rep(1,12),rep(0,4)));
perf2 <- performance(pred2,"tpr","fpr");
auc2 <- performance(pred2,'auc')@y.values;
pred3 <- prediction(c(NCG[c(13:15,29:31)],NCG[c(16,17,32,33)]),c(rep(1,6),rep(0,4)));
perf3 <- performance(pred3,"tpr","fpr");
auc3 <- performance(pred3,'auc')@y.values;

##P-value evaluation
PI1 <- wilcox.test(NCG[c(1:7,18:22)],NCG[c(13:15,29:31)],alternative = "greater");
PI2 <- wilcox.test(NCG[c(1:7,18:22)],NCG[c(16,17,32,33)],alternative = "greater");
PI3 <- wilcox.test(NCG[c(13:15,29:31)],NCG[c(16,17,32,33)],alternative = "greater")
