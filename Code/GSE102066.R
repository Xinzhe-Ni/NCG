##Data processing
source('Processdata.R');
load("hprdAsigH-13Jun12.Rd");
counts <- read.table("GSE102066_normalized_counts.experiment.DataMatrix.txt",sep="\t",header=T);
data <- Processdata(counts,hprdAsigH.m);

##Calculation of SR
source('CompMaxSR.R');
source('CompNS.R');
source('CompS.R');
source('CompSRana.R');
source('DoIntegPPI.R');
source('DoSCENT.R');
int.o <- DoIntegPPI(data$exp,data$adj);
maxSR <- CompMaxSR(int.o$adjMC);
SR.v <- vector(); 
for(s in 1:ncol(data$exp)){
  SR.v[s] <- CompSRana(int.o$expMC[,s],int.o$adjMC,maxSR=maxSR)$sr;
}

##Calculation of NCG
source('CompECC.R');
source('CompNCG.R');
load("hs_km.Rda");
ECC <- CompECC(int.o$adjMC);
NCG <- CompNCG(ECC,int.o$expMC,km);

##Construction of lineage trajectory
source('DoSCENT.R');
source('DoSCENTalt.R');
pheno.v <- c(rep(1,80),rep(2,78),rep(3,85),rep(4,80),rep(5,79),rep(6,81));
#day0,day1,day5,day7,day10,day30
scent.o_SR <- DoSCENT(data$exp,sr.v=SR.v,pheno.v);
scent.o_NCG <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred <- prediction(c(NCG[1:158],NCG[403:483]),c(rep(1,158),rep(0,81)));
perf <- performance(pred,"tpr","fpr");
auc <- performance(pred,'auc')@y.values;

##P-value evaluation
PI <- wilcox.test(NCG[1:158],NCG[403:483],alternative = "greater")
