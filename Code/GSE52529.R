##Data processing
source('Processdata.R');
load("hprdAsigH-13Jun12.Rd");
counts <- read.table("GSE52529_fpkm_matrix.txt",sep="\t",header=T);
data <- Processdata(counts,hprdAsigH.m);
#Gene ID Conversion
r <- as.character(rownames(data$exp));
geneIDselect <- mapIds(org.Hs.eg.db,keys=r,column="SYMBOL",keytype="ENSEMBL",multiVals="first");
rownames(data$exp) <- geneIDselect;
data$exp <- data$exp[!is.na(rownames(data$exp)),];

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
pheno.v <- c(rep(1,96),rep(2,96),rep(3,96),rep(4,84));
#0h,24h,48h,72h
scent.o_NCG <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred <- prediction(c(NCG[1:96],NCG[289:372]),c(rep(1,96),rep(0,84)));
perf <- performance(pred,"tpr","fpr");
auc <- performance(pred,'auc')@y.values;

##P-value evaluation
PI <- wilcox.test(NCG[1:96],NCG[289:372],alternative = "greater")
