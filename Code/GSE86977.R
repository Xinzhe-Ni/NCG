##Data processing
require("preprocessCore");
require("AnnotationDbi");
require("org.Hs.eg.db");
load("hprdAsigH-13Jun12.Rd");
counts <- read.csv("GSE86977_UMI_20K.2684.csv",header=T);
counts <- counts[!duplicated(counts[,1]),];
row.names(counts) <- counts[,1];
counts <- counts[,-1];
counts[counts > 16] <- 16;
ncounts <- normalize.quantiles(as.matrix(counts),copy=FALSE);
ncounts[ncounts < 1] <- 0.1;
#Gene ID Conversion
r <- rownames(hprdAsigH.m);
geneIDselect <- select(org.Hs.eg.db, keys=r,columns="SYMBOL", keytype="ENTREZID");
rownames(hprdAsigH.m) <- geneIDselect[,2];
colnames(hprdAsigH.m) <- geneIDselect[,2];
data <- list(exp=ncounts,adj=hprdAsigH.m);

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
pheno.v <- c(rep(1,40),rep(2,504),rep(3,278),rep(4,595),rep(5,502),rep(6,765));
#day0,day12,day19,day26,day40,day54
scent.o_NCG <- DoSCENTalt(data$exp,sr.v=NCG,pheno.v);

##AUC evaluation
require("ROCR");
pred1 <- prediction(c(NCG[1:40],NCG[823:1417]),c(rep(1,40),rep(0,595)));
perf1 <- performance(pred1,"tpr","fpr");
auc1 <- performance(pred1,'auc')@y.values;
pred2 <- prediction(c(NCG[1:40],NCG[1920:2684]),c(rep(1,40),rep(0,765)));
perf2 <- performance(pred2,"tpr","fpr");
auc2 <- performance(pred2,'auc')@y.values;
pred3 <- prediction(c(NCG[823:1417],NCG[1920:2684]),c(rep(1,595),rep(0,765)));
perf3 <- performance(pred3,"tpr","fpr");
auc3 <- performance(pred3,'auc')@y.values;

##P-value evaluation
PI1 <- wilcox.test(NCG[1:40],NCG[823:1417],alternative = "greater");
PI2 <- wilcox.test(NCG[1:40],NCG[1920:2684],alternative = "greater");
PI3 <- wilcox.test(NCG[823:1417],NCG[1920:2684],alternative = "greater")
