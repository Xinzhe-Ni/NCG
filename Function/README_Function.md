This is a description of functions in NCG program. It needs to be emphasized that "DoIntegPPI" and "DoSCENT" are from github: https://github.com/aet21/SCENT and "DoSCENTalt" is an altered version of "DoSCENT". Signaling entropy in SCENT program is replaced by NCG. Here are the details:



### Processdata

**description**: Processes all datasets except Yao1 including quantile normalization, log2-transformation and other steps, and performs a gene ID conversion on PPI network.

**usage**: Processdata(counts, PPI)

**arguments**: 

counts: The scRNA-seq data matrix with rows labeling genes and columns labeling single cells.

PPI: The adjacency matrix of a user-given PPI network with rownames and
colnames labeling a gene ID.

**value**:

exp: The scRNA-seq data matrix after preprocessing.

adj: The adjacency matrix of PPI network after gene ID conversion.



### DoIntegPPI

**description**: This function finds the common genes between the scRNA-Seq data matrix and the genes present in the PPI network, and constructs the maximally connected subnetwork and reduced expression matrix for the computation of signaling entropy.

**usage**: DoIntegPPI(exp.m, ppiA.m)

**arguments**: 

exp.m: The scRNA-Seq data matrix with rows labeling genes and columns labeling single cells.

ppiA.m: The adjacency matrix of a user-given PPI network with rownames and
colnames labeling a gene ID.

**value**:

expMC: Reduced expression matrix with genes in the maximally connected subnetwork.

adjMC: Adjacency matrix of the maximally connected subnetwork.



### CompECC

**description**: For an adjacency matrix of protein-protein interaction network, computes the Edge Clustering Coefficient matrix, whose size is the same as PPI network matrix.

**usage**: CompECC(PPI)

**arguments**: 

PPI: The adjacency matrix of a user-given PPI network with rownames and
colnames labeling a gene ID.

**value**:

ECC:  The Edge Clustering Coefficient matrix computed by PPI network.



### CompNCG

**description**: This function finds the common genes between the scRNA-Seq data matrix, ECC matrix and gene-to-gene functional similarity matrix, and computes the NCG value for each cell.

**usage**: CompNCG(ECC, exp, km)

**arguments**:

ECC: The output value of function CompECC.

exp: The output value of function Processdata.

km: The pre-compiled pairwise Kappa similarity matrix on Gene Ontology of human genes.

**value**: 

NCG_nor: The normalized single-cell potency measure computing by Edge Clustering Coefficient, scRNA-seq and Gene Ontology similarity scores.



### DoSCENT

**description**: Main user function implement SCENT. This function infers the discrete potency states of single cells, its distribution across the single cell
population, identifies potency-coexpression clusters of single cells, called landmarks, and finally infers the dependencies of these landmarks which can aid in recontructing lineage trajectories in time course scRNA-Seq experiments.



### DoSCENTalt

**description**: An altered version of DoSCENT. This function deletes the logarithmic step of signaling entropy to do a Gaussian fitting and keeps other steps unchanged.

**usage**: DoSCENTalt(exp.m, NCG, pheno.v = NULL, mixmod = NULL, maxPS = 5, pctG = 0.01, kmax = 9, pctLM = 0.05, pcorTH = 0.1)

**arguments**:

exp.m: Normalized single-cell RNA-Seq data matrix, with rows labeling genes and columns labeling single cells.

NCG: The output value of function CompNCG.

pheno.v: A phenotype vector for the single cells, of same length and order as the columns of exp.m.

mixmod: Specifies whether the Gaussian mixture model to be fit assumes
components to have different (default) or equal variance. In the latter
case, use mixmod=c("E").

maxPS: Maximum number of potency states to allow, when inferring discrete
potency states of single cells. Default value is 5.

pctG: Percentage of all genes in \code{exp.m} to select from each principal
component in an SVD/PCA of \code{exp.m}. The union set of all selected genes is then used for clustering. Default value is 0.01.

kmax: Maximum number of co-expression clusters to allow when performing
clustering. Default value is 9. Larger values are not allowed.

pctLM: Percentage of total number of single cells to allow as a minimum size
for selecting interesting landmarks i.e. potency-coexpression clusters
of single cells. Default value is 0.05.

pcorTH: Threshold for calling significant partial correlations. Default value is
0.1. Usually, single-cell experiments profile large number of cells, so
0.1 is a sensible threshold.

**value**:

potS: Inferred discrete potency states for each single cell. It is indexed so that the index increases as the NCG of the state decreases.

distPSPH: Table giving the distribution of single-cells across potency states and phenotypes.

prob: Table giving the probabilities of each potency state per phenotype value.

hetPS: The normalised NCG of potency per phenotype value.

cl: The co-expression clustering index for each single cell.

pscl: The potency coexpression clustering label for each single cell.

distPSCL: The distribution of single cell numbers per potency state and coexpression cluster.

medLM: A matrix of medoids of gene expression for the selected landmarks.

srPSCL: The average NCG of single cells in each potency coexpression cluster.

srLM: The average NCG of single cells in each landmark.

distPHLM: Table giving the distribution of single cell numbers per phenotype and landmark.

cellLM: Nearest landmark for each single cell.

cellLM2: A vector specifying the nearest and next-nearest landmark for each single cell.

adj: Weighted adjacency matrix between landmarks with entries giving the number of single cells mapping closest to the two landmarks.

pcorLM: Partial correlation matrix of landmarks as estimated from the expression medoids.

<<<<<<< HEAD
netLM: Adjacency matrix of landmarks specifying which partial correlations are significant.
=======
netLM: Adjacency matrix of landmarks specifying which partial correlations are significant.
>>>>>>> 952cc1d5a5bc1180d555430e1a7baf7c17709157
