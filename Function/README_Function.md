This is a description of functions in NCG program. It needs to be emphasized that CompMaxSR, CompNS, CompS, CompSRana, DoIntegPPI and DoSCENT are from github: https://github.com/aet21/SCENT. Signaling entropy in SCENT program is replaced by NCG. Here are the details:



CompECC: For an adjacency matrix of protein-protein interaction network, computes the Edge Clustering Coefficient matrix, whose size is the same as PPI network matrix.



CompMaxSR: For a given maximally connected unweighted network, specified by an adjacency matrix, there is a maximum possible signaling entropy rate, which this function computes.



CompNC: This function finds the common genes between the scRNA-Seq data matrix, ECC matrix and gene-to-gene functional similarity matrix, and computes the NC value for each cell.



CompNCG: This function finds the common genes between the scRNA-Seq data matrix, ECC matrix and gene-to-gene functional similarity matrix, and computes the NCG value for each cell.



CompNS: Computes the normalized local signaling entropy for a gene. This is an internal function which the user does not need to invoke.



CompS: Computes the local signaling entropy for a gene. This is an internal function which the user does not need to invoke.



CompSRana: This is the main user function for computing signaling entropy of single cells. It takes as input the gene expression profile of a single cell
and the adjacency matrix of a connected network. These inputs will be
typically the output of the \code{DoIntegPPI} function.



DoIntegPPI: This function finds the common genes between the scRNA-Seq data matrix and the genes present in the PPI network, and constructs the maximally connected subnetwork and reduced expression matrix for the computation of signaling entropy.



DoSCENT: Main user function implement SCENT. This function infers the discrete potency states of single cells, its distribution across the single cell
population, identifies potency-coexpression clusters of single cells, called landmarks, and finally infers the dependencies of these landmarks which can aid in recontructing lineage trajectories in time course scRNA-Seq experiments.



DoSCENTalt: An altered version of DoSCENT. This function deletes the logarithmic step of signaling entropy to do a Gaussian fitting and keeps other steps unchanged.



Processdata: Processes all datasets except Yao1 including quantile normalization, log2-transformation and other steps, and performs a gene ID conversion on PPI network.
