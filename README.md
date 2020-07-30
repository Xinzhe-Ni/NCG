# NCG
NCG is a single-cell potency model to predict single-cell differentiation and discriminate pluripotent and non-pluripotent cells.
NCG was proposed by following steps: 
(1) computation of Edge Clustering Coefficient (ECC) based on PPI network, 
(2) computation of gene-to-gene functional similarities using GO terms, 
(3) Assignment of gene-to-gene functional similarities as weights to ECC, 
(4) Combination of weighted ECC with scRNA-Seq data.