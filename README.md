# Weighted Gene Co-expression Network Analysis 
Weighted Gene Co-expression Network Analysis (WGCNA) is a commonly used unsupervised method to cluster genes based on their expression profiles. 
This repository contains one R script (wgcna.R) that performs weighted correlation network analysis on COVID-19 data.

## Approach
### Step 1: Construct a weighted correlation matrix
•	Obtain a gene expression matrix (count matrix) with rows as genes or probe ids and columns as samples. The values can be intensity value in case of microarray data or quantifications (raw counts) in case of RNA-Seq data. The count matrix is used to find correlation between different genes and finally we obtain a correlation matrix. 
•	Similarity matrix is the measure to find correlation between genes. It can be unsigned or signed. It’s better to use signed matrix since it will provide the information whether the gene is upregulated (positive correlation) or downregulated (negative correlation).
•	Based on this information, the negative ones can be clustered in one separate group and can be further studied to find their role.
•	Adjacency matrix can be obtained either by hard thresholding or soft thresholding.
•	Hard thresholding results in the loss of expression information. It results in unweighted network because the adjacency matrix indicates whether or not a pair of genes are connected. It does not provide any information about how strong or weak is the connection between the pair of genes.
•	Soft thresholding is preferred because it results in a weighted network. It emphasizes more on stronger associations (larger correlation coefficients) and suppress low correlations (noise); raising the correlation to a power will reduce the noise of the correlations in the adjacency matrix.
•	WGCNA R package has a function called “pickSoftThreshold()” to pick the appropriate power term to generate a scale free network (few hubs; non-uniform distribution).

### Step 2: Identify modules
•	Identify modules by clustering genes using the weighted correlation matrix.
•	To cluster genes into network modules use a network proximity measure and the proximity measure used here is the Topological Overlap Measure (TOM). If 2 genes have high correlation, it will have high TOM and low dissimilarity.
•	Merge modules with similar expression profiles into more meaningful modules. This can be done by calculating pairwise Eigengene correlations.
•	Construct the tree (cluster dendrogram), and add a line at the height of 0.25. This height corresponds to a correlation of 75%. Any modules that have correlation greater than 75% are related and hence can be merged.

### Step 3: Correlate modules with phenotypes/traits
•	To correlate or quantify the association between the expression profile and a particular trait of interest, calculate the correlation of the trait with previously identified module Eigengenes.
•	Once the gene significance or the correlation values and the corresponding p-values for all the modules and the traits are available, create a graphical representation like a heatmap which shows the module-trait relationships. This can be helpful to visualize the significance of these correlations with the traits. The row represents the module eigengene and the column represents the traits. Each cell contains a p-value and a correlation value.
•	In addition to correlating modules with external phenotypes or traits, it might also be of interest to identify driver genes in those modules that correlate with the phenotype or trait of interest.


## Study Design
Data: RNA-Seq data of Peripheral Blood Mononuclear Cells (PBMCs) from a group of 17 COVID-19 subjects and 17 healthy controls (GSE152418).
Total samples = 34
NCBI GEO > GSE152418 > Download the supplementary file GSE152418_p20047_Study1_RawCounts.txt.gz (http)

## Biological Questions that are to be addressed:
•	What are the genes or clusters of genes (modules) significantly associated with COVID-19 individuals?
•	What are the genes that are significantly associated with severe COVID-19 cases?

    
## Requirements:
- R (v4.2.1)
- R libraries and packages:
•	BiocManager
•	GO.db
•	impute
•	devtools
•	WGCNA (v1.71)
•	DESeq2 (v1.32.0)
•	GEOquery (v2.60.0)
•	tidyverse (v1.3.1)
•	CorLevelPlot (v0.99.0)
•	gridExtra (v2.3)
•	ggplot2
•	annotables


## Result
1. Cluster Dendrogram plot: There are fewer colours in the merged than the unmerged which indicates that there were a lot of similar modules in the unmerged which were merged together in one module.
   The merged modules are used for the further analysis.

2. Heatmap: The level of significance is associated with the number of asterisks. Three asterisk means that the module has high significance to that trait of interest.
   So, answering the first part of the question as to which are the genes or the cluster of genes that have significant association with COVID-19 individuals.
   Module grey, turquoise, red, salmon and yellow have high significance and they are significantly associated with the diseased state.
   Similarly, module turquoise is significantly associated with severe disease condition.

3. Associating module eigengenes (continuous variable) with diseases severity (categorical variable):
   The idea here is not to find correlation but to identify association/relation between eigengene modules and a categorical trait of interest and to evaluate whether they are statistically significant or not.
   The other alternatives to correlation method to finding significant association between module eigengenes and categorical traits is by:
    •	performing t-test assuming that the variance is equal.
    •	using linear models
   All these methods will yield a similar result.

4. The highly connected intramodular hub genes in the turquoise module was identified by calculating the correlation of the module eigengenes and the gene expression profiles.
   The module membership measures were obtained by calculating the correlation first and then obtained the p-values.
   Looking at the module membership measures and the module membership measures p-values we can identify which are the genes having high membership measures in the modules of interest.
   Looking at the p-values we can identify which are the genes having significantly high module membership in the modules of interest.
   Hence, answering the second part of the question as to what are the genes that are associated with severe COVID-19 cases.

5. The Ensembl IDs of top 25 genes were converted into Gene IDs.

## Conclusion
Hence, the genes from our module of interest (turquoise module) are identified.
The top 25 genes that are associated with severe COVID-19 cases are:
 1. Cell division cycle 42 (CDC42)                             
 2. Aurora kinase A (AURKA)                                  
 3. Cell division cycle 45 (CDC45)                             
 4. Interleukin 21 receptor (IL21R)                           
 5. Signal peptidase complex subunit 2 (SPCS2)                 
 6. Cyclin dependent kinase inhibitor 2C (CDKN2C)               
 7. PRKR interacting protein 1 (PRKRIP1)                         
 8. Signal peptidase complex subunit 3 (SPCS3)                 
 9. Kinesin light chain 4 (KLC4)                             
10. ASXL transcriptional regulator 2 (ASXL2)                   
11. Essential meiotic structure-specific endonuclease 1 (EME1)
12. MRN complex interacting protein (MRNIP)                   
13. Small nuclear ribonucleoprotein U11/U12 subunit 25 (SNRNP25)
14. Solute carrier family 25 member 46 (SLC25A46)                
15. Cation channel sperm associated 2 (CATSPER2)                  
16. Kelch like family member 26 (KLHL26)                       
17. BUB1 mitotic checkpoint serine/threonine kinase (BUB1)   
18. SMAD family member 1 (SMAD1)                              
19. BCL6 corepressor (BCOR)                                  
20. Nuclear receptor binding protein 2 (NRBP2)                 
21. Endoplasmic reticulum oxidoreductase 1 alpha (ERO1A)      
22. Transcription factor Dp-1 (TFDP1)                         
23. Small nucleolar RNA, C/D box 59A (SNORD59A)                  
24. Leucine rich repeat containing 37 member A2 (LRRC37A2)       
25. TEC                          

These genes can be further looked into to study the underlying mechanisms and pathways that contribute to the severity of the disease COVID-19.

  

   

