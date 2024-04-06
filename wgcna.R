# 27-03-2024
# Script to perform weighted gene co-expression network analysis (WGCNA)
# setwd("C:/Users/salon/OneDrive/Desktop/WGCNA")

# Load the libraries
library(BiocManager)
BiocManager::install("GO.db", install_args = "-t 1800")
library(GO.db)
install.packages("WGCNA", install_args = "-t 900", dependencies = T)
BiocManager::install("impute")
library(impute)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(devtools)
devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
library(gridExtra)
library(ggplot2)

allowWGCNAThreads()   # Allow multi-threading (optional)

# Step 1: Fetch Data ---------------------------------------------
data <- read.delim("GSE152418_p20047_Study1_RawCounts.txt", header = T)

# Get metadata
gse <- getGEO(GEO = 'GSE152418', GSEMatrix = T)
metadata <- pData(phenoData(gse[[1]]))
metadata <- metadata[,c(1,2,46:50)]

# Prepare data
data <- data %>% 
  gather(key = 'sample', value = 'counts', -ENSEMBLID) %>%
  mutate(sample = gsub('\\.', '-', sample)) %>%
  inner_join(., metadata, by = c('sample' = 'title')) %>%
  select(1,3,4) %>%
  spread(key = 'geo_accession', value ='counts') %>% 
  column_to_rownames(var = 'ENSEMBLID')

# Step 2: Quality Control - outlier detection -----------------------------
# Detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# Remove the detected outlier genes
data <- data[gsg$goodGenes == T,]

# Detect outlier sample - Hierarchial clustering - method 1
htree <- hclust(dist(t(data)), method = 'average')
plot(htree)

# Detect outlier sample - PCA - method 2
pca <- prcomp(t(data))
pca_data <- data.frame(pca$x)
pca_var <- pca$sdev^2
pca_var_percentage <- round((pca_var)/(sum(pca_var))*100, digits = 2)
ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label=rownames(pca_data))+
  labs(title = "GEO IDS of RNA-Seq Data of Peripheral Blood Mononuclear Cells (PBMCs) from 34 Individuals",
       x = paste0('PC1: ', pca_var_percentage[1], ' %'),
       y = paste0('PC2: ', pca_var_percentage[2], ' %')) +
  theme(plot.title = element_text(size=10, face="bold.italic"))

## If there are any batch effects observed, correct them before moving ahead

# Exclude the outlier samples
samples_to_be_excluded <- c('GSM4614993', 'GSM4614995', 'GSM4615000')
counts_data <- data[, !(colnames(data) %in% samples_to_be_excluded)]


# Step 3: Normalization ------------------------------------------
# Create a DeSeq2 dataset
colData <- metadata %>% 
  filter(!row.names(.) %in% samples_to_be_excluded) 

# Fixing column names in column data
names(colData)
names(colData) <- gsub(':ch1', '', names(colData)) 
names(colData) <- gsub(' ', '_', names(colData))

# making sure the row names in colData matches to column names in counts_data
all(row.names(colData) %in% colnames(counts_data))

# are they in the same order?
all(row.names(colData) == colnames(counts_data))

# Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,    # Not specifying model because
                              colData = colData,          # we need deseq dataset to perform
                              design = ~ 1)               # variance stabilizing transformation

# Remove all genes with counts < 15 in more than 75% of the samples (31*0.75 = 23.25)
# Suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]

# Perform variance stabilizing transformation
dds75_vst <- vst(dds75)

# Get normalized counts
norm_counts <- assay(dds75_vst) %>% 
  t()   # Transpose the norm_counts since we want rows as the samples and columns as the genes

# Step 4: Network Construction --------------------------------------
# Choose a set of soft thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm_counts, powerVector = power, 
                         networkType = 'signed', verbose = 3)
sft_data <- sft$fitIndices

# Viualization to pick ideal power (max R square value, min mean connectivity)

a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = 'red')+
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2')+
  theme_classic()

a2 <- ggplot(sft_data, aes(Power, mean.k., label=Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  labs(x='Power', y='Mean Connectivity')+
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Convert norm_counts matrix to numeric
norm_counts[] <- sapply(norm_counts, as.numeric)
soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# Memory estimate w.r.t block size
bwnet <- blockwiseModules(norm_counts, maxBlockSize = 6500, TOMType = 'signed', 
                          power = soft_power, mergeCutHeight = 0.25, #Threshold that we want to merge the similar modules
                          numericLabels = F, randomSeed = 1234, verbose = 3)

cor <- temp_cor

# Step 5: Module Eigengenes ----------------------------------
module_eigengenes <- bwnet$MEs
head(module_eigengenes)

# Get number of genes for each module
table(bwnet$colors)

#Convert labels to colors for plotting
merged_colors <- labels2colors(bwnet$colors)
unmerged_colors <- labels2colors(bwnet$unmergedColors)

# Plot the dendrogram and the module colors before and after merging underneath

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(unmerged_colors[bwnet$blockGenes[[1]]], merged_colors[bwnet$blockGenes[[1]]]), 
                    c("Unmerged", "Merged"), dendroLabels = F, addGuide = T,
                    hang = 0.03, guideHang = 0.05)


bwnet$colors[bwnet$blockGenes[[1]]]

# 6A. Relate modules to traits --------------------------------------------------

# Create traits file - binarize categorical variables
traits <- colData %>%
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)

# Binarize categorical variables
colData$severity <- factor(colData$severity, levels = c('Healthy', 'Convalescent', 
                                                        'Moderate', 'ICU', 'Severe'))
severity_output <- binarizeCategoricalColumns(colData$severity, minCount = 1)

# Combine traits data with severity output
traits <- cbind(traits,severity_output)

# Define number of genes and samples
nGenes <- ncol(norm_counts)
nSamples <- nrow(norm_counts)

module_trait_corr <- cor(module_eigengenes, traits, use = 'p')
module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSamples)

# Visualize module-trait-correlation as a heatmap
heatmap_data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap_data)  

heatmap_data <- heatmap_data %>% 
 column_to_rownames(var = 'Row.names') 

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[19:23], # columns with trait names
             y = names(heatmap_data)[1:18],  # columns with eigengen names
             col = c('blue', 'skyblue', 'white', 'pink', 'red'))

# Extract genes from the modules that are significantly associated with the diseased state
# and significantly associated with severe disease condition
module_gene_mapping <- as.data.frame(bwnet$colors)
severe_covid19_genes <- module_gene_mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames() %>% 
  as.data.frame()


# 6B. Intramodular analysis: Identifying driver genes ----------------------------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the 
# correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module_membership_measures <- cor(module_eigengenes, norm_counts, use = 'p')
module_membership_measures_pvals <- corPvalueStudent(module_membership_measures, nSamples)

module_membership_measures_pvals[1:10, 1:10]

# Calculate the gene significance and the associated p-values
gene_signif_corr <- cor(norm_counts, traits$data.Severe.vs.all, use ='p') 
gene_signif_corr_pvals <- corPvalueStudent(gene_signif_corr, nSamples)

top_25 <- gene_signif_corr_pvals %>% 
  as.data.frame() %>% 
  rename(., p_values = V1) %>% 
  arrange(p_values) %>% 
  head(25) 


# Convert ensembl_ids to gene ids 
library(annotables)
severe_covid19_genes <- grch38 %>% 
  filter(ensgene %in% row.names(top_25))
 



