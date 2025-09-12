# Setup ------

# load libraries
# - AUCell: gene set scoring per cell (AUC-based)
# - MouseGastrulationData: example single-cell datasets from mouse gastrulation
# - SingleR: automated cell type annotation by similarity to a reference
# - bluster: clustering utilities that wrap common clustering algorithms
# - scater: SingleCellExperiment utilities & plotting helpers
# - scran: preprocessing & marker detection helpers
# - pheatmap: simple heatmap visualization for contingency tables
# - GSEABase: representation of gene sets (GeneSet, GeneSetCollection)



# Data retrieval -----

# load WTChimeraData processed dataset (sample 5) and store as a SingleCellExperiment (sce)
# - WTChimeraData(): convenience wrapper from MouseGastrulationData that returns
#   an experiment-like object. Here we request 'processed' counts (already
#   filtered/normalized-ish depending on the package maintainers' preprocessing).
# - We print the object to inspect dimensions and metadata.

# select 1000 cells for a quicker interactive run (set.seed for reproducibility)
# - sampling reduces compute time for this tutorial while preserving structure.


# Pre-processing ----- 
# normalize counts
# - logNormCounts from scater computes size-factor based normalization followed by log-transformation
# - it creates an assay named "logcounts" in the SingleCellExperiment object
# and run PCA
# - runPCA uses the 'logcounts' by default; here we do not pre-select HVGs so it will
#   by default use a set of most variable genes (often implemented inside runPCA)
# - PCA reduces dimensionality for visualization and neighborhood-based operations


# Clustering -----
# simple clustering based on PCA, using clusterCells function
# - BLUSPARAM defines graph construction + clustering algorithm settings.
# - cluster.fun = "louvain" uses Louvain community detection on a nearest-neighbor graph.
# - colLabels(sce) stores cluster ids for each cell in the SingleCellExperiment object in colData(sce)$label
# - use table() to show cluster sizes

# calculate dimensionality reduction (UMAP), based on PCA
# - runUMAP embeds cells into 2D for visualization using the PCA coordinates by default
# - this is deterministic if you set seeds / use deterministic UMAP settings (not done here)
# and plot the clusters on the UMAP, colored by label (cluster)


# Challenge 1 -----
# Our clusters look semi-reasonable, but what if we wanted to make them less granular?
# Look at the help documentation for ?clusterCells and ?NNGraphParam to find out what we’d need to change to get fewer, larger clusters.


# Marker Gene Detection ------

# change rownames to SYMBOL so downstream functions print gene symbols instead of Ensembl IDs
# - rowData(sce)$SYMBOL is supplied by the example dataset; real data may require mapping
# use scoreMarkers() to identify marker genes
# - scoreMarkers computes differential expression statistics for each cluster vs others.
# - It returns a list where each element corresponds to a cluster and contains a table
#   of marker stats such as log-fold change, AUC, p-values, etc.

# explore the output of scoreMarkers() 


# plot few marker genes for cluster 1

# order the marker genes based on mean.AUC, largest on top
# - mean.AUC is a measure of how well a gene distinguishes the cluster from others
# select top 6 marker genes
# plot their expression across clusters
# - plotExpression will show per-cell expression grouped by the 'x' column (here 'label')

# Challenge 2 -----
# Looking at the last plot, what clusters are most difficult to distinguish from cluster 1? 
# Now re-run the UMAP plot from the previous section. Do the difficult-to-distinguish clusters make sense?


# Cell type annotation -----
# Assigning cell labels from reference data -----
# based on the highest Spearman rank correlations between reference and sample cells
# - SingleR compares each query cell to reference cell profiles and assigns the label
#   with the highest similarity. It produces raw scores and a pruned label set.

# get reference atlas, sample 29 from EmbryoAtlasData()

# subset to 1000 cells

# check which cell types are present in reference

# pre-process reference (same normalization strategy as query)

# some cleaning, remove lowly abundant cell types
# - keep only reference cell types that have at least 10 cells (arbitrary threshold)

# change ref rownames to SYMBOL to match query

# restrict to genes present in both datasets (sce and ref)
# - SingleR requires the same features in test and ref; intersect to be safe

# convert sparse arrays to dense (for SingleR)
# - SingleR expects matrices; if your assays are sparse (e.g. dgCMatrix) convert first.
# - For very large datasets converting to dense may blow memory; use care in practice.

# run SingleR 
# - test: query matrix (cells in columns), ref: reference matrix
# - labels: vector of labels for the reference columns (ref$celltype)
# - SingleR returns per-cell scores and final/pruned labels (res$labels, res$pruned.labels)


# visualize SingleR scores
# one can see similarity to other cell types via heatmap of scores

# check cell annotation distribution in the clusters
# - produce a contingency table of SingleR annotations vs clustering labels
# - log2 scaling + pseudocount for visualization stability

# check also the present labels in our dataset (in sce: celltype.mapped)
# - this cross-checks SingleR results against any pre-existing mapped labels


# Challenge 3 -----
# Assign the SingleR annotations as a column in the colData for the query object sce.


# Assigning cell labels from gene sets -----

# Let's create a set of markers for our reference dataset
# based on pairwise Wilcoxon tests to find genes upregulated in each cell type
# getTopMarkers aggregates pairwise results and returns top genes per group

# we will convert the list of marker genes to a GeneSetCollection
# - GeneSetCollection is the container expected by many gene-set scoring functions

# we will use AUCell package for per-cell gene set scoring -> compute rankings
# - AUCell_buildRankings ranks genes within each cell by expression (required input for AUCell_calcAUC)
# - as.matrix(counts(sce)) uses the raw counts; AUCell recommends using counts for ranking

# calculate AUC scores for each gene set in each cell
# - cell.aucs is an AUCell-class object containing AUCs in an assay-like structure

# transpose to get a regular matrix: cells x gene-sets



# assign the max scoring identity to each cell 
# - max.col finds the index of the maximum score in each row (cell)
# compare to previously mapped celltype (optional)

# QC of the AUCs 
# - explore thresholds and diagnostic histograms for the first 9 gene sets
# - AUCell_exploreThresholds can optionally assign thresholds (assign = TRUE)

# Challenge 4 -----
# Inspect the diagnostics for the next nine cell types. Do they look okay?


# Exercises -----

## Exercise 1: Clustering -----
# try Leiden clustering with a specified resolution and objective


## Exercise 2: Reference marker genes -----
# score markers on the reference dataset (grouped by celltype)

# plot labels on pre-computed UMAP 


# Repetitive work -> write a helper function to extract top marker symbols
# using marker data.frame and number of top markers as two input parameters
# output should be a vector with rownames of the input data.frame with length 
# matching the number of top markers asked in the input

# proportion of overlap between two top-100 marker lists

## Extension Challenge 1: Group pair comparisons -----
# Why do you think marker genes are found by aggregating pairwise comparisons rather than iteratively comparing each cluster to all other clusters?

## Extension Challenge 2: Parallelizing SingleR -----
# SingleR can be computationally expensive. How do you set it to run in parallel?


## Extension Challenge 3: Critical inspection of diagnostics -----
# The first set of AUCell diagnostics don’t look so good for some of the examples here. Which ones? Why?


# Session Info ----
# sessionInfo() or modern version: sessioninfo::session_info()
