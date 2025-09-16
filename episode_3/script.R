# 3 Cell type annotation -----

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

library(AUCell)
library(MouseGastrulationData)
library(SingleR)
library(bluster)
library(scater)
library(scran)
library(pheatmap)
library(GSEABase)


# Data retrieval -----

# load WTChimeraData processed dataset (sample 5) and store as a SingleCellExperiment (sce)
# - WTChimeraData(): convenience wrapper from MouseGastrulationData that returns
#   an experiment-like object. Here we request 'processed' counts (already
#   filtered/normalized-ish depending on the package maintainers' preprocessing).
# - We print the object to inspect dimensions and metadata.
sce <- WTChimeraData(samples = 5, type = "processed")
sce

# select 1000 cells for a quicker interactive run (set.seed for reproducibility)
# - sampling reduces compute time for this tutorial while preserving structure.
set.seed(123)
ind <- sample(ncol(sce), 1000)
sce <- sce[,ind]
sce


# Pre-processing ----- 
# normalize counts
# - logNormCounts from scater computes size-factor based normalization followed by log-transformation
# - it creates an assay named "logcounts" in the SingleCellExperiment object
sce <- logNormCounts(sce)
# and run PCA
# - runPCA uses the 'logcounts' by default; here we do not pre-select HVGs so it will
#   by default use a set of most variable genes (often implemented inside runPCA)
# - PCA reduces dimensionality for visualization and neighborhood-based operations
sce <- runPCA(sce)
sce


# Clustering -----
# simple clustering based on PCA, using clusterCells function
# - BLUSPARAM defines graph construction + clustering algorithm settings.
# - cluster.fun = "louvain" uses Louvain community detection on a nearest-neighbor graph.
# - colLabels(sce) stores cluster ids for each cell in the SingleCellExperiment object in colData(sce)$label
# - use table() to show cluster sizes
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))
table(colLabels(sce))
table(colData(sce)$label)

# calculate dimensionality reduction (UMAP), based on PCA
# - runUMAP embeds cells into 2D for visualization using the PCA coordinates by default
# - this is deterministic if you set seeds / use deterministic UMAP settings (not done here)
sce <- runUMAP(sce, dimred = "PCA")
# and plot the clusters on the UMAP, colored by label (cluster)
plotReducedDim(sce, "UMAP", color_by = "label")


# Challenge 1 -----
# Our clusters look semi-reasonable, but what if we wanted to make them less granular?
# Look at the help documentation for ?clusterCells and ?NNGraphParam to find out what we’d need to change to get fewer, larger clusters.
# - Typical parameters that change granularity: 'k' for nearest neighbors, resolution-like parameters
# - Increasing k or using a lower resolution often produces fewer, larger clusters.

sce$clust2 <- clusterCells(sce, use.dimred = "PCA",
                           BLUSPARAM = NNGraphParam(cluster.fun = "louvain",
                                                    k = 30))
# - Here we increased k from the default (often 10-20) to 30; this enlarges the "neighborhood"
#   and tends to merge small clusters into larger ones.
plotReducedDim(sce, "UMAP", color_by = "clust2")


# Marker Gene Detection ------

# change rownames to SYMBOL so downstream functions print gene symbols instead of Ensembl IDs
# - rowData(sce)$SYMBOL is supplied by the example dataset; real data may require mapping
sce
rownames(sce) <- rowData(sce)$SYMBOL
sce
# use scoreMarkers() to identify marker genes
# - scoreMarkers computes differential expression statistics for each cluster vs others.
# - It returns a list where each element corresponds to a cluster and contains a table
#   of marker stats such as log-fold change, AUC, p-values, etc.
markers <- scoreMarkers(sce)

# explore the output of scoreMarkers() 
markers
markers[[1]]


# plot few marker genes for cluster 1
c1_markers <- markers[[1]]

# order the marker genes based on mean.AUC, largest on top
# - mean.AUC is a measure of how well a gene distinguishes the cluster from others
ord <- order(c1_markers$mean.AUC, decreasing = TRUE)
# select top 6 marker genes
top.markers <- head(rownames(c1_markers[ord,]))
# plot their expression across clusters
# - plotExpression will show per-cell expression grouped by the 'x' column (here 'label')
plotExpression(sce, 
               features = top.markers, 
               x        = "label", 
               color_by = "label")

# Challenge 2 -----
# Looking at the last plot, what clusters are most difficult to distinguish from cluster 1? 
# Now re-run the UMAP plot from the previous section. Do the difficult-to-distinguish clusters make sense?
plotReducedDim(sce, "UMAP", color_by = "label")


# Cell type annotation -----
# Assigning cell labels from reference data -----
# based on the highest Spearman rank correlations between reference and sample cells
# - SingleR compares each query cell to reference cell profiles and assigns the label
#   with the highest similarity. It produces raw scores and a pruned label set.

# get reference atlas, sample 29 from EmbryoAtlasData()
ref <- EmbryoAtlasData(samples = 29)
ref

# subset to 1000 cells
set.seed(123)
ind <- sample(ncol(ref), 1000)
ref <- ref[,ind]
ref

# check which cell types are present in reference
tab <- sort(table(ref$celltype), decreasing = TRUE)
tab

# pre-process reference (same normalization strategy as query)
ref <- logNormCounts(ref)

# some cleaning, remove lowly abundant cell types
# - keep only reference cell types that have at least 10 cells (arbitrary threshold)
abu.ct <- names(tab)[tab >= 10]
ind <- ref$celltype %in% abu.ct
ref <- ref[,ind] 
ref

# change ref rownames to SYMBOL to match query
rownames(ref) <- rowData(ref)$SYMBOL
ref

# restrict to genes present in both datasets (sce and ref)
# - SingleR requires the same features in test and ref; intersect to be safe
shared_genes <- intersect(rownames(sce), rownames(ref))
sce <- sce[shared_genes,]
ref <- ref[shared_genes,]

# convert sparse arrays to dense (for SingleR)
# - SingleR expects matrices; if your assays are sparse (e.g. dgCMatrix) convert first.
# - For very large datasets converting to dense may blow memory; use care in practice.
sce.mat <- as.matrix(assay(sce, "logcounts"))
ref.mat <- as.matrix(assay(ref, "logcounts"))

# run SingleR 
# - test: query matrix (cells in columns), ref: reference matrix
# - labels: vector of labels for the reference columns (ref$celltype)
# - SingleR returns per-cell scores and final/pruned labels (res$labels, res$pruned.labels)
res <- SingleR(test = sce.mat, 
               ref = ref.mat,
               labels = ref$celltype)
res


# visualize SingleR scores
# one can see similarity to other cell types via heatmap of scores
plotScoreHeatmap(res)

# check cell annotation distribution in the clusters
# - produce a contingency table of SingleR annotations vs clustering labels
tab <- table(anno = res$pruned.labels, cluster = colLabels(sce))
# - log2 scaling + pseudocount for visualization stability
pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))

# check also the present labels in our dataset (in sce: celltype.mapped)
# - this cross-checks SingleR results against any pre-existing mapped labels
tab <- table(res$pruned.labels, sce$celltype.mapped)
pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))


# Challenge 3 -----
# Assign the SingleR annotations as a column in the colData for the query object sce.
# - This stores the assigned labels directly alongside other per-cell metadata
sce$SingleR_label <- res$pruned.labels


# Assigning cell labels from gene sets -----

# Let's create a set of markers for our reference dataset
# based on pairwise Wilcoxon tests to find genes upregulated in each cell type
wilcox.z <- pairwiseWilcox(ref, ref$celltype, lfc = 1, direction = "up")
# getTopMarkers aggregates pairwise results and returns top genes per group
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, 
                           pairwise = FALSE, n = 50)
lengths(markers.z) # show how many markers for each group

# we will convert the list of marker genes to a GeneSetCollection
# - GeneSetCollection is the container expected by many gene-set scoring functions
all.sets <- lapply(names(markers.z), 
                   function(x) GeneSet(markers.z[[x]], setName = x))
all.sets <- GeneSetCollection(all.sets)
all.sets

# we will use AUCell package for per-cell gene set scoring -> compute rankings
# - AUCell_buildRankings ranks genes within each cell by expression (required input for AUCell_calcAUC)
# - as.matrix(counts(sce)) uses the raw counts; AUCell recommends using counts for ranking
rankings <- AUCell_buildRankings(as.matrix(counts(sce)),
                                 plotStats = FALSE, verbose = FALSE)

# calculate AUC scores for each gene set in each cell
# - cell.aucs is an AUCell-class object containing AUCs in an assay-like structure
cell.aucs <- AUCell_calcAUC(all.sets, rankings)

# transpose to get a regular matrix: cells x gene-sets
results <- t(assay(cell.aucs))

head(results)


# assign the max scoring identity to each cell 
# - max.col finds the index of the maximum score in each row (cell)
new.labels <- colnames(results)[max.col(results)]
# compare to previously mapped celltype (optional)
tab <- table(new.labels, sce$celltype.mapped)
tab[1:4,1:4]

# QC of the AUCs 
# - explore thresholds and diagnostic histograms for the first 9 gene sets
# - AUCell_exploreThresholds can optionally assign thresholds (assign = TRUE)
par(mfrow = c(3,3))
AUCell_exploreThresholds(cell.aucs[1:9], plotHist = TRUE, assign = TRUE) 

# Challenge 4 -----
# Inspect the diagnostics for the next nine cell types. Do they look okay?
par(mfrow = c(3,3))
AUCell_exploreThresholds(cell.aucs[10:18], plotHist = TRUE, assign = TRUE) 


# Exercises -----

## Exercise 1: Clustering -----
# The Leiden algorithm is similar to the Louvain algorithm, 
# but it is faster and has been shown to result in better connected 
# communities. Modify the above call to clusterCells to carry out 
# the community detection with the Leiden algorithm instead. 
# Visualize the results in a UMAP plot.
# Hint: The NNGraphParam constructor has an argument cluster.args.
#       This allows to specify arguments passed on to the cluster_leiden 
#       function from the igraph package. Use the cluster.args argument to 
#       parameterize the clustering to use modularity as the objective function 
#       and a resolution parameter of 0.5.
arg_list <- list(objective_function = "modularity",
                 resolution_parameter = .5)
sce$leiden_clust <- clusterCells(sce, use.dimred = "PCA",
                                 BLUSPARAM = NNGraphParam(cluster.fun = "leiden", 
                                                          cluster.args = arg_list))

plotReducedDim(sce, "UMAP", color_by = "leiden_clust")

## Exercise 2: Reference marker genes -----
# Identify the marker genes in the reference single cell 
# experiment, using the celltype labels that come with the 
# dataset as the groups. Compare the top 100 marker genes of 
# two cell types that are close in UMAP space. 
# Do they share similar marker sets?
markers <- scoreMarkers(ref, groups = ref$celltype)
markers

# plot labels on pre-computed UMAP (ref object must have UMAP computed beforehand)
# - this helps to visually inspect where reference cell types lie in reduced space
plotReducedDim(ref, dimred = "umap", color_by = "celltype") 

# Repetitive work -> write a helper function to extract top marker symbols
# using marker data.frame and number of top markers as two input parameters
# output should be a vector with rownames of the input data.frame with length 
# matching the number of top markers asked in the input
order_marker_df <- function(m_df, n = 100) {
    ord <- order(m_df$mean.AUC, decreasing = TRUE)
    rownames(m_df[ord,][1:n,])
}

x <- order_marker_df(markers[["Erythroid2"]])

y <- order_marker_df(markers[["Erythroid3"]])

# proportion of overlap between two top-100 marker lists
length(intersect(x,y)) / 100

## Extension Challenge 1: Group pair comparisons -----
# Why do you think marker genes are found by aggregating pairwise 
# comparisons rather than iteratively comparing each cluster to 
# all other clusters?
# - Short answer: pairwise comparisons are better at finding genes that distinguish a target group from each individual alternative group.
#   Aggregation (e.g. taking genes consistently up in all pairwise comparisons) reduces false positives driven by a single contrasting group.

## Extension Challenge 2: Parallelizing SingleR -----
# SingleR can be computationally expensive. How do you set it to run in parallel?
# - Use BiocParallel::MulticoreParam (or other BPPARAM backends) and pass via BPPARAM argument
library(BiocParallel)
my_bpparam <- MulticoreParam(workers = 4)
res2 <- SingleR(test = sce.mat, 
                ref = ref.mat,
                labels = ref$celltype,
                BPPARAM = my_bpparam)


## Extension Challenge 3: Critical inspection of diagnostics -----
# The first set of AUCell diagnostics don’t look so good for some of the examples here. Which ones? Why?
# - Look for multimodal or very flat histograms: these indicate poor enrichment or lack of a clear "on" population.
# - Also check whether the top-scoring gene sets are biologically plausible for the cell types in your sample.


# Session Info ----
# sessionInfo() or modern version: sessioninfo::session_info()
sessionInfo()
 
