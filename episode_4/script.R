# 4 Multi-sample analyses -----

# Setup and data exploration -----

# load libraries
# - MouseGastrulationData: example single-cell datasets from mouse gastrulation
# - batchelor: batch correction methods for single-cell data (e.g. fastMNN)
# - edgeR: differential expression and abundance analysis (negative binomial)
# - scater: utilities for quality control, visualization, preprocessing
# - ggplot2: general-purpose plotting library for visualizations
# - scran: preprocessing, clustering, HVG detection, pseudo-bulk utilities
# - pheatmap: simple heatmap visualization for counts or summary tables
# - scuttle: basic single-cell utilities (normalization, QC, etc.)
library(MouseGastrulationData)
library(batchelor)
library(edgeR)
library(scater)
library(ggplot2)
library(scran)
library(pheatmap)
library(scuttle)


# Data retrieval and pre-processing -----

# load WTChimeraData processed dataset (samples 5:10) and store as a SingleCellExperiment (sce)
# - Sample 5: E8.5 injected cells (tdTomato+, pool 3)
# - Sample 6: E8.5 host cells (tdTomato-, pool 3)
# - Sample 7: E8.5 injected cells (tdTomato+, pool 4)
# - Sample 8: E8.5 host cells (tdTomato-, pool 4)
# - Sample 9: E8.5 injected cells (tdTomato+, pool 5)
# - Sample 10: E8.5 host cells (tdTomato-, pool 5)
#
# - WTChimeraData() returns processed SingleCellExperiment objects (already filtered)
# - Inspect dimensions and metadata; colData() contains per-cell metadata including sample, pool, celltype
sce <- WTChimeraData(samples = 5:10, type = "processed")
sce
colData(sce)

# remove problematic cells (doublets, stripped nuclei)
# - flagged in celltype.mapped metadata
# - exclude these from downstream analysis
drop <- sce$celltype.mapped %in% c("stripped", "Doublet")
sce <- sce[,!drop]
sce

# randomly select 50% of the cells from each sample (downsampling for speed)
# - use tapply to split by sample, then randomly sample 50% per sample
# - set.seed ensures reproducibility
set.seed(29482)
idx <- unlist(tapply(colnames(sce), sce$sample, function(x) {
    perc <- round(0.50 * length(x))
    sample(x, perc)
}))
sce <- sce[,idx]
sce


# normalize the data
# - logNormCounts applies size-factor normalization followed by log-transform
sce <- logNormCounts(sce)

# select HVGs
# - modelGeneVar decomposes variance into technical vs biological components
# - block = sce$sample ensures variance modeling accounts for sample-specific differences
dec <- modelGeneVar(sce, block = sce$sample)
chosen.hvgs <- dec$bio > 0

# calculate dimensionality reductions: PCA and t-SNE
# - PCA on top HVGs (or bio>0 genes), reduce dimensionality for denoising
# - runTSNE using PCA coordinates as input
sce <- runPCA(sce, subset_row = chosen.hvgs, ntop = 1000)
sce <- runTSNE(sce, dimred = "PCA")

# convert sample id to factor for easier plotting
sce$sample <- as.factor(sce$sample)

# plot t-SNE coloured by sample
plotTSNE(sce, colour_by = "sample")

# define a custom palette with distinct colors for many cell types
# - adapted from pals::polychrome palette
# - ensures sufficient visual distinction when plotting many categories
color_vec <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
               "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", 
               "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D", 
               "#B10DA1", "#3B00FB", "#1CBE4F", "#FA0087", "#333333", "#F7E1A0", 
               "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", 
               "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#FBE426")

plotTSNE(sce, colour_by = "celltype.mapped") +
    scale_color_manual(values = color_vec) +
    theme(legend.position = "bottom")

# Challenge 1 -----
# It seems like samples 5 and 6 are separated off from the others in gene 
# expression space. Given the group of cells in each sample, why might this 
# make sense versus some other pair of samples? What is the factor presumably 
# leading to this difference?


# Correcting Batch Effects -----
# correct sample-specific effects with correctExperiments using FastMNN
# - batch = sample column
set.seed(10102)
merged <- correctExperiments(
    sce, 
    batch = sce$sample, 
    subset.row = chosen.hvgs,
    PARAM = FastMnnParam(
        merge.order = list(
            list(1,3,5), # WT replicates
            list(2,4,6)  # tdTomato replicates
        )
    )
)
merged <- runTSNE(merged, dimred = "corrected")

# plot corrected t-SNE coloured by batch (sample)
plotTSNE(merged, colour_by = "batch")

# plot corrected t-SNE coloured by cell type
plotTSNE(merged, colour_by = "celltype.mapped") +
    scale_color_manual(values = color_vec) +
    theme(legend.position = "bottom")

# Challenge 2 -----
# True or False: after batch correction, no batch-level information is present. 

# Differential Expression -----

# Pseudo-bulk Samples -----
# Aggregate counts across cell types and biological samples to form pseudo-bulk replicates
# - aggregateAcrossCells groups by celltype.mapped and sample metadata
# - each column corresponds to a unique combination
summed <- aggregateAcrossCells(
    merged, 
    id = colData(merged)[,c("celltype.mapped", "sample")]
)
summed


# Differential Expression Analysis -----
# Use edgeR to model counts with negative binomial distribution
# Focus on Mesenchyme cells
current <- summed[, summed$celltype.mapped == "Mesenchyme"]
y <- DGEList(counts(current), samples = colData(current))
y

# filter out pseudo-bulk samples with <10 contributing cells
discarded <- current$ncells < 10
y <- y[,!discarded]
summary(discarded)

# filter out low-expression genes (noise)
keep <- filterByExpr(y, group = current$tomato)
y <- y[keep,]
summary(keep)

# normalization factors (TMM)
y <- calcNormFactors(y)
y$samples

# quality control: MA plots per sample
par(mfrow = c(2,3))
for (i in seq_len(ncol(y))) {
    plotMD(y, column = i)
}

# densities of logCPM distributions
par(mfrow = c(1,1))
plotDensities(cpm(y, log=TRUE))

# MDS plot to check grouping by tomato condition
par(mfrow = c(1,1))
limma::plotMDS(cpm(y, log = TRUE), 
               col = ifelse(y$samples$tomato, "red", "blue"))


# design matrix
# - factors: pool (batch), tomato injection (condition)
design <- model.matrix(~factor(pool) + factor(tomato),
                       data = y$samples)
design

# estimate dispersion parameters (common, trended, tagwise)
y2 <- estimateDisp(y, design)
summary(y2$trended.dispersion)
plotBCV(y2)

# fit quasi-likelihood GLM
fit <- glmQLFit(y, design, robust = TRUE, abundance.trend = TRUE)

# check prior distribution parameters
summary(fit$s2.prior)
summary(fit$df.prior)

# plot genewise dispersion estimates
plotQLDisp(fit)

# test for differential expression due to tomato injection
res <- glmQLFTest(fit, coef = ncol(design))

# summarize DE results
summary(decideTests(res))

# top 10 DE genes
topTags(res)

# automate DE across all cell types
summed.filt <- summed[,summed$ncells >= 10]

de.results <- pseudoBulkDGE(
    summed.filt, 
    label = summed.filt$celltype.mapped,
    design = ~factor(pool) + tomato,
    coef = "tomatoTRUE",
    condition = summed.filt$tomato 
)

cur.results <- de.results[["Allantois"]]
cur.results[order(cur.results$PValue),]

# Challenge 3 -----
# Clearly some of the results have low p-values. 
# What about the effect sizes? What does logFC stand for?


# Differential Abundance -----

# test differences in cell abundances (counts of cells per sample)
abundances <- table(merged$celltype.mapped, merged$sample) 
abundances <- unclass(abundances) 
abundances

# match sample metadata
extra.info <- colData(merged)[match(colnames(abundances), merged$sample),]

# build DGEList with cell counts
y.ab <- DGEList(abundances, samples = extra.info)

# design matrix with pool and tomato
design <- model.matrix(~factor(pool) + factor(tomato), y.ab$samples)

# estimate dispersion (not needed since edgeR 4.0)
# y.ab <- estimateDisp(y.ab, design, trend = "none")

# fit GLM for DA analysis
fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
res.ab <- glmQLFTest(fit.ab, coef = ncol(design))

summary(decideTests(res.ab))
topTags(res.ab, n=Inf)

# Background On Compositional Effect -----
# If total counts differ between conditions, relative abundances may shift spuriously.

# Normalize counts with TMM for compositional bias
y.ab2 <- calcNormFactors(y.ab)
y.ab2$samples$norm.factors

# estimate dispersion (not needed since edgeR 4.0)
# y.ab2 <- estimateDisp(y.ab2, design, trend = "none")

fit.ab2 <- glmQLFit(y.ab2, design, robust = TRUE, abundance.trend = FALSE)
res.ab2 <- glmQLFTest(fit.ab2, coef = ncol(design))

summary(decideTests(res.ab2))
topTags(res.ab2, n = Inf)

# Testing Against a log-fold Change Threshold -----
# Apply glmTreat to require abs(logFC) > 1 for significance
res.ab.lfc <- glmTreat(fit.ab, coef = ncol(design), lfc = 1)
summary(decideTests(res.ab.lfc))
topTags(res.ab.lfc, n = Inf)


# Exercise 1 -----
# Use the pheatmap package to create a heatmap of the abundances table. 
# Does it comport with the model results?
# 
# Hint: You can simply hand pheatmap() a matrix as its only argument. 
#       pheatmap() has a million options you can adjust, but the defaults are 
#       usually pretty good. Try to overlay sample-level information with the 
#       annotation_col argument for an extra challenge.

pheatmap(y.ab$counts)

anno_df <- y.ab$samples[,c("tomato", "pool")]
anno_df$pool = as.character(anno_df$pool)
anno_df$tomato <- ifelse(anno_df$tomato,
                         "tomato+",
                         "tomato-")
pheatmap(y.ab$counts,
         annotation_col = anno_df)

pheatmap(log1p(y.ab$counts),
         annotation_col = anno_df)

pheatmap(log1p(y.ab$counts), scale="row",
         annotation_col = anno_df)

# Exercise 2 -----
# Try re-running the pseudobulk DGE without the pool factor in the design specification. 
# Compare the logFC estimates and the distribution of p-values for the Erythroid3 cell type.
#
# Hint: After running the second pseudobulk DGE, you can join the two DataFrames 
#       of Erythroid3 statistics using the merge() function. You will need to create 
#       a common key column from the gene IDs.

de.results2 <- pseudoBulkDGE(
    summed.filt, 
    label = summed.filt$celltype.mapped,
    design = ~tomato,
    coef = "tomatoTRUE",
    condition = summed.filt$tomato 
)

eryth1 <- de.results$Erythroid3
eryth2 <- de.results2$Erythroid3

eryth1$gene <- rownames(eryth1)
eryth2$gene <- rownames(eryth2)

comp_df <- merge(eryth1, eryth2, by = 'gene')
comp_df <- comp_df[!is.na(comp_df$logFC.x),]

ggplot(comp_df, aes(logFC.x, logFC.y)) + 
    geom_abline(lty = 2, color = "grey") +
    geom_point() 

pval_df <- reshape(comp_df[,c("gene", "PValue.x", "PValue.y")],
                   direction = "long", 
                   v.names = "Pvalue",
                   timevar = "pool_factor",
                   times = c("with pool factor", "no pool factor"),
                   varying = c("PValue.x", "PValue.y"))

ggplot(pval_df, aes(Pvalue)) + 
    geom_histogram(boundary = 0,
                   bins = 30) + 
    facet_wrap("pool_factor")

# Extension Challenge 1: Group effects -----
# Having multiple independent samples in each experimental group is always helpful, 
# but it is particularly important when it comes to batch effect correction. Why?


# Session Info ----
# sessionInfo() or modern version: sessioninfo::session_info()
sessionInfo()
