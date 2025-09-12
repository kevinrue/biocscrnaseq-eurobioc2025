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

# remove problematic cells (doublets, stripped nuclei)
# - flagged in celltype.mapped metadata
# - exclude these from downstream analysis

# randomly select 50% of the cells from each sample (downsampling for speed)
# - use tapply to split by sample, then randomly sample 50% per sample
# - set.seed ensures reproducibility


# normalize the data
# - logNormCounts applies size-factor normalization followed by log-transform

# select HVGs
# - modelGeneVar decomposes variance into technical vs biological components
# - block = sce$sample ensures variance modeling accounts for sample-specific differences/batch

# calculate dimensionality reductions: PCA and t-SNE
# - PCA on top HVGs (or bio>0 genes), reduce dimensionality for denoising
# - runTSNE using PCA coordinates as input

# convert sample id to factor for easier plotting

# plot t-SNE coloured by sample

# define a custom palette with distinct colors for many cell types
# - adapted from pals::polychrome palette
# - ensures sufficient visual distinction when plotting many categories
color_vec <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
               "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", 
               "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D", 
               "#B10DA1", "#3B00FB", "#1CBE4F", "#FA0087", "#333333", "#F7E1A0", 
               "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", 
               "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#FBE426")

# Challenge 1 -----
# It seems like samples 5 and 6 are separated from the others in expression space.
# Given the groups of cells, why might this be? Likely due to tomato+ vs tomato- condition
# rather than just biological replicate. The factor is the tdTomato injection.


# Correcting Batch Effects -----
# correct sample-specific effects with correctExperiments using FastMNN
# - batch = sample column

# plot corrected t-SNE coloured by batch (sample)

# plot corrected t-SNE coloured by cell type

# Challenge 2 -----
# True or False: after batch correction, no batch-level information is present.
# In practice, some residual batch effect may remain; correction reduces but doesn’t erase.


# Differential Expression -----

# Pseudo-bulk Samples -----
# Aggregate counts across cell types and biological samples to form pseudo-bulk replicates
# - aggregateAcrossCells groups by celltype.mapped and sample metadata
# - each column corresponds to a unique combination


# Differential Expression Analysis -----
# Use edgeR to model counts with negative binomial distribution
# Focus on Mesenchyme cells

# filter out pseudo-bulk samples with <10 contributing cells

# filter out low-expression genes (noise)

# normalization factors (TMM)

# quality control: MA plots per sample

# densities of logCPM distributions

# MDS plot to check grouping by tomato condition


# design matrix
# - factors: pool (batch), tomato injection (condition)

# estimate dispersion parameters (common, trended, tagwise)

# fit quasi-likelihood GLM

# check prior distribution parameters

# plot genewise dispersion estimates

# test for differential expression due to tomato injection

# summarize DE results

# top 10 DE genes

# automate DE across all cell types



# Challenge 3 -----
# Low p-values suggest significance, but effect size (logFC) matters.
# logFC = log2 fold-change in expression between tomato+ vs tomato-.


# Differential Abundance -----

# test differences in cell abundances (counts of cells per sample)

# match sample metadata

# build DGEList with cell counts

# design matrix with pool and tomato

# fit GLM for DA analysis


# Background On Compositional Effect -----
# If total counts differ between conditions, relative abundances may shift spuriously.

# Normalize counts with TMM for compositional bias



# Testing Against a log-fold Change Threshold -----
# Apply glmTreat to require abs(logFC) > 1 for significance


# Exercise 1 -----
# Heatmap of abundances table with pheatmap.
# Can overlay sample metadata (tomato, pool) as column annotations.





# Exercise 2 -----
# Re-run pseudoBulkDGE without pool factor; compare to original.
# Merge results for Erythroid3 and compare logFC estimates and p-values.








# Extension Challenge 1: Group effects -----
# Multiple independent replicates per group are crucial for reliable correction.
# Without replicates, batch effect and condition are confounded.


# Session Info ----
# sessionInfo() or modern version: sessioninfo::session_info()
