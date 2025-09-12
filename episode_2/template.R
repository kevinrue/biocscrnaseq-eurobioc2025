#################################
# Setup and experimental design #
#################################

## Load the following libraries
## - MouseGastrulationData
## - DropletUtils
## - ggplot2
## - EnsDb.Mmusculus.v79
## - scuttle
## - scater
## - scran
## - scDblFinder



## Use WTChimeraData() to load the 'raw' data for the fifth sample.
## Assign the output to a new object called 'sce'.



## The raw data above came back as a list of one element.
## Extract that first element and reassign it to 'sce'.



## Display the 'sce' object.



######################
# Droplet processing #
######################

## Use barcodeRanks() to compute the barcode rank statistics for the raw counts in the 'sce' object.



## Remove rows with duplicated values in the 'rank' column (for plotting speed later)



## Make a ggplot figure with:
## - the rank of each barcode along the x-axis
## - the total UMI for each barcode along the y-axis



## Update the plot to add horizontal lines indicating the knee and inflection points computed by barcodeRanks().
## That means:
## Create a data.frame (called 'line_df') with the following columns:
## - 'cutoff' which contains the names of the two lines (i.e., 'inflection' and 'knee')
## - 'value' which contains the position where each line will be positioned
## Copy the ggplot code above, using 'geom_hline()' to add the lines.



## CHALLENGE ##
## What is the median number of total counts in the raw data?


##############################
# Testing for empty droplets #
##############################

## use emptyDrops() to test for empty droplets using the raw counts of the 'sce' object.
## assign the result to a new object named 'e.out'
## Tip: Set the random seed to ensure reproducible results! (e.g., 100)



## Use summary() on the 'FDR' column of 'e.out' to see how many barcodes have an FDR equal or lower than 0.001.



## Subset 'sce' to the barcodes with FDR equal or lower than 0.001 and reassign it to 'sce'.



###################
# Quality control #
###################

## Use mapIds() to map get the chromosome name for each feature in 'sce'.
## Assign the result to a new object named 'chr.loc'.



## Use 'chr.loc' to identify which features are located on the mitochondrial chromosome.
## Assign the result to a new object called 'is.mito'.
## Hint: the mitochondrial chromosome is named 'MT'.



## Use perCellQCMetrics() to compute cell-wise quality control metrics.
## Use 'is.mito' to compute an additional set of quality control metrics on those features.
## Assign the result to a new object named 'df'.



## Add 'df' to the column-wise metadata of 'sce'.



## Use table() to see how many barcodes have a UMI sum strictly below 10,000.



## Use table() to see how many barcodes have a mitochondrial fraction strictly above 10%.
## Hint: the mitochondrial fraction is reported as numbers between 0 and 100.



## Use summary() to see the distribution of the number of detected features per barcode.



## Use summary() to see the distribution of mitochondrial fraction per barcode.



## Use perCellQCFilters() to generate a DataFrame of logical values
## indicating which cells should be discarded based on the metrics in 'df'.
## Use sub.fields = "subsets_Mito_percent" to include the mitochondrial fraction.
## Assign the result to a new object named 'reasons'.



## Add a new column named 'discard' in the column-wise metadata
## The values should be a copy of the 'discard' column in 'reasons'.



## CHALLENGE ##
## Maybe our sample preparation was poor and we want the QC to be more strict.
## How could we change the set the QC filtering to use 2.5 MADs as the threshold for outlier calling?


####################
# Diagnostic plots #
####################

## Use plotColData() to visualise
## the distribution of UMI sum along the Y axis
## coloured by the 'discard' column of cell-wise metadata
## Optionally, set the title of the plot to 'Total count'



## Same as above, except
## visualise the number of detected features along the Y axis
## set the title to 'Detected features'



## Same as above, except
## visualise the mitochondrial fraction along the Y axis
## set the title to 'Mito percent'



## Subset 'sce' to the barcodes that are NOT marked for discarding
## and reassign it to 'sce'.



#################
# Normalization #
#################

## Use librarySizeFactors() to compute library size factors for the 'sce' object.
## Assign the result to a new object named 'lib.sf'.



## Use summary() to see the distribution of values per barcode.



## Similarly, use ggplot() and geom_histogram() to see the distribution as a plot.
## Hint: Store the size factors in a data.frame named 'sf_df', in a column named 'size_factor'.



##################################
# Normalization by deconvolution #
##################################

## Use quickCluster() to compute clusters of cells with similar expression profiles.
## Assign the result to a new object named 'clust'.
## Tip: Set the random seed to ensure reproducible results! (e.g., 100)



## Use table() to visualise the number of barcodes assigned to each cluster.



## Use pooledSizeFactors() to compute size factors by deconvolution.
## Use 'clust' to provide the precomputed cluster assignments.
## Assign the result to a new object named 'deconv.sf'.



## Use summary() to visualise the distribution of deconvoluted size factors per barcode.



## Use ggplot() to visualise the standard size factors and the deconvoluted size factors.
## Colour the points by the cluster assignments in 'clust'.
## Hints:
## - Add 'clust' and 'deconv.sf' to 'sf_df', in columns named 'deconv_sf' and 'clust', respectively.
## - Apply log10 transformation to both x- and y- axes.



## Use sizeFactors() to assign the deconvoluted size factors to the 'sce' object.



## Use logNormCounts() to normalise the counts in 'sce' using the deconvoluted size factors.



## CHALLENGE ##
## Fill in the blanks for normalization that uses simpler library size factors instead of deconvolution.


#####################
# Feature Selection #
#####################

## Use modelGeneVar() to model and decompose the gene-wise biological and technical variance.
## Assign the output to a new object called 'dec.sce'.



## Use metadata() to extract the fitted trend from the 'dec.sce' object.
## Assign the output to a new object called 'fit.sce'.



## Use ggplot() to visualise the variance and mean expression of each gene.
## Overlay the fitted trend of from 'fit.sce'.
## Hint: Create a data.frame named 'mean_var_df' with two columns named 'mean' and 'var'
## to store the 'mean' and 'var' elements of 'fit.sce'.



## Use getTopHVGs()to fetch the 1,000 most variable genes,
## based on the variance modelling statistics stored in 'dec.sce'.



## CHALLENGE ##
## Imagine you have data that were prepared by three people with varying level of experience,
## which leads to varying technical noise. How can you account for this blocking structure when selecting HVGs?



############################
# Dimensionality Reduction #
############################

## Use RunPCA() to compute a principal component analysis on 'sce',
## using only the set of variable genes selected above.
## Reassign the result to 'sce'.



## Use ggplot() to visualise the variance explained by each principal component.
## Hint: Create a new data.frame named 'pct_pca_df' with two columns named 'PC' and 'pct_var'
## to store the name of each principal component and the variance it explains, respectively.



## Use plotPCA() to visualise the cells along the first two principal components,
## coloured by UMI sum.



## Use plotReducedDim() to visualise each pair of principal components for the first three PCs.



## Use runTSNE() to compute a t-SNE layout based on the PCA results.
## Tip: Set the random seed to ensure reproducible results! (e.g., 100)



## Use plotTSNE() to visualise the t-SNE layout.



## Use runUMAP() to compute a UMAP layout based on the PCA results.
## Tip: Set the random seed to ensure reproducible results! (e.g., 111)



## Use plotUMAP() to visualise the UMAP layout.



## CHALLENGE ##
## Re-run the UMAP for the same sample starting from the pre-processed data (i.e. not type = "raw").
## What looks the same? What looks different?



##########################
# Doublet identification #
##########################

## Use computeDoubletDensity() to compute doublet scores for each cell in 'sce',
## using only the set of variable genes selected earlier.
## Tip: Set the argument 'dims' to 50, to match the number of principal components computed earlier.
## Assign the output to a new object called 'dbl.send'.


## Use summary() to see the distribution of doublet density values per barcode.



## Add the doublet density values to the cell-wise metadata of 'sce',
## in a new column named 'DoubletScore'.



## Use plotTSNE() to visualise the t-SNE layout coloured by doublet density.



## Use doubletThresholding() to classify cells into singlets and doublets,
## using the doublet density scores computed above.
## Assign the output to a new object called 'dbl.calls'.



## Use summary() to display the number of singlets and doublets called.



## Add the doublet calls to the cell-wise metadata of 'sce',
## in a new column named 'doublet'.



## Use plotColData() to visualise the distribution of doublet scores as a violin plot
## coloured by doublet status.



## Use plotTSNE() to visualise the t-SNE layout coloured by doublet call.



## Use plotColData() to visualise:
## - UMI sum on the x-axis
## - number of genes detected on the y-axis
## - coloured by doublet score



## Same as above, except
## - colour by doublet call



## EXERCISE 1 - NORMALISATION ##
## Here we used the deconvolution method implemented in scran based on a previous clustering step.
## Use the pooledSizeFactors to compute the size factors without considering a preliminary clustering.
## Compare the resulting size factors via a scatter plot.
## How do the results change?
## What are the risks of not including clustering information?



## EXERCISE 2 - PBMC DATA ##
## The DropletTestFiles package includes the raw output from Cell Ranger
## of the peripheral blood mononuclear cell (PBMC) dataset from 10X Genomics,
## publicly available from the 10X Genomics website.
## Repeat the analysis of this vignette using those data.

## The hint demonstrates how to identify, download, extract, and read
## the data starting from the help documentation of ?DropletTestFiles::listTestFiles,
## but try working through those steps on your own for extra challenge
## (they’re useful skills to develop in practice).



## EXTENSION CHALLENGE 1: SPIKE-INS ##
## Some sophisticated experiments perform additional steps so that they can estimate size factors from so-called “spike-ins”. Judging by the name, what do you think “spike-ins” are, and what additional steps are required to use them?



## EXTENSION CHALLENGE 2: BACKGROUND RESEARCH ##
## Run an internet search for some of the most highly variable genes we identified in the feature selection section.
## See if you can identify the type of protein they produce or what sort of process they’re involved in.
## Do they make biological sense to you?



## EXTENSION CHALLENGE 3: REDUCED DIMENSIONALITY REPRESENTATIONS ##
## Can dimensionality reduction techniques provide a perfectly accurate representation of the data?



###########
# THE END #
###########
