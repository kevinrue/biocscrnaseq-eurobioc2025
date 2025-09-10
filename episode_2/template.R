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

