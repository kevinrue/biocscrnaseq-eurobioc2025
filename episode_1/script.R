################
# Bioconductor #
################

## Use install.packages() to install the package 'ggplot2' from CRAN.

install.packages("ggplot2")
# loading the package
library(ggplot2)
library("ggplot2")

## Use install.packages() to install the package 'BiocManager' from CRAN.

install.packages("BiocManager")

## Use BiocManager::install() to install the package 'SingleCellExperiment' from Bioconductor.

BiocManager::install("SingleCellExperiment")

## Navigate to <https://bioconductor.org/packages/release/bioc/html/scater.html>
## Then, copy paste instructions from web page to install the package 'scater'.

if (!require("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

BiocManager::install("scater")

## Use BiocManager::install() to update your packages to the latest version available.

BiocManager::install()

#################################
# The SingleCellExperiment class #
#################################

## Load the libraries 'SingleCellExperiment' and 'MouseGastrulationData'.

library(SingleCellExperiment)
library(MouseGastrulationData)

## Use WTChimeraData() to load the fifth sample.
## Assign the output to a new object called 'sce'.

?WTChimeraData
help(WTChimeraData)

sce <- WTChimeraData(samples = 5)

sce
print(sce)

## CHALLENGE ##
## Get the data for a different sample from WTChimeraData (other than the fifth one).

WTChimeraData(samples = 1)

## Display the names of the assays in the 'sce' object.

assayNames(sce)
names(assays(sce))

## Display the values in the first three rows and first three columns
## of the counts matrix of the 'sce' object.

counts(sce)[1:3, 1:3]

## Display the first three rows and first four columns of the column-wise metadata
## of the 'sce' object.

colData(sce)[1:3, 1:4]

## Display the first three rows and first two columns of the row-wise metadata
## of the 'sce' object.

rowData(sce)[1:3, 1:2]

## Use the '$' operator to add a new column named 'my_sum' to the column-wise metadata
## of the 'sce' object.
## The values in that column should be the sum for each column of the counts matrix.

sce$stage[1:3]

sce$my_sum <- colSums(counts(sce))

colData(sce)[1:3, ]

## CHALLENGE ##
## Add a column of gene-wise metadata to the rowData.

rowData(sce)$row_sum <- rowSums(counts(sce))

rowData(sce)[1:3,]

## Display the list of dimensionality reduction results present in the 'sce' object.

reducedDims(sce)

reducedDim(sce, "pca.corrected.E8.5")[1:3,]

## Load the 'scater' library.
## Then, use plotReducedDim() to visualise the 'pca.corrected.E8.5' layout
## coloured by the 'stage.mapped' column-wise metadata.

library(scater)
?plotReducedDim
plotReducedDim(sce, dimred = 'pca.corrected.E8.5', colour_by = "stage.mapped")

## EXERCISE 1 ##
## Create a SingleCellExperiment object “from scratch”.
## That means:
## Start from a matrix (either randomly generated or with some fake data in it)
## and add one or more columns as colData.

count_matrix <- matrix(1:20, ncol = 10)

SingleCellExperiment(assays = list(counts = count_matrix))

# see also:
# library(DropletUtils)
# read10xCounts()

## EXERCISE 2 ##
## Combine two SingleCellExperiment objects.
## The MouseGastrulationData package contains several datasets.
## Download sample 6 of the chimera experiment. Use the cbind function to combine the new data with the sce object created before.

sce_5 <- WTChimeraData(samples = 5)
sce_6 <- WTChimeraData(samples = 6)

cbind(sce_5, sce_6)

###########
# THE END #
###########

