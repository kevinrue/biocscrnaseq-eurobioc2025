################
# Bioconductor #
################

## Use install.packages() to install the package 'ggplot2' from CRAN



## Use install.packages() to install the package 'BiocManager' from CRAN



## Use BiocManager::install() to install the package 'SingleCellExperiment' from Bioconductor



## Navigate to <https://bioconductor.org/packages/release/bioc/html/scater.html>
## Then, copy paste instructions from web page to install the package 'scater'



## Use BiocManager::install() to update your packages to the latest version available



#################################
# The SingleCellExperiment class #
#################################

## Load the libraries 'SingleCellExperiment' and 'MouseGastrulationData'



## Use 'WTChimeraData()' to load the fifth sample.
## Assign the output to a new object called 'sce'



## CHALLENGE ##
## Get the data for a different sample from WTChimeraData (other than the fifth one).



## Display the names of the assays in the 'sce' object.



## Display the values in the first three rows and first three columns
## of the counts matrix of the 'sce' object.



## Display the first three rows and first four columns of the column-wise metadata
## of the 'sce' object.



## Display the first three rows and first two columns of the row-wise metadata
## of the 'sce' object.



## Use the '$' operator to add a new column named 'my_sum' to the column-wise metadata
## of the 'sce' object.
## The values in that column should be the sum for each column of the counts matrix.



## CHALLENGE ##
## Add a column of gene-wise metadata to the rowData.


