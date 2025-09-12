#################################
# Out of memory representations #
#################################

## Load the library 'TENxBrainData'



## Use TENxBrainData20k() to load a SingleCellExperiment for this episode.
## Assign the output to a new object called 'sce.brain'.



## Display 'sce.brain' and briefly inspect it.



## Use counts() to display the counts matrix of 'sce.brain'.
## Pay particular attention to the class of the object.



## Use object.size() to display the memory footprint of the counts matrix of 'sce.brain'.



## Combine file.info() and path() to display the size of the file containing the counts data on disk.



## Use counts() to extract the counts matrix of 'sce.brain' again,
## but this time assign the result to a new object named 'tmp'.



## Use log2() to apply log2-transformation to the counts after adding 1 to avoid taking the log of 0.
## Reassign the result to 'tmp'.



## Display 'tmp' and pay particular attention to its class.



## Load the library 'scater'



## Use grepl() and the pattern '^mt-' to identify  mitochondrial features.
## Assign the result to a new object named 'is.mito'.



## Use perCellQCMetrics()  to calculate QC metrics for 'sce.brain',
## giving the subset of mitochondrial genes to compute additional QC metrics based on that subset.
## Assign the result to a new object named 'qcstats'.


###################
# Parallelization #
###################

## Load the library 'BiocParallel'



## Use MulticoreParam() to create a parameter object for parallelisation.
## Give it a single 'worker'.
## Assign the result to a new object named 'param'.



## Use bpalpply() to apply sqrt() on the vector of numbers c(4, 9, 16, 25),
## using the parallelisation object created above.



## Load the following libraries
## - MouseGastrulationData
## - scran



## Use WTChimeraData() to load the 'processed' data for the fifth sample.
## Assign the output to a new object called 'sce'.



## Use logNormCounts() to compute log-normalised counts for 'sce'.
## Reassign the result to 'sce'.



## Use modelGeneVar() to model gene variance for 'sce',
## using MulticoreParam() to parallelise calculations across two workers.
## Assign the result to a new object named 'dec.mc'.



## Same as above except:
## Try using SnowParam to parallelise calculations across two workers.
## Assign the result to a new object named 'dec.snow'.



## CHALLENGE ##
## How do you turn on progress bars with parallel processing?




#######################
# Fast approximations #
#######################

## Load the library 'bluster'



## Use runPCA() to compute principal components for 'sce'.



## Use clusterCells() to assign cluster labels for each cell in 'sce',
## using the PCA results computed above,
## and using a NNGraphParam() parameter object specifying 'louvain' clustering.
## Assign the cluster labels to the 'sce' object using colLabels().



## Load the following libraries
## - scran
## - BiocNeighbors



## Use clusterCells() as above except
## use an AnnoyParam() parameter object to use approximate nearest neighbour search.
## Assign the result to a new object named 'clusters'.



## Use table() to compare the cluster assignments
## - using exact nearest neighbour search, stored in colLabels()
## - using approximate nearest neighbour search, stored in 'clusters'


