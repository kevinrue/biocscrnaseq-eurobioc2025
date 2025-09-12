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



## 