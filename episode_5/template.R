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



## Use pairwiseRand() to measure the similarity between the two sets
## of cluster assignments.
## Assign the result to a new object named 'rand'.
## Hint: set the argument 'mode' to 'index' to return the Rand index.



## Copy paste the code from the online workshop materials
## to demonstrate a case where approximate neighbour search deteriorates.
## <https://carpentries-incubator.github.io/bioc-scrnaseq/large_data.html#fast-approximations>
## Tip: Set the random seed to ensure reproducible results! (e.g., 1000)



################################
# Singular value decomposition #
################################

## Load the following libraries
## - scater
## - BiocSingular


## Use runPCA() with RandomParam() to compute 20 principal components for 'sce'.
## Tip: Set the random seed to ensure reproducible results! (e.g., 101000)
## Assign the result to a new object named 'r.out'.



## Combine str() and reduceDim() to display information about the PCA results computed above.



## Repeat the last two steps except
## - Use IrlbaParam() instead of RandomParam()
## Tip: Set the random seed to ensure reproducible results! (e.g., 101001)



## CHALLENGE ##
## The uncertainty from approximation error is sometimes aggravating.
## “Why can’t my computer just give me the right answer?”
## One way to alleviate this feeling is to quantify the approximation error on a small test set like the sce we have here.
## Using the ExactParam() class, visualize the error in PC1 coordinates compared to the RSVD results.



################################################################
# Interoperability with popular single-cell analysis ecosytems #
################################################################

## Load the library 'Seurat'



## Use BiocManager::install() to install the R package in the GitHub repository 'satijalab/seurat-data'



## Load the library 'SeuratData'



## Use InstallData() to install the dataset 'pbmc3k'



## Use LoadData() to load the 'pbmc3k.final' version of the dataset 'pbmc3k'.
## Assign the result to a new object named 'pbmc'.



## Use UpdateSeuratObject() to update the object for compatibility
## with the version of Seurat active in your R session.
## Reassign the result to 'pbmc'.



## Use as.SingleCellExperiment() to convert 'pbmc' from a Seurat object
## to a Bioconductor SingleCellExperiment object.
## Assign the result to a new object named 'pbmc.sce'.



## Use WTChimeraData() to load the 'processed' data for the fifth sample.
## Assign the output to a new object called 'sce'.



## Use as.matrix() to convert the first assay (the counts matrix)
## to a regular matrix.



## Use logNormCounts() to compute log-normalised counts.
## Reassign the result to 'sce'.



## Use as.Seurat() to convert 'sce' from a Bioconductor SingleCellExperiment object
## to a Seurat object.
## Assign the result to a new object named 'sobj'.



##########
# Scanpy #
##########

## Load the library 'zellkonverter'



## Use system.file() to get the full file path to the file 'krumsiek11.h5ad'
## supplied in the 'zellkonverter' package.



## Use readH5AD() to import the file above as a SingleCellExperiment object.



## Use tempfile() to generate a temporary file name with the suffix '.h5ad'.



## Use writeH5AD() to export 'sce' to the temporary file above.


## EXERCISE 1 - OUT OF MEMORY REPRESENTATION ##
## Write the counts matrix of the wild-type chimera mouse gastrulation dataset to an HDF5 file.
## Create another counts matrix that reads the data from the HDF5 file.
## Compare memory usage of holding the entire matrix in memory as opposed to holding the data out of memory.



## EXERCISE 2 - PARALLELISATION ##
## Perform a PCA analysis of the wild-type chimera mouse gastrulation dataset
## using a multicore backend for parallel computation.
## Compare the runtime of performing the PCA either in serial execution mode,
## in multicore execution mode with 2 workers, and in multicore execution mode with 3 workers.



###########
# THE END #
###########
