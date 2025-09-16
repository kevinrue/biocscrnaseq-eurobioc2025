###################.
# HCA Metadata ####
###################.

# Collect metadata from sample database url


# Get a view of the first 10 columns in the metadata 
# with glimpse()


# (Using the pipe operator)


# Tally the tissue types across datasets:


# Do the same for the assay types:

#...............#
# Challenge #####
#...............#
# Look through the full list of metadata column names. 
# Do any other metadata columns jump out as interesting 
# to you for your work?




###################################.
# Downloading single cell data ####
###################################.

# Identify cells meeting the following criteria:
#     - African ethnicity
#     - 10x assay
#     - lung parenchyma tissue
#     - CD4 cells



# Now we can use get_single_cell_experiment():


# Get data scaled to counts per million:


# Get data on only specific genes:


#  Return a Seurat object 


# Save your SingleCellExperiment 


#########################################.
# Exercise 1: Basic counting + piping ####
#########################################.
# Use count and arrange to get the number of cells per tissue 
# in descending order.


#########################################.
# Exercise 2: Tissue & type counting ####
#########################################.
# count() can group by multiple factors by simply adding another 
# grouping column as an additional argument. Get a tally of the 
# highest number of cell types per tissue combination. What tissue 
# has the most numerous type of cells?



################################################.
# Exercise 3: Comparing metadata categories ####
################################################.
# Spot some differences between the tissue and tissue_harmonised 
# columns. Use count to summarise.



##############################################.
# Exercise 4: Highly specific cell groups ####
##############################################.
# Now that we are a little familiar with navigating the metadata, 
# let’s obtain a SingleCellExperiment of 10X scRNA-seq counts of 
# cd8 tem lung cells for females older than 80 with COVID-19. Note: 
# Use the harmonized columns, where possible.



###########.
# THE END #
###########.

