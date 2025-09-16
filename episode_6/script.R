###################.
# HCA Metadata ####
###################.

# Install CuratedAtlasQueryR from Bioconductor
BiocManager::install("CuratedAtlasQueryR")

# Load libraries CuratedAtlasQueryR and dplyr
library(CuratedAtlasQueryR)
library(dplyr)

# Collect metadata from sample database url
metadata <- get_metadata(remote_url = CuratedAtlasQueryR::SAMPLE_DATABASE_URL)

# Get a view of the first 10 columns in the metadata 
# with glimpse()

metadata |> 
  select(1:10) |> 
  glimpse()

# (Using the pipe operator)
mtcars |> 
  filter(cyl != 4) |> 
  summarise(avg_disp = mean(disp),
            .by = cyl)
summarise(filter(mtcars, cyl != 4), mean_disp = mean(disp),
          .by = cyl)

# Tally the tissue types across datasets:
metadata |> 
  distinct(tissue, dataset_id) |> 
  count(tissue) |> 
  arrange(-n)

# Do the same for the assay types:
metadata |> 
  distinct(assay, dataset_id) |> 
  count(assay)

#...............#
# Challenge #####
#...............#
# Look through the full list of metadata column names. 
# Do any other metadata columns jump out as interesting 
# to you for your work?
colnames(metadata)



###################################.
# Downloading single cell data ####
###################################.

# Identify cells meeting the following criteria:
#     - African ethnicity
#     - 10x assay
#     - lung parenchyma tissue
#     - CD4 cells

sample_subset <- metadata |> 
  filter(
    ethnicity == "African" &
    grepl("10x", assay) &
    tissue == "lung parenchyma" &
    grepl("CD4", cell_type)
      
  )

# Now we can use get_single_cell_experiment():
single_cell_counts <- sample_subset |> 
  get_single_cell_experiment()

# Get data scaled to counts per million:
sample_subset |> 
  get_single_cell_experiment(assays = "cpm")


# Get data on only specific genes:
sample_subset |> 
  get_single_cell_experiment(assays = "cpm",
                             features = c("PUM1", "MAFB"))


#  Return a Seurat object 
sample_subset |> 
  get_seurat()


# Save your SingleCellExperiment 
single_cell_counts |> 
  HDF5Array::saveHDF5SummarizedExperiment("single_cell_counts")

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

