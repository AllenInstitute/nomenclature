
################################################################
## Variable definitions

# `taxonomy_id` is the name of the taxonomy in the format: <CS><YYMMDD><T>, where:
#    CS stands for "cell set" and CT stands for "cell type"
#    YYMMDD represents a 6 digit date format (Y=year, M=month, D=day)
#    T is a 1-digit taxonomy counter, which allows up to 10 taxonomies on the same date

# IMPORTANT NOTE: To keep taxonomy IDs unique, please select a taxonomy_id NOT IN THIS TABLE:
#    https://docs.google.com/spreadsheets/d/10gYNyOhc0YHOYKjgsvLfumf65CqLeiVDCE3Wrxz4Txo/edit?usp=sharing 
#    and add your taxonomy_ID.  Strategies for better tracking of taxonomy IDs are currently 
#    under consideration. 

taxonomy_id <- "CS1910121"

# `first_label` is a named vector, where 
#    the values correspond to labels (e.g., Neuron) and 
#    the names correspond to the FIRST cluster label in the tree where that label should be used
#    NOTE: this code assumes that all clusters of the same label will be in a single block in the dendrogram
first_label <- setNames(
    c("Neuron",            "Non-neuron"),
    c("Inh L1 LAMP5 NDNF", "Astro L1-6 FGFR3 ETNPPL"))
	
	
################################################################
## Read in dendrogram

# We have provided a dendrogram called "dend" as an example.  Any dendrogram of cell types in the 
# "dendrogram" R format, where the desired cluster aliases are in the "labels" field of the dendrogram
# will work for this code.  Other formats might work and will try to be forced into dendrogram format.
load("dend_humanMTG.RData")

# Attempt to format dendrogram if the input is in a different format
dend <- as.dendrogram(dend)

# Plot the dendrogram to test
# plot(dend,main="You should see your desired cell type names on the base of this plot")
	
	
################################################################
## Load required libraries and functions

library(dplyr)
library(dendextend)
library(ggplot2)
library(jsonlite)
library(data.table)
source("required_scripts.R")  # Additional required files
options(stringsAsFactors = FALSE)


################################################################
## Assign the nomenclature!

# This is done as a single line of code, but can also be run line by line if desired --
#   in this case, see the "requiredScripts.R" file

nomenclature_information <- build_nomenclature_table(dend, first_label, taxonomy_id)


################################################################
## Save the initial dendrogram and nomenclature table

pdf("initial_dendrogram.pdf",height=8,width=15)
plot_dend(nomenclature_information$initial_dendrogram, node_size=3)
dev.off()

write.csv(nomenclature_information$cell_set_information,"nomenclature_table_initial.csv",row.names=FALSE)


################################################################
##   TEXT DISCUSSING MANUAL ANNOTATION OF ALIAS FIELDS HERE   ##
################################################################

print("Add alt aliases to the table and then rename 'nomenclature_table.csv'.")


################################################################
## Read in the updated nomenclature

updated_nomenclature <- read.csv("nomenclature_table_humanMTG.csv")


################################################################
## Update the dendrogram and plot the results

updated_dendrogram <- update_dendrogram_with_nomenclature(nomenclature_information$initial_dendrogram,updated_nomenclature)

pdf("updated_dendrogram.pdf",height=8,width=15)
plot_dend(updated_dendrogram, node_size=3)
dev.off()


################################################################
## Save the dendrogram in various formats

# Save as an R data object
save(updated_dendrogram, file="updated_dendrogram.RData")

# Convert to a list
# NOTE: Only some features of dendrogram can be converted to a list.  If this function 
#       crashes, the "omit_names" variable may need to be updated
dend_list <- dend_to_list(updated_dendrogram, omit_names = c("markers","markers.byCl","class"))

# Save as a json file
dend_JSON <- toJSON(dend_list, complex = "list", pretty = TRUE)
out <- file("dend.json", open = "w")
writeLines(dend_JSON, out)
close(out)


################################################################
## Read in cell meta-data and output a cell mapping table

# Read in metadata and collect correct columns for sample name and cell set accession id
metadata  <- read.csv("cell_metadata.csv")
samples   <- metadata$sample_name

# OPTION 1: COLUMN FOR ACCESSION ID ALREADY EXISTS
# cell_id <- metadata$cell_type_accession_label

# OPTION 2: NEED TO GENERATE COLUMN FROM DENDROGRAM LABELS
label_col <- "cluster_label"  # Column name with dendrogram labels
cell_id   <- updated_nomenclature[match(metadata[,label_col],updated_nomenclature$cell_set_alias),"cell_set_accession"]
cell_id[is.na(cell_id)] = "none"

# Perform the mapping
mapping   <- cell_set_mapping_from_dendrogram(updated_dendrogram,samples,cell_id)

# Output a csv of the mapping
fwrite(mapping,"cell_mapping.csv")
