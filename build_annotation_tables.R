
################################################################
## Variable definitions

# `taxonomy_name` is the name of the taxonomy in the format: <CS><YYMMDD><T>, where:
#    CS stands for "cell set" and CT stands for "cell type"
#    YYMMDD represents a 6 digit date format (Y=year, M=month, D=day)
#    T is a 1-digit taxonomy counter, which allows up to 10 taxonomies on the same date
taxonomy_name <- "CS1910121"

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
source("required_scripts.R")  # Additional required files
options(stringsAsFactors = FALSE)


################################################################
## Assign the nomenclature!

# This is done as a single line of code, but can also be run line by line if desired --
#   in this case, see the "requiredScripts.R" file

nomenclature_information <- build_nomenclature_table(dend)


################################################################
## Save the initial dendrogram and nomenclature table

pdf("initial_dendrogram.pdf",height=8,width=15)
plot_dend(nomenclature_information$initial_dendrogram, node_size=3)
dev.off()

write.csv(nomenclature_information$cell_set_information,"nomenclature_table.csv",row.names=FALSE)


################################################################
##   TEXT DISCUSSING MANUAL ANNOTATION OF ALIAS FIELDS HERE   ##
################################################################


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
