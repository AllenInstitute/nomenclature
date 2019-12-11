This script will allow easy application of the [Allen Institute nomenclature convention](http://confluence.corp.alleninstitute.org/pages/viewpage.action?pageId=69045259) to an existing hierarchically-structured taxonomy of cell types.  It takes an R `dendrogram` object, along with a couple of critical pieces of meta-data, as input, and outputs a table with nomenclature and a labeled dendrogram to allow for manual annotation of addition cell type and cell set aliases.  Finally, after manual annotation, the resulting dendrogram can be saved in multiple formats.  

# Workspace set-up

Prior to running the scripts, the follow steps need to be performed.
1. Install R on your computer
2. Install these libraries: `dplyr`, `dendextend`, and `ggplot2`.  
3. *(Optional)* Install the `jsonlite` library if you want to save the final dendrogram in json format.
4. Download `required_scripts.R` to your working directory.
5. *(Optional)* Download `dend_humanMTG.RData` and `nomenclature_table_humanMTG.csv` if you want to use our example from human MTG and 5 additional cortical areas (see our [Allen Institute Transcriptomics Explorer](http://celltypes.brain-map.org/rnaseq/human)).
6. *(Optional)* Download `build_annotation_tables.R`, which is an r script of this page.


# Build the nomenclature

### Load required libraries 

``` r
library(dplyr)
library(dendextend)
library(ggplot2)
library(jsonlite)  # optional
```

### Load the accessory scripts.

```
source("required_scripts.R")  # Additional required files
options(stringsAsFactors = FALSE)
```

### Define the variables `taxonomy_name` and `first_label`.  

`taxonomy_name` is the name of the taxonomy in the format: <CS><YYMMDD><T>, where:
* CS stands for "cell set" and CT stands for "cell type"
* YYMMDD represents a 6 digit date format (Y=year, M=month, D=day)
* T is a 1-digit taxonomy counter, which allows up to 10 taxonomies on the same date
`first_label` is a named vector, where 
* the values correspond to labels (e.g., Neuron) and 
* the names correspond to the FIRST cluster label in the tree where that label should be used
* *NOTE: this code assumes that all clusters of the same label will be in a single block in the dendrogram*

``` r
taxonomy_name <- "CS1910121"

first_label <- setNames(
    c("Neuron",            "Non-neuron"),
    c("Inh L1 LAMP5 NDNF", "Astro L1-6 FGFR3 ETNPPL"))
``` 	
	
### Read in dendrogram

We have provided a dendrogram called "dend" as an example.  Any dendrogram of cell types in the "dendrogram" R format, where the desired cluster aliases are in the "labels" field of the dendrogram will work for this code.  Other formats might work and will try to be forced into dendrogram format.

``` r
# REPLACE THIS LINE OF CODE WITH CODE TO READ IN YOUR DENDROGRAM, IF NEEDED
load("dend_humanMTG.RData")

# Attempt to format dendrogram if the input is in a different format
dend <- as.dendrogram(dend)

# Uncomment this line if you'd like to plot the dendrogram to test it's format
# plot(dend,main="You should see your desired cell type names on the base of this plot")
```
	
### Assign the nomenclature!

The entire script for assigning all nomenclature is done as a single line of code.  If you'd prefer to run it line by line (for example if your data is in slightly different format), see the `build_nomenclature_table` function in the `required_scripts.R` file.  This function has been reasonably-well commented to try and explain how each section works.  

``` r
nomenclature_information <- build_nomenclature_table(dend)
```

### Save the initial dendrogram and nomenclature table

These are the final files, prior to any manual annotations described below.  `nomenclature_table.csv` contains four annotation files-- --which are described in detail [on our website, here](http://confluence.corp.alleninstitute.org/pages/viewpage.action?pageId=69045259), as well as an `original_label` column, which can be used to match the rows of this table with the nodes in `initial_dendrogram.pdf`.  *Note: for convenience, we are using the term "cell_set" to refer to both "cell_type" and "cell_set" in this table.  Leaf nodes should be considered cell types as described in the nomenclature documentation. 

``` r
pdf("initial_dendrogram.pdf",height=8,width=15)
plot_dend(nomenclature_information$initial_dendrogram, node_size=3)
dev.off()

write.csv(nomenclature_information$cell_set_information,"nomenclature_table.csv",row.names=FALSE)
```

# Manual annotation of aliases

By default all cell types ("leaf" nodes) have exactly one alias and all cell sets (internal nodes) do not have any aliases.  This step is where you can update the `nomenclature_table.csv` file to add these aliases.  Unless you note an error, **no not update this file except to fill blanks in the `cell_set_alias` and `cell_set_alt_alias` collumns.**  An example file for our human MTG data set (`nomenclature_table_humanMTG.csv`) is provided as an example.  More detail about aliases is provided at the link above, but in short:

* cell_set_alias - nodes representing a useful collection of cell types can be manually tagged with an alias by matching the node label shown on the plotted dendrogram with the `original_label` column.  For example, a node containing all of the Pvalb+ GABA-ergic neuron types might be labeled "Pvalb".  Cell types already have labels, which should be left alone.
* cell_type_alt_alias - this slot represents any additional names that you'd like to save for a given cell type.  These could include **common usage names** (such as "Chanelier cells") or other aliases that will allow matching between taxonomies.  Typically this column is used only for cell types and not internal nodes.

If desired, other columns can also be added to this table, and those columns will be appended to the dendrogram object in the code below as well.  Once this file has been completed, save the result as a csv file, and continue with the code below, using that file name as input.  

# Update and save the dendrogram

### Read in the updated nomenclature

``` r
updated_nomenclature <- read.csv("nomenclature_table_humanMTG.csv")
```

### Update the dendrogram and plot the results

This code will take the information from the table above and add it to the initial dendrogram object.  When plotted the only visible difference will be that the new cell set alias names (if any) will show up to replace the n## labels from the initial plot.  

``` r
updated_dendrogram <- update_dendrogram_with_nomenclature(nomenclature_information$initial_dendrogram,updated_nomenclature)

pdf("updated_dendrogram.pdf",height=8,width=15)
plot_dend(updated_dendrogram, node_size=3)
dev.off()
```

### Save the dendrogram in various formats

Plots only show a small fraction of the data available in these dendrogram objects; to see the rest the dendrogram needs to be saved.  We fine both the R "dendrogram" format, as well as the "json" format useful for different applications at the Allen Institute and code for saving data as both are presented below.

``` r 
# Save as an R data object
save(updated_dendrogram, file="updated_dendrogram.RData")
```

(This section can be skipped if json format is not needed.)

``` r
# Convert to a list
# NOTE: Only some features of dendrogram can be converted to a list.  If this function 
#       crashes, the "omit_names" variable may need to be updated
dend_list <- dend_to_list(updated_dendrogram, omit_names = c("markers","markers.byCl","class"))

# Save as a json file
dend_JSON <- toJSON(dend_list, complex = "list", pretty = TRUE)
out <- file("dend.json", open = "w")
writeLines(dend_JSON, out)
close(out)
```
