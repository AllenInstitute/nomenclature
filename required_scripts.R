######################################################################################
## MAIN NOMENCLATURE FUNCTION

build_nomenclature_table <- function(
   dend, 
   first_label = setNames("All",labels(dend)[1]), 
   taxonomy_id = "CS1910121")
{

# `first_label` is a named vector, where 
#    the values correspond to labels (e.g., Neuron) and 
#    the names correspond to the FIRST cluster label in the tree where that label should be used
#    NOTE: this code assumes that all clusters of the same label will be in a single block in the dendrogram
#    Here is an example
#first_label <- setNames(
#    c("Neuron",            "Non-neuron"),
#    c("Inh L1 LAMP5 NDNF", "Astro L1-6 FGFR3 ETNPPL"))

################################################################
## Update the dendrogram labels

dend <- overwrite_dend_node_labels(dend)$dend
dend_start <- dend

################################################################
## Initialize the data frame

anno <- data.frame(original_label = dend %>% get_nodes_attr("label"),
                   cell_set_accession = "",
				   cell_set_label = "",
				   cell_set_alias = "",
				   cell_set_alt_alias = "")
rownames(anno) <- dend %>% get_nodes_attr("label")


################################################################
## Update cell_set_alias for cell types

anno[labels(dend),"cell_set_alias"] <- labels(dend)


################################################################
## Update cell_set_accession

# Define cluster_id
num_clusters    <- length(labels(dend))
is_leaf         <- is.element(rownames(anno),labels(dend))
anno$cluster_id <- 0
anno[labels(dend),"cluster_id"] <- 1:num_clusters
anno[!is_leaf,"cluster_id"] <- (num_clusters+1):dim(anno)[1]

# How many digits to use for cell set number (either 3 or as few as possible above 3)
cs_digits = max(3,nchar(dim(anno)[1]))

# Define cell_set_accession
anno$cell_set_accession <- paste0(taxonomy_id,substr(10^cs_digits+anno$cluster_id,2,100)) 


################################################################
## Define cell_set_label

# For cell types
first_label <- first_label[intersect(labels(dend),names(first_label))]
labs <- c(which(is.element(labels(dend),names(first_label))),num_clusters+1)
for (i in 1:length(first_label)){
  lb  <- labs[i]:(labs[i+1]-1)
  num <-   substr(10^cs_digits+anno[is_leaf,]$cluster_id[lb]-labs[i]+1,2,100)
  anno[is_leaf,]$cell_set_label[lb] <- paste(first_label[i],num)
}

# For internal nodes (cell sets)
value <- dend %>% get_nodes_attr("label")
names(value) <- dend %>% get_nodes_attr("label")
value[labels(dend)] <- anno[labels(dend),]$cell_set_label
value <- get_cell_set_designation(dend,value)
value[1] <- "All cells"
anno$cell_set_label <- as.character(value)


################################################################
## Reorganize anno table

anno <- anno[,c("cell_set_accession","original_label","cell_set_label",
                "cell_set_alias","cell_set_alt_alias")]
rownames(anno) <- NULL
anno <- anno[order(anno$cell_set_accession),]
anno$taxonomy_id = taxonomy_id

################################################################
## Update the dendrogram with this new information

dend_out <- update_dendrogram_with_nomenclature(dend,anno)

################################################################
## Return anno and dend
list(cell_set_information = anno, initial_dendrogram = dend_start, updated_dendrogram = dend_out)

}


######################################################################################
## ADDITIONAL NOMENCLATURE FUNCTIONS

update_dendrogram_with_nomenclature <- function(dend, cell_set_information, 
  current_label = "original_label", new_label = "cell_set_alias"){

  # Add all relevant information
  for (cn in colnames(cell_set_information)){
    value <- cell_set_information[,cn]
    names(value) <- cell_set_information[,current_label]
    dend <- add_attr_to_dend(dend,value,cn)
  }

  # Make the label match cell_set_alias
  value <- cell_set_information[,new_label]
  names(value) <- cell_set_information[,current_label]
  dend <- add_attr_to_dend(dend,value,"label")

  # Return dendrogram
  dend
}


cell_set_mapping_from_dendrogram <- function(dend,samples,cell_id,
  mapping=NULL,continue=TRUE){  # DO NOT SET THESE LAST TWO VARIABLES

  # Add the mapping for the current cell set
  cell_set_id <- (dend %>% get_nodes_attr("cell_set_accession"))[1]
  kp <- is.element(cell_id,dend %>% get_leaves_attr("cell_set_accession"))
  if(sum(kp)>0){
    newMapping  <- cbind(samples[kp],cell_set_id)
    mapping     <- rbind(mapping,newMapping)
  }
  
  # Recursively work through the dendrogram adding all mappings
  if(length(dend)>1) 
    for (i in 1:length(dend)){
       mapping <- cell_set_mapping_from_dendrogram(dend[[i]],samples,cell_id,mapping,continue)
  } else if(continue) {
    continue <- FALSE
    mapping  <- cell_set_mapping_from_dendrogram(dend,samples,cell_id,mapping,continue)
  }
  
  # Return results
  colnames(mapping) <- c("sample_name","cell_set_accession")
  mapping
}


######################################################################################
## Support functions

get_dendrogram_value <- function(dend,value, sep=" "){
  labs <- as.character(value[labels(dend)])
  clas <- as.character(unclass(sapply(labs, function(x) strsplit(x,sep)[[1]][1])))
  nums <- as.character(unclass(sapply(labs, function(x) strsplit(x,sep)[[1]][2])))
  if(length(unique(clas))>1){
   val <- paste(unique(clas),collapse = "/")
  } else {
   val <- paste(clas[1], paste(unique(range(nums)),collapse="-"))
  }
  return(val)
 }

 
get_cell_set_designation <- function(dend, value, sep=" ") {
 if (length(dend)>1) {
   for(i in 1:length(dend)){
     value[attr(dend[[i]],"label")] <- get_dendrogram_value(dend[[i]],value)
     value = get_cell_set_designation(dend[[i]], value=value)
   }
 }
 return(value)
}


add_attr_to_dend <- function(dend, value, attribute="label") {
  # Value must be of length nnodes(dend) and be named as such
  if(!is.na(value[attr(dend,"label")]))
    attr(dend, attribute) <- value[attr(dend,"label")]
  if (length(dend)>1) {
    for(i in 1:length(dend)){
      dend[[i]]=add_attr_to_dend(dend[[i]], value=value, attribute)
    }
  }
  return(dend)
}



plot_dend <- function (dend, dendro_data = NULL, node_size = 1, r = c(-0.1, 1)) 
{
    require(dendextend)
    require(ggplot2)
    if (is.null(dendro_data)) {
        dendro_data = as.ggdend(dend)
        dendro_data$nodes$label = get_nodes_attr(dend, "label")
        dendro_data$nodes = dendro_data$nodes[is.na(dendro_data$nodes$leaf), 
            ]
    }
    node_data = dendro_data$nodes
    label_data <- dendro_data$labels
    segment_data <- dendro_data$segments
    if (is.null(node_data$node_color)) {
        node_data$node_color = "black"
    }
    ggplot() + geom_text(data = node_data, aes(x = x, y = y, 
        label = label, color = node_color), size = node_size, 
        vjust = 1) + geom_segment(data = segment_data, aes(x = x, 
        xend = xend, y = y, yend = yend), color = "gray50") + 
        geom_text(data = label_data, aes(x = x, y = -0.01, label = label, 
            color = col), size = node_size, angle = 90, hjust = 1) + 
        scale_color_identity() + theme_dendro() + scale_y_continuous(limits = r)
}


overwrite_dend_node_labels <- function (dend, n = 1, lab = labels(dend)) 
{
    if (!is.element(attr(dend, "label"),lab)) {
        attr(dend, "label") = paste0("n", n)
        n = n + 1
    }
    if (length(dend) > 1) {
        for (i in 1:length(dend)) {
            tmp = overwrite_dend_node_labels(dend[[i]], n, lab)
            dend[[i]] = tmp[[1]]
            n = tmp[[2]]
        }
    }
    return(list(dend = dend, n))
}


dend_to_list <- function(dend, omit_names = c("markers","markers.byCl","class")) {

  # NOTE: If this function crashes, the "omit_names" variable may need to be updated to
  #       exclude additional variables in the dendrogram object
  node_attributes <- as.data.frame(attributes(dend)[!is.element(names(attributes(dend)),omit_names)])
  node_attributes <- unique(node_attributes[,names(node_attributes) != "names"])
  if("leaf" %in% names(node_attributes)) {
    return(list(leaf_attributes = node_attributes))
  } else {
    y <- dend
    attributes(y) <- NULL
    class(y) <- "list"
    children <- y
    
    dend <- list(node_attributes = node_attributes,
                 children = children)
    
    if(length(dend$children) > 1) {
      for(i in 1:length(dend$children)) {
        dend$children[[i]] <- dend_to_list(dend$children[[i]])
      }
    }
    return(dend)
  }
  
}

