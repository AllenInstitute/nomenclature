######################################################################################
## MAIN NOMENCLATURE FUNCTION

build_nomenclature_table <- function(
   dend, 
   first_label  = setNames("All",labels(dend)[1]), 
   taxonomy_id  = "CCN201910121",
   taxonomy_author = "Unspecified", 
   taxonomy_citation = "", 
   structure    = "neocortex",
   ontology_tag = "UBERON:0001950")
{

# dend is the dendrogram variable and must be included (all other variables can be inferred)
# `first_label` is a named vector, where 
#    the values correspond to labels (e.g., Neuron) and 
#    the names correspond to the FIRST cluster label in the tree where that label should be used -OR-
#       the index correspoding the desired cluster label
#    NOTE: this code assumes that all clusters of the same label will be in a single block in the dendrogram
#    Here is an example
# first_label <- setNames(
#    c("Neuron",            "Non-neuron"),
#    c("Inh L1 LAMP5 NDNF", "Astro L1-6 FGFR3 ETNPPL"))

#`taxonomy_id` is the name of the taxonomy in the format: <CCN><YYYYMMDD><T>, where:
# * CCN stands for "Common Cell type Nomenclature"
# * YYYYMMDD represents an 8 digit date format (Y=year, M=month, D=day)
# * T is a 1-digit taxonomy counter, which allows up to 10 taxonomies on the same date  
  
#`taxonomy_author` is the name of a point person for this taxonomy (e.g., someone who is 
#   responsible for it's content).  This could either be the person who built the taxonomy, 
#   the person who uploaded the data, or the first or corresponding author on a relevant manuscript.  
  
#`taxonomy_citation` is a the citation or permanent data identifier corresponding to the taxonomy 
#   (or "" if there is no associated citation).  Ideally the DOI for the publication will be used, 
#   or alternatively some other permanent link.  
  
# `structure` is the brain or body structure from which the taxonomy was derived
# `ontology_tag` is the identifier of the above structure from a specific ontology (or "none")
# Note that `structure` and `ontology_tag` can be manually defined separately for each cell set at a later step.

  
################################################################
## Update the dendrogram labels

dend <- overwrite_dend_node_labels(dend)$dend
dend_start <- dend


################################################################
## Initialize the data frame

anno <- data.frame(original_label = dend %>% get_nodes_attr("label"),
                   cell_set_accession = "",
                   cell_set_label = "",
                   cell_set_preferred_alias = "",
                   cell_set_aligned_alias = "",
                   cell_set_additional_aliases = ""
                   )
rownames(anno) <- dend %>% get_nodes_attr("label")


################################################################
## Update cell_set_preferred_alias for cell types

anno[labels(dend),"cell_set_preferred_alias"] <- labels(dend)


################################################################
## Update cell_set_accession

# Define cluster_id
num_clusters    <- length(labels(dend))
is_leaf         <- is.element(rownames(anno),labels(dend))
anno$cluster_id <- 0
#anno[labels(dend),"cluster_id"] <- 1:num_clusters  # Replaced with next line, which avoids error
anno$cluster_id[rownames(anno) %in% labels(dend] <- 1:num_clusters   #replaces cluster_id=0 with sequential id
anno[!is_leaf,"cluster_id"] <- (num_clusters+1):dim(anno)[1]
cluster_id      <- anno$cluster_id

# How many digits to use for cell set number (as few as possible)
cs_digits = max(1,nchar(dim(anno)[1]))

# Define cell_set_accession
anno$cell_set_accession <- paste0(taxonomy_id,"_",anno$cluster_id) #substr(10^cs_digits+anno$cluster_id,2,100)) 
anno$cell_set_accession <- gsub("CCN","CS",anno$cell_set_accession)


################################################################
## Define cell_set_label
# NOTE: cell_set_label is deprecated in the CCN, but is still required for this script to run properly

# First convert from dendrogram index to label if needed
if(!is.na(suppressWarnings(as.numeric(names(first_label))))){
  index <- as.numeric(names(first_label))
  first_label <- first_label[as.character(intersect(1:length(labels(dend)),index))]
  index <- as.numeric(names(first_label))
  names(first_label) <- labels(dend)[index]
}

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
## Append the cell_set_structure

anno$cell_set_structure    <- structure
anno$cell_set_ontology_tag <- ontology_tag


################################################################
## Reorganize anno table

anno <- anno[,c("cell_set_accession","original_label","cell_set_label",
                "cell_set_preferred_alias","cell_set_aligned_alias","cell_set_additional_aliases",
                "cell_set_structure","cell_set_ontology_tag")]
rownames(anno) <- NULL
anno <- anno[order(cluster_id),]

anno$cell_set_alias_assignee <- paste0(taxonomy_author)
anno$cell_set_alias_citation <- paste0(taxonomy_citation)
anno$taxonomy_id             <- taxonomy_id

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
  current_label = "original_label", new_label = "cell_set_preferred_alias"){

  # Add all relevant information
  for (cn in colnames(cell_set_information)){
    value <- cell_set_information[,cn]
    names(value) <- cell_set_information[,current_label]
    dend <- add_attr_to_dend(dend,value,cn)
  }

  # Make the label match cell_set_preferred_alias
  value <- cell_set_information[,new_label]
  names(value) <- cell_set_information[,current_label]
  dend <- add_attr_to_dend(dend,value,"label")

  # Return dendrogram
  dend
}


cell_assignment_from_dendrogram <- function(dend,samples,cell_id,
  mapping=data.frame(sample_name=samples),continue=TRUE){  # DO NOT SET THESE LAST TWO VARIABLES
  
  # Add the mapping for the current cell set
  cell_set_id <- (dend %>% get_nodes_attr("cell_set_accession"))[1]
  kp <- is.element(cell_id,dend %>% get_leaves_attr("cell_set_accession"))
  if(sum(kp)>0){
    mapping[,cell_set_id] <- ifelse(kp,1,0)
  }
  
  # Recursively work through the dendrogram adding all mappings
  if(length(dend)>1) 
    for (i in 1:length(dend)){
       mapping <- cell_assignment_from_dendrogram(dend[[i]],samples,cell_id,mapping,continue)
  } else if(continue) {
    continue <- FALSE
    mapping  <- cell_assignment_from_dendrogram(dend,samples,cell_id,mapping,continue)
  }
  
  # Return results
  mapping
}



cell_assignment_from_groups_of_cell_types <- function(updated_nomenclature,cell_id,mapping){

  ## Determine the relevant cell sets to annotate 
  used_ids      <- intersect(updated_nomenclature$cell_set_accession,colnames(mapping))
  missed_ids    <- setdiff(updated_nomenclature$cell_set_accession,colnames(mapping))
  missed_labels <- updated_nomenclature$cell_set_label[match(missed_ids,updated_nomenclature$cell_set_accession)]
  used_labels   <- setdiff(updated_nomenclature$cell_set_label,missed_labels)
  
  used_class       <- as.character(lapply(used_labels,function(x) strsplit(x," ")[[1]][1]))
  missed_class     <- as.character(lapply(missed_labels,function(x) strsplit(x," ")[[1]][1]))
  cell_type_labels <- missed_labels[is.element(missed_class,used_class)]

  ## Find the corresponding cell types for those cell sets
  labsL <- list()
  cs_digits <- max(1,nchar(dim(updated_nomenclature)[1]))
  if (length(missed_labels)>0){
    for (i in 1:length(missed_labels)){
      m <- eval(parse(text=paste("c(",gsub("-",":",gsub(missed_class[i],"",missed_labels[i])),")")))
      m <- substr(10^cs_digits+m,2,100)
      labsL[[missed_labels[i]]] <- paste(missed_class[i],m)
    }
  }
  
  ## Now annotate them!
  print("Cell sets added to table:")
  for (cl in cell_type_labels){
    ids <- updated_nomenclature$cell_set_accession[match(labsL[[cl]],updated_nomenclature$cell_set_label)]
    kp  <- is.element(cell_id,ids)
    if(sum(kp)>0){
      cell_set_id <- updated_nomenclature$cell_set_accession[match(cl,updated_nomenclature$cell_set_label)]
      mapping[,cell_set_id] <- ifelse(kp,1,0)
      print(cell_set_id)
    }
  }
  
  # Return results
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

