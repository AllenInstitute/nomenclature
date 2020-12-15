make_preferred_alias_unique <- function(nomenclature){
  # This function converts the preferred_alias vector to a unique vector
  #  combining preferred alias and cell set label
  preferred_alias_in <- nomenclature$cell_set_preferred_alias
  lus <- function(x) length(unique(setdiff(x,"")))
  children    <- split_column(nomenclature,"child_cell_set_accessions")
  numChildren <- apply(children,1,lus)
  preferred_alias_in[numChildren==max(numChildren)] = "All cells"
  label <- gsub(paste0(gsub("CCN","CS",nomenclature$taxonomy_id[1]),"_"),"", 
                nomenclature$cell_set_accession)
  
  preferredAlias <- NULL
  for (i in 1:length(preferred_alias_in)) {
    y <- unique(children[i,])
    isDescendent <- apply(children,1,function(x,y) lus(y)==lus(intersect(x,y)),y)
	isDescendent[i] = TRUE
	isDescendent[preferred_alias_in==""] = FALSE
	isDescendent[numChildren==max(numChildren)] = TRUE
	kp <- which(isDescendent&(numChildren==min(numChildren[isDescendent])))
	if(length(y)==1) kp = i  # To account for no children case
	sa <- ifelse(i==kp[1],""," subset")
	pa <- paste0(preferred_alias_in[kp[1]],sa," (",label[i],")")
    preferredAlias <- c(preferredAlias,pa)
  }
  preferredAlias
}


make_sample_csv <- function(mapping, uri, dataset, accessionId, 
                            preferredAlias, ExternalId=NULL,
							useIds = accessionId){
  # This function returns a table of the correct format for sample.csv
  out <- NULL
  mapping <- mapping[,accessionId]
  for (i in 1:length(accessionId)) if(is.element(accessionId[i],useIds)){
    kp  <- mapping[,i]>0
    tmp <- data.frame(
	  Name = rownames(mapping)[kp],
	  Uri = uri[kp],
	  DataSet = dataset,
	  CellSetPreferredAlias = preferredAlias[i],
	  CellSetProbability = mapping[kp,i]
	)
	if(!is.null(ExternalId)) 
	  tmp$ExternalId = ExternalId[kp]
	out <- rbind(out,tmp)
  }
  if(!is.element(colnames(out),"ExternalId"))
    out$ExternalId = 1:dim(out)[1]
  rownames(out) = NULL
  cn = c("ExternalId","Name","Uri","DataSet","CellSetPreferredAlias",
         "CellSetProbability")
  out[,cn]
}


split_column <- function(nomenclature, column_name, split_char = "\\|", expand=FALSE){
  # This function splits a column in the nomenclature table by a defined character
  # --- Split out the values in the correct column
  splitColumn = strsplit(nomenclature[,column_name],split_char)
  len <- -1
  for (i in 1:length(splitColumn)) len=max(len,length(splitColumn[[i]]))
  # --- Convert it to a matrix and fill in missing values
  out <- matrix("",nrow=length(splitColumn),ncol=len)
  for (i in 1:length(splitColumn)){
    tmp = splitColumn[[i]]
	if(length(tmp)>1) for (j in 2:length(tmp)) if(tmp[j]=="")
	  tmp[j] = tmp[j-1]
	if(length(tmp)>0){
      out[i,1:length(tmp)] = tmp
	  if(expand&(length(tmp)<len)){
	    out[i,(length(tmp)+1):len] = tmp[length(tmp)]
	  }
	}
  }
  # --- Return results
  out
}


find_cell_set_family <- function(preferredAlias,nomenclature){
  lus <- function(x) length(unique(setdiff(x,"")))
  children    <- split_column(nomenclature,"child_cell_set_accessions")
  numChildren <- apply(children,1,lus)
  familyPairs <- NULL
  for (i in 1:length(preferredAlias)) if(numChildren[i]!=max(numChildren)){
    y <- unique(children[i,])
	if(length(y)==1) y = nomenclature$cell_set_accession[i]
    isDescendent <- apply(children,1,function(x,y) lus(y)==lus(intersect(x,y)),y)
	isDescendent[i] = FALSE
	kp <- which(isDescendent&(numChildren==min(numChildren[isDescendent])))
	familyPairs <- rbind(familyPairs,c(preferredAlias[kp[1]],preferredAlias[i]))
  }
  colnames(familyPairs) = c("ParentAlias","ChildAlias")
  familyPairs
}


retrieve_alias_information <- function(preferredAlias,nomenclature){
  # Get info for aliases
  lus <- function(x) length(unique(setdiff(x,"")))
  assignees <- split_column(nomenclature,"cell_set_alias_assignee",expand=TRUE)
  citations <- split_column(nomenclature,"cell_set_alias_citation",expand=TRUE)
  aliases   <- split_column(nomenclature,"cell_set_additional_alias")
  numAliases <- apply(aliases,1,lus)
  out <- NULL
  for (i in 1:length(numAliases)) if(numAliases[i]>=1) 
   for (j in 1:numAliases[i]){
    out <- rbind(out,c(preferredAlias[i],aliases[i,j],assignees[i,j+1],citations[i,j+1]))
  }
  colnames(out) <- c("CellSetPreferredAlias","AliasName","AliasAssignedBy","AliasCitation")
  
  # Get info for cell_set_labels and return both
  outLabel <- data.frame(
    CellSetPreferredAlias = preferredAlias,
	AliasName = nomenclature$cell_set_label,
	AliasAssignedBy = "CCN",
	AliasCitation = "github.com/AllenInstitute/nomenclature"
  )
  rbind(as.data.frame(out),outLabel)
}


create_zip_archive = function(prefix = "cell_set_archive"){
  # Make a manifest file
  manifest = '{
"Sample": "Sample.csv",
"CellSet": "CellSet.csv",
"CellSetToCellSet": "CellSetToCellSet.csv",
"CellSetToCellTypeAlias": "CellSetToCellTypeAlias.csv"
}'
  manifestFile = paste0(prefix,"_manifest.json")
  cat(manifest, file = manifestFile)
  
  # Add all relevant files to a zip file with desired filename
  files = c(manifestFile,"CellSetToCellSet.csv","CellSet.csv",
            "CellSetToCellTypeAlias.csv","Sample.csv")
  zip(paste0(prefix,".zip"), files=files)
}


