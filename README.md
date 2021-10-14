# nomenclature

### Overview

This repository contains code to generate the standardized **Common Cell type Nomenclature (CCN)** from an R "dendrogram", which will be presented in our manuscript soon to be resubmitted to eLife.  The initial version of the CCN is described in detail [on our website, here](https://portal.brain-map.org/explore/classes/nomenclature), with and updated version **now published [on eLife](https://elifesciences.org/articles/59928)**.  We will (hopefully) update our website once the CCN is finalized.  If you want to apply nomenclature of the format initially described, please use the initial release (not recommended).

### Application of CCN to a taxonomy

#### To get started, download the `scripts` and `data` folders, and the follow the directions in `build_annotation_tables.rmd`.  An HTML version of this file is accessible **[AT THIS LINK - build_annotation_tables.nb.html](http://htmlpreview.github.io/?https://github.com/AllenInstitute/nomenclature/blob/master/scripts/build_annotation_tables.nb.html)**.

Please visit [the Cell Taxonomy category of our community forum](https://community.brain-map.org/c/cell-taxonomies) to provide suggestions about the CCN itself.

**We are currently saving used taxonomy_id [IN THIS TABLE](https://docs.google.com/spreadsheets/d/10gYNyOhc0YHOYKjgsvLfumf65CqLeiVDCE3Wrxz4Txo/edit?usp=sharing).  Please choose a unique one and update the table accordingly if you use the CCN.**  A cell type taxonomy service is currently under construction at the Allen Institute, which will will replace this non-ideal system.  CCN output for several of the taxonomies in this table are included in [**Supplementary File 1**](https://github.com/AllenInstitute/nomenclature/raw/master/data/Supplementary_File_1_from_Miller_et_al_2020.zip) from [our paper on eLife](https://elifesciences.org/articles/59928).

### Conversions to other formats

The CCN is intended to be a standardized format for researchers to save taxonomies resulting from any cell typing study.  To facilitate ingest of these taxonomies into relevant databases and ontologies, we intend to provide conversion scripts (as available) from the CCN format to formats required by these external tools.  Currently, the following conversions are available:  
1. **Allen Institute Cell Type Taxonomy Service (CTTS)**: This is a tool for saving taxonomies in a database that is *currently only available to Allen Institute employees*.  To convert to this format download the three "convert_to_CTTS" from [this folder](https://github.com/AllenInstitute/nomenclature/tree/master/scripts/conversion_scripts) and follow the identical directions in the [.Rmd](https://github.com/AllenInstitute/nomenclature/blob/master/scripts/conversion_scripts/convert_to_CTTS.rmd) or [.html](http://htmlpreview.github.io/?https://github.com/AllenInstitute/nomenclature/blob/master/scripts/conversion_scripts/convert_to_CTTS.nb.html) files.
2. **Cell Annotation Platform (CAP)**: This is a tool being built as part of the Human Cell Atlas (HCA) for annotating taxonomies.  Scripts for injest into this service are in process.  Current content of CAP is available [here](http://celltype.info/).  

If you would like to see scripts for conversions to other file formats, please email Jeremy Miller (jeremym _at_ alleninstitute _dot_ org).

### License

The license for this package is available on Github at: https://github.com/AllenInstitute/allen_institute_nomenclature/blob/master/LICENSE

### Level of Support

This tool will only be updated to correct errors and to reflect updates in the nomenclature schema.

### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/allen_institute_nomenclature/blob/master/CONTRIBUTION
