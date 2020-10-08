#' Generate data files required for shiny app
#'
#' Generate data files required for shiny app. Five files will be generated, 
#' namely (i) the shinycell config \code{prefix_conf.rds}, (ii) the gene 
#' mapping object config \code{prefix_gene.rds}, (iii) the single-cell gene 
#' expression \code{prefix_gexpr.h5}, (iv) the single-cell metadata 
#' \code{prefix_meta.rds} and (v) the defaults for the Shiny app 
#' \code{prefix_def.rds}. A prefix is specified in each file to allow for the 
#' use of multiple single-cell datasets in one Shiny app. Note that both 
#' \code{makeShinyFiles} and \code{makeShinyCodes} functions are ran when 
#' running the wrapper function \code{makeShinyApp}.
#'
#' @param obj input single-cell data object. Both Seurat objects (v3+) and 
#'   SingleCellExperiment objects are accepted.
#' @param scConf shinycell config data.table
#' @param gex.assay assay in single-cell data object to use for plotting 
#'   gene expression. Default is to either use the "RNA" assay for Seurat 
#'   objects (v3+) or the "logcounts" assay for SingleCellExperiment objects
#' @param gex.slot slot in single-cell assay to plot. This is only used 
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param gene.mapping specify a mapping to convert gene identifiers. A 
#'   named vector must be supplied where \code{names(gene.mapping)} correspond 
#'   to the actual gene identifiers in the gene expression matrix whereas 
#'   \code{gene.mapping} correspond to new identifiers to map to. Partial 
#'   mapping is allowed where the mapping is only provided for some gene 
#'   identifiers, the remaining gene identifiers will remain unchanged and a 
#'   warning message will be presented. Default is \code{NA} which is not to 
#'   perform any mapping
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to 
#'   show in bubbleplot / heatmap
#' @param default.dimred character vector specifying two default dim. reduction 
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r
#'
#' @examples
#' makeShinyFiles(seu, scConf, gex.assay = "RNA", gex.slot = "data",
#'                shiny.prefix = "sc1", shiny.dir = "shinyApp/",
#'                default.gene1 = "GATA3", default.gene2 = "DNMT3L",
#'                default.multigene = c("ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
#'                                      "DPPA5","SLC7A2","GATA3","KRT19"),
#'                default.dimred = c("UMAP_1", "UMAP_2")
#'
#' @export
makeShinyFiles <- function(
  obj, scConf, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"), 
  gene.mapping = NA, shiny.prefix, shiny.dir = "shinyApp/",
  default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
  default.dimred = c("UMAP_1", "UMAP_2")){
  ### Preprocessing and checks
  # Generate defaults for gex.assay / gex.slot
  if(class(obj)[1] == "Seurat"){
    if(is.na(gex.assay[1])){gex.assay = "RNA"}
    gex.matrix = slot(obj@assays[[gex.assay[1]]], gex.slot[1])
    sc1meta = data.table(sampleID = rownames(obj@meta.data), obj@meta.data)
  } else if (class(obj)[1] == "SingleCellExperiment"){
    if(is.na(gex.assay[1])){gex.assay = "logcounts"}
    gex.matrix = assay(obj, gex.assay[1])
    sc1meta = data.table(sampleID = rownames(obj@colData), obj@colData)
  } else {
    stop("Only Seurat or SingleCellExperiment objects are accepted!")
  }
  gex.ident = rownames(gex.matrix)
  
  # Perform gene.mapping if provided
  if(!is.na(gene.mapping[1])){
    # Check if gene.mapping is partial or not
    if(!all(gex.ident %in% names(gene.mapping))){
      warning("Mapping for some gene identifiers are not provided!")
      tmp1 = gex.ident[gex.ident %in% names(gene.mapping)]
      tmp1 = gene.mapping[tmp1]
      tmp2 = gex.ident[!gex.ident %in% names(gene.mapping)]
      names(tmp2) = tmp2
      gene.mapping = c(tmp1, tmp2)
    } 
    gene.mapping = gene.mapping[gex.ident]
    names(gene.mapping) = NULL
    rownames(gex.matrix) = gene.mapping
  } else {
    gene.mapping = rownames(gex.matrix)
  }
  
  # Check that default.gene1 / default.gene2 / default.multigene exist
  if(!default.gene1[1] %in% gene.mapping){
    stop("default.gene1 does not exist in gene expression!")
  }
  if(!default.gene2[1] %in% gene.mapping){
    stop("default.gene2 does not exist in gene expression!")
  }
  if(!all(default.multigene %in% gene.mapping)){
    stop("default.multigene does not exist in gene expression!")
  }
  
  
  
  ### Actual object generation
  # Make XXXmeta.rds and XXXconf.rds (updated with dimred info)
  sc1conf = scConf
  sc1conf$dimred = FALSE
  sc1meta = sc1meta[, c("sampleID", as.character(sc1conf$ID)), with = FALSE]
  # Extract dimred and append to both XXXmeta.rds and XXXconf.rds...
  if(class(obj)[1] == "Seurat"){
    for(iDR in names(obj@reductions)){
      # Extract dimred and append to sc1meta
      drMat = obj@reductions[[iDR]]@cell.embeddings
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = drMat[sc1meta$sampleID, ]          # Ensure ordering
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)            
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }
  } else if (class(obj)[1] == "SingleCellExperiment"){
    for(iDR in names(obj@reducedDims)){
      # Extract dimred and append to sc1meta
      drMat = obj@reducedDims[[iDR]]
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = drMat[sc1meta$sampleID, ]          # Ensure ordering
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)            
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }
  }
  sc1conf$ID = as.character(sc1conf$ID)     # Remove levels
  
  # Make XXXgexpr.h5
  if(!dir.exists(shiny.dir)){dir.create(shiny.dir)}
  filename = paste0(shiny.dir, "/", shiny.prefix, "gexpr.h5")
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",  dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = dim(gex.matrix), maxdims = dim(gex.matrix)),
    chunk_dims = c(1,dim(gex.matrix)[2]))
  sc1gexpr.grp.data[, ] <- as.matrix(gex.matrix[, sc1meta$sampleID])
  sc1gexpr$close_all()
  
  # Make XXXgenes.rds
  sc1gene = seq(dim(gex.matrix)[1])
  names(sc1gene) = rownames(gex.matrix)
  
  # Make XXXdef.rds (list of defaults)
  if(all(default.dimred %in% sc1conf[dimred == TRUE]$ID)){
    default.dimred[1] = sc1conf[ID == default.dimred[1]]$UI
    default.dimred[2] = sc1conf[ID == default.dimred[2]]$UI
  } else if(all(default.dimred %in% sc1conf[dimred == TRUE]$UI)) {
    default.dimred = default.dimred    # Nothing happens
  } else {
    # Try to guess... and give a warning
    guess = gsub("[0-9]", "", default.dimred[1])
    if(length(grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case=TRUE)) >= 2){
      default.dimred = sc1conf[dimred == TRUE]$UI[
        grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case = TRUE)[1:2]]
    } else {
      nDR = length(sc1conf[dimred == TRUE]$UI)
      default.dimred = sc1conf[dimred == TRUE]$UI[(nDR-1):nDR]
    }
    warning(paste0("default.dimred not found, switching to ", 
                   default.dimred[1], " and ", default.dimred[1]))
  }
  # Note that we stored the display name here
  sc1def = list()
  sc1def$meta1 = scConf[default == 1]$UI    # Use display name
  sc1def$meta2 = scConf[default == 2]$UI    # Use display name 
  sc1def$gene1 = default.gene1[1]           # Actual == Display name
  sc1def$gene2 = default.gene2[1]           # Actual == Display name
  sc1def$genes = default.multigene          # Actual == Display name
  sc1def$dimred = default.dimred            # Use display name 
  scConf = scConf[, -"default", with = FALSE]
  
  
  
  ### Saving objects
  saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  saveRDS(sc1meta, file = paste0(shiny.dir, "/", shiny.prefix, "meta.rds"))
  saveRDS(sc1gene, file = paste0(shiny.dir, "/", shiny.prefix, "gene.rds"))
  saveRDS(sc1def, file = paste0(shiny.dir, "/", shiny.prefix, "def.rds"))
  return(sc1conf)
}


