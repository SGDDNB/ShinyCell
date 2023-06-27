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
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment 
#'   data or input file path for h5ad / loom files
#' @param scConf shinycell config data.table
#' @param gex.assay assay in single-cell data object to use for plotting 
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay, 
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts", 
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'     \item{loom files}: "matrix" or any assay in "layers",
#'       default is "matrix"
#'   }
#' @param gex.slot slot in single-cell assay to plot. This is only used 
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param gene.mapping specifies whether to convert human / mouse Ensembl gene 
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols. 
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is 
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users 
#'   can supply a named vector where \code{names(gene.mapping)} correspond 
#'   to the actual gene identifiers in the gene expression matrix and 
#'   \code{gene.mapping} correspond to new identifiers to map to
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to 
#'   show in bubbleplot / heatmap
#' @param default.dimred character vector specifying the two default dimension 
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param chunkSize number of genes written to h5file at any one time. Lower 
#'   this number to reduce memory consumption. Should not be less than 10
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate hdf5r
#'
#' @examples
#' makeShinyFiles(seu, scConf, gex.assay = "RNA", gex.slot = "data",
#'                shiny.prefix = "sc1", shiny.dir = "shinyApp/",
#'                default.gene1 = "GATA3", default.gene2 = "DNMT3L",
#'                default.multigene = c("ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
#'                                      "DPPA5","SLC7A2","GATA3","KRT19"),
#'                default.dimred = c("UMAP_1", "UMAP_2"))
#'
#' @export
makeShinyFiles <- function(
  obj, scConf, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"), 
  gene.mapping = FALSE, shiny.prefix = "sc1", shiny.dir = "shinyApp/",
  default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
  default.dimred = NA, chunkSize = 500, markers.all=FALSE, markers.top20=FALSE, de.genes=FALSE, gene.ranks=FALSE){
  ### Preprocessing and checks
  # Generate defaults for gex.assay / gex.slot
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    if(is.na(gex.assay[1])){gex.assay = "RNA"}
    gex.matdim = dim(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
    gex.rownm = rownames(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
    gex.colnm = colnames(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
    # defGenes = obj@assays[[gex.assay[1]]]@var.features[1:10]
    defGenes = Seurat::VariableFeatures(obj)[1:10]
    if(is.na(defGenes[1])){
      warning(paste0("Variable genes for seurat object not found! Have you ",
                     "ran `FindVariableFeatures` or `SCTransform`?"))
      defGenes = gex.rownm[1:10]
    }
    sc1meta = data.table(sampleID = rownames(obj@meta.data), obj@meta.data)
    
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    if(is.null(colnames(obj)[1])){
      colnames(obj) = paste0("cell_", seq(ncol(obj)))
    }    # Populate cell IDs if they are not present
    if(is.na(gex.assay[1])){gex.assay = "logcounts"}
    gex.matdim = dim(SummarizedExperiment::assay(obj, gex.assay[1]))
    gex.rownm = rownames(SummarizedExperiment::assay(obj, gex.assay[1]))
    gex.colnm = colnames(SummarizedExperiment::assay(obj, gex.assay[1]))
    defGenes = gex.rownm[1:10]
    sc1meta = SingleCellExperiment::colData(obj)
    sc1meta = data.table(sampleID = rownames(sc1meta), 
                         as.data.frame(sc1meta))
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    if(is.na(gex.assay[1])){gex.assay = "X"}
    # Can just check X since inpH5$layers should have same dimensions
    ad <- import("anndata", convert = FALSE)
    sp <- import('scipy.sparse', convert = FALSE)
    inpH5 = ad$read_h5ad(obj)
    gex.matdim = rev(unlist(py_to_r(inpH5$X$shape)))  
    gex.rownm = py_to_r(inpH5$var_names$values)
    gex.colnm = py_to_r(inpH5$obs_names$values)
    defGenes = gex.rownm[1:10]
    sc1meta = data.table(sampleID = gex.colnm)
    sc1meta = cbind(sc1meta, data.table(py_to_r(inpH5$obs$values)))
    colnames(sc1meta) = c("sampleID", py_to_r(inpH5$obs$columns$values))
    for(i in colnames(sc1meta)[-1]){
      sc1meta[[i]] = unlist(sc1meta[[i]])   # unlist and refactor
      if(as.character(inpH5$obs[i]$dtype) == "category"){
        sc1meta[[i]] = factor(sc1meta[[i]], levels = 
                                py_to_r(inpH5$obs[i]$cat$categories$values))
      }
    } 

  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    if(is.na(gex.assay[1])){gex.assay = "matrix"}
    # Can just check matrix since inpLM[["layers"]] should have same dimensions
    inpLM = H5File$new(obj, mode = "r+")
    gex.matdim = rev(inpLM[["matrix"]]$dims)
    gex.rownm = inpLM[["row_attrs"]][["Gene"]]$read()
    for(i in unique(gex.rownm[duplicated(gex.rownm)])){
      gex.rownm[gex.rownm == i] = paste0(i, "-", seq(sum(gex.rownm == i)))
    } # make unique gene names
    gex.colnm = inpLM[["col_attrs"]][["CellID"]]$read()
    defGenes = gex.rownm[1:10]
    cellIdx = which(inpLM[["col_attrs"]]$names == "CellID")
    sc1meta = data.table(sampleID = gex.colnm)
    for(i in inpLM[["col_attrs"]]$names[-cellIdx]){
      tmp = inpLM[["col_attrs"]][[i]]$read()
      if(length(tmp) == nrow(sc1meta)){sc1meta[[i]] = tmp}
    }
     
  } else {
    stop("Only Seurat/SCE objects or h5ad/loom file paths are accepted!")
  }
  
  # Perform gene.mapping if specified (also map defGenes)
  if(gene.mapping[1] == TRUE){
    if(sum(grepl("^ENSG000", gex.rownm)) >= sum(grepl("^ENMUSG000", gex.rownm))){
      tmp1 = fread(system.file("extdata", "geneMapHS.txt.gz", 
                               package = "ShinyCell"))
    } else {
      tmp1 = fread(system.file("extdata", "geneMapMM.txt.gz", 
                               package = "ShinyCell"))
    }
    gene.mapping = tmp1$geneName
    names(gene.mapping) = tmp1$geneID
  }
  # Check if gene.mapping is partial or not
  if(gene.mapping[1] == FALSE){
    gene.mapping = gex.rownm      
    names(gene.mapping) = gex.rownm    # Basically no mapping
  } else {
    if(!all(gex.rownm %in% names(gene.mapping))){
      # warning("Mapping for some gene identifiers are not provided!")
      tmp1 = gex.rownm[gex.rownm %in% names(gene.mapping)]
      tmp1 = gene.mapping[tmp1]
      tmp2 = gex.rownm[!gex.rownm %in% names(gene.mapping)]
      names(tmp2) = tmp2
      gene.mapping = c(tmp1, tmp2)
    } 
    gene.mapping = gene.mapping[gex.rownm]
  }
  defGenes = gene.mapping[defGenes]
  
  # Check default.gene1 / default.gene2 / default.multigene
  default.gene1 = default.gene1[1]
  default.gene2 = default.gene2[1]
  if(is.na(default.gene1)){default.gene1 = defGenes[1]}
  if(is.na(default.gene2)){default.gene2 = defGenes[2]}
  if(is.na(default.multigene[1])){default.multigene = defGenes}
  if(default.gene1 %in% gene.mapping){
    default.gene1 = default.gene1
  } else if(default.gene1 %in% names(gene.mapping)){
    default.gene1 = gene.mapping[default.gene1]
  } else {
    warning("default.gene1 doesn't exist in gene expression, using defaults...")
    default.gene1 = defGenes[1]
  }
  if(default.gene2 %in% gene.mapping){
    default.gene2 = default.gene2
  } else if(default.gene2 %in% names(gene.mapping)){
    default.gene2 = gene.mapping[default.gene2]
  } else {
    warning("default.gene2 doesn't exist in gene expression, using defaults...")
    default.gene2 = defGenes[2]
  }
  if(all(default.multigene %in% gene.mapping)){
    default.multigene = default.multigene
  } else if(all(default.multigene %in% names(gene.mapping))){
    default.multigene = gene.mapping[default.multigene]
  } else {
    warning(paste0("default.multigene doesn't exist in gene expression, ", 
                   "using defaults..."))
    default.multigene = defGenes
  }
  

  ### Actual object generation
  # Make XXXmeta.rds and XXXconf.rds (updated with dimred info)
  sc1conf = scConf
  sc1conf$dimred = FALSE
  sc1meta = sc1meta[, c("sampleID", as.character(sc1conf$ID)), with = FALSE]
  # Factor metadata again
  for(i in as.character(sc1conf[!is.na(fID)]$ID)){
    sc1meta[[i]] = factor(sc1meta[[i]],
                          levels = strsplit(sc1conf[ID == i]$fID, "\\|")[[1]])
    levels(sc1meta[[i]]) = strsplit(sc1conf[ID == i]$fUI, "\\|")[[1]]
    sc1conf[ID == i]$fID = sc1conf[ID == i]$fUI
  }
  # Extract dimred and append to both XXXmeta.rds and XXXconf.rds...
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    for(iDR in names(obj@reductions)){
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
    # SCE Object
    for(iDR in SingleCellExperiment::reducedDimNames(obj)){
      drMat = SingleCellExperiment::reducedDim(obj, iDR)
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      if(is.null(colnames(drMat))){
        colnames(drMat) = paste0(iDR, seq(ncol(drMat)))
      }
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
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    for(iDR in py_to_r(inpH5$obsm_keys())){
      drMat = py_to_r(inpH5$obsm[iDR])
      tmpName = gsub("pca", "pc", gsub("X_", "", iDR))
      tmpName = paste0(tmpName, "_", 1:ncol(drMat))
      colnames(drMat) = tmpName
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }
    
  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    nDR = inpLM[["col_attrs"]]$names[
      grep("pca|tsne|umap", inpLM[["col_attrs"]]$names, ignore.case = TRUE)]
    for(iDR in nDR){
      drMat = t(inpLM[["col_attrs"]][[iDR]]$read())
      colnames(drMat) = paste0(iDR, "_", 1:ncol(drMat))
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
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
    space = H5S$new("simple", dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1,gex.matdim[2]))
  chk = chunkSize
  while(chk > (gex.matdim[1]-8)){
    chk = floor(chk / 2)     # Account for cases where nGene < chunkSize
  } 
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    for(i in 1:floor((gex.matdim[1]-8)/chk)){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        slot(obj@assays[[gex.assay[1]]], gex.slot[1])[((i-1)*chk+1):(i*chk),])
    }
    sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
      slot(obj@assays[[gex.assay[1]]], gex.slot[1])[(i*chk+1):gex.matdim[1],])
    
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    for(i in 1:floor((gex.matdim[1]-8)/chk)){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        SummarizedExperiment::assay(obj, gex.assay[1])[((i-1)*chk+1):(i*chk),])
    }
    sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
      SummarizedExperiment::assay(obj, gex.assay[1])[(i*chk+1):gex.matdim[1],])
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    if(gex.assay == "X"){
      scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$X)))
    } else {
      scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$layers[[gex.assay]])))
    }
    for(i in 1:floor((gex.matdim[1]-8)/chk)){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        scGEX[((i-1)*chk+1):(i*chk),])
    }
    sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
      scGEX[(i*chk+1):gex.matdim[1],])

  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    if(gex.assay == "matrix"){
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- t(
          inpLM[["matrix"]][, ((i-1)*chk+1):(i*chk)])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- t(
        inpLM[["matrix"]][, (i*chk+1):gex.matdim[1]])
    } else {
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- t(
          inpLM[["layers"]][[gex.assay]][, ((i-1)*chk+1):(i*chk)])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- t(
        inpLM[["layers"]][[gex.assay]][, (i*chk+1):gex.matdim[1]])
    }
  }

  # sc1gexpr.grp.data[, ] <- as.matrix(gex.matrix[,])
  sc1gexpr$close_all()
  if(!isTRUE(all.equal(sc1meta$sampleID, gex.colnm))){
    sc1meta$sampleID = factor(sc1meta$sampleID, levels = gex.colnm)
    sc1meta = sc1meta[order(sampleID)]
    sc1meta$sampleID = as.character(sc1meta$sampleID)
  }
  
  # Make XXXgenes.rds
  sc1gene = seq(gex.matdim[1])
  names(gene.mapping) = NULL
  names(sc1gene) = gene.mapping
  sc1gene = sc1gene[order(names(sc1gene))]
  sc1gene = sc1gene[order(nchar(names(sc1gene)))]
  
  # Make XXXdef.rds (list of defaults)
  if(all(default.dimred %in% sc1conf[dimred == TRUE]$ID)){
    default.dimred[1] = sc1conf[ID == default.dimred[1]]$UI
    default.dimred[2] = sc1conf[ID == default.dimred[2]]$UI
  } else if(all(default.dimred %in% sc1conf[dimred == TRUE]$UI)) {
    default.dimred = default.dimred    # Nothing happens
  } else {
    warn = TRUE
    if(is.na(default.dimred[1])){
      default.dimred = "umap"
      warn = FALSE
    }
    # Try to guess... and give a warning
    guess = gsub("[0-9]", "", default.dimred[1])
    if(length(grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case=TRUE)) >= 2){
      default.dimred = sc1conf[dimred == TRUE]$UI[
        grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case = TRUE)[1:2]]
    } else {
      nDR = length(sc1conf[dimred == TRUE]$UI)
      default.dimred = sc1conf[dimred == TRUE]$UI[(nDR-1):nDR]
    }
    if(warn){
      warning(paste0("default.dimred not found, switching to ", 
                     default.dimred[1], " and ", default.dimred[1]))
    } # Warn if user-supplied default.dimred is not found
  }
  # Note that we stored the display name here
  sc1def = list()
  sc1def$meta1 = sc1conf[default == 1]$UI   # Use display name
  sc1def$meta2 = sc1conf[default == 2]$UI   # Use display name 
  sc1def$gene1 = default.gene1              # Actual == Display name
  sc1def$gene2 = default.gene2              # Actual == Display name
  sc1def$genes = default.multigene          # Actual == Display name
  sc1def$dimred = default.dimred            # Use display name 
  tmp = nrow(sc1conf[default != 0 & grp == TRUE])
  if(tmp == 2){
    sc1def$grp1 = sc1def$meta1
    sc1def$grp2 = sc1def$meta2
  } else if(tmp == 1){
    sc1def$grp1 = sc1conf[default != 0 & grp == TRUE]$UI
    if(nrow(sc1conf[default == 0 & grp == TRUE]) == 0){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    }
  } else {
    sc1def$grp1 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    if(nrow(sc1conf[default == 0 & grp == TRUE]) < 2){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[2]
    }
  }
  sc1conf = sc1conf[, -c("fUI", "default"), with = FALSE]
  
  
  
  ### Saving objects
  #saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  saveRDS(sc1meta, file = paste0(shiny.dir, "/", shiny.prefix, "meta.rds"))
  saveRDS(sc1gene, file = paste0(shiny.dir, "/", shiny.prefix, "gene.rds"))
  saveRDS(sc1def,  file = paste0(shiny.dir, "/", shiny.prefix, "def.rds"))

  sc1conf$extra_tabs <- FALSE

  
  if(class(obj)[1]=='Seurat') {
    if(markers.all==TRUE) { 
      if(!is.null(obj@misc$markers$presto$all)){
        m_all <- obj@misc$markers$presto$all
        names(m_all)[names(m_all) == 'group'] <- 'cluster'
        sc1conf$extra_tabs[1] = TRUE
        saveRDS(m_all, file=paste0(shiny.dir, "/", shiny.prefix, "m_all.rds"))
      }
      else {
        print("Warning: 'markers.all' not found in Seurat (structure expected: seurat@misc$markers$presto$all);\n \
          corresponding Shiny tab will be created but with an error message instead of what is expected...\n")
        #sc1conf$extra_tabs[1] = FALSE
      }
    }

    if(markers.top20==TRUE) {
      if(!is.null(obj@misc$markers$presto$top_20)) {
        m_t20 <- obj@misc$markers$presto$top_20
        names(m_t20)[names(m_t20) == 'group'] <- 'cluster'
        sc1conf$extra_tabs[2] = TRUE
        saveRDS(m_t20, file=paste0(shiny.dir, "/", shiny.prefix, "m_t20.rds"))
      }
      else {
        print("Warning: 'markers.top20' not found in Seurat (structure expected: seurat@misc$markers$presto$top_20);\n \
          corresponding Shiny tab will be created but with an error message instead of what is expected...\n")
        #sc1conf$extra_tabs[2] = FALSE
      }
    }

    if(de.genes==TRUE) { 
      if(!is.null(obj@misc$DE_genes$libra$overall)) {
        de_genes <- obj@misc$DE_genes$libra$overall
        sc1conf$extra_tabs[3] = TRUE
        saveRDS(de_genes, file=paste0(shiny.dir, "/", shiny.prefix, "de_genes.rds"))
      }
      else {
        print("Warning: 'de.genes' not found in Seurat (structure expected: seurat@misc$DE_genes$libra$overall);\n \
          corresponding Shiny tab will be created but with an error message instead of what is expected...\n")
        #sc1conf$extra_tabs[3] = FALSE
      }
    }

    if(gene.ranks==TRUE) {
      if(!is.null(obj@misc$gene_ranks$aucell$all)) {
        gene_ranks <- obj@misc$gene_ranks$aucell$all
        sc1conf$extra_tabs[4] = TRUE
        saveRDS(gene_ranks, file=paste0(shiny.dir, "/", shiny.prefix, "gene_ranks.rds"))
      }
      else {
        print("Warning: 'gene.ranks' not found in Seurat (structure expected: seurat@misc$gene_ranks$aucell$all);\n \
          corresponding Shiny tab will be created but with an error message instead of what is expected...\n")
        #sc1conf$extra_tabs[4] = FALSE
      }
    }
  }

  

  #print(sc1conf)
  saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  return(sc1conf)
}


