#' doSCVI
#' Function to find differentially expressed genes using scVI
#' 
#' @import reticulate
#' @param seuObj Seurat obj where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' 
#' @examples 
#' library(Cepo)
#' data(cellbench)
#' seu = Seurat::as.Seurat(cellbench)
#' seu@active.ident = as.factor(seu$celltype)
#' cty = as.factor(seu$celltype)
#' res <- doSCVI(seu, labels)

doSCVI <- function(seuObj, cellTypes, n_epochs=400){
    
    require(reticulate) 
    conda_list()
    use_condaenv(condaenv = 'scvi-env', required = TRUE)
    
    sc <- import('scanpy', convert = FALSE)
    scvi <- import('scvi', convert = FALSE)
    scvi$settings$progress_bar_style = 'tqdm'
    
    seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = nrow(seuObj))
    
    adata <- sc$AnnData(
        X   = t(as.matrix(GetAssayData(seuObj,slot='counts'))), #scVI requires raw counts
        obs = seuObj[[]],
        var = GetAssay(seuObj)[[]]
    )
    
    # run seteup_anndata
    scvi$data$setup_anndata(adata)
    
    # create the model
    #model = scvi$model$SCVI(adata, use_cuda = TRUE)
    model = scvi$model$SCVI(adata)
    
    # train the model
    model$train(max_epochs = as.integer(n_epochs))
    
    # get the latent represenation
   #latent = model$get_latent_representation()
   #
   ## put it back in our original Seurat object
   #latent <- as.matrix(latent)
   #rownames(latent) = colnames(seuObj)
   #seuObj[['scvi']] <- CreateDimReducObject(embeddings = latent, 
   #                                         key = "scvi_", assay = DefaultAssay(seuObj))
   #
   #seuObj <- FindNeighbors(seuObj, dims = 1:10, reduction = 'scvi')
   #seuObj <- FindClusters(seuObj, resolution =1)
   #
   #seuObj <- RunUMAP(seuObj, dims = 1:10, reduction = 'scvi', n.components = 2)
    
    #DimPlot(seuObj, reduction = "umap", pt.size = 3)
    
    # DE 
    # we need pandas in order to put the seurat clusters into our original anndata
    pd <- import('pandas', convert = FALSE)
    cp <- import('copy', convert = FALSE)
    
    cty <- droplevels(as.factor(cellTypes))
    
    # workaround for adding seurat_clusters to our anndata
    #new_obs <- c(adata$obs, pd$DataFrame(seuObj[['celltype']]))

    tt <- list()
    for (i in 1:nlevels(cty)) {
        
        print(levels(cty)[[i]])
        adata_tmp = cp$copy(adata)
        tmp_celltype = as.integer(ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
        tmp_celltype = ifelse(cellTypes == levels(cellTypes)[i], 1, 0)
        tmp_celltype = data.frame(cellType = tmp_celltype)
        rownames(tmp_celltype) = colnames(seuObj)
        
        new_obs <- c(adata_tmp$obs, pd$DataFrame(tmp_celltype))
        new_obs <- pd$concat(new_obs, axis = 1)
        
        adata_tmp$obs <- new_obs
        
        DE <- model$differential_expression(adata_tmp, 
                                            groupby='cellType', 
                                            group1 = 1, group2 = 0)
        
        tt[[i]] = reticulate::py_to_r(DE)
        
    }
    names(tt) <- levels(cty)
    return(tt)
}
