#' doLR
#' Function to find differentially expressed genes using logistic regression
#' 
#' @importFrom Seurat FindMarkers
#' @param seuObj Seurat obj where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @examples 
#' library(Cepo)
#' data(cellbench)
#' seu = Seurat::as.Seurat(cellbench)
#' seu@active.ident = as.factor(seu$celltype)
#' cty = as.factor(seu$celltype)
#' res <- doLR(seu, labels)

#' doLR
#' Function to find differentially expressed genes using logistic regression
#' 
#' @importFrom Seurat FindMarkers
#' @param seuObj Seurat obj where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @examples 
#' library(Cepo)
#' data(cellbench)
#' seu = Seurat::as.Seurat(cellbench)
#' seu@active.ident = as.factor(seu$celltype)
#' cty = as.factor(seu$celltype)
#' res <- doLR(seu, labels)

doLR <- function(seuObj, cellTypes){
    cty <- droplevels(as.factor(cellTypes))
    
    tt <- list()
    for (i in 1:nlevels(cty)) {
        res = Seurat::FindMarkers(seuObj, 
                                  ident.1 = levels(cty)[[i]], 
                                  test.use = "LR", 
                                  min.pct = 0, 
                                  min.cells.feature = 0,
                                  min.cells.group = 0,
                                  logfc.threshold = 0)
        ######## get statistics outside of function
        #stats = res$p_val_adj
        #stats = sapply(1:length(res$avg_logFC), function(y) {
        #    if (res$avg_logFC[[y]] > 0) { stats[[y]]} else { 1 - stats[[y]] }
        #})
        #names(stats) = rownames(res)
        #stats = -log10(stats)
        #stats[is.infinite(stats)] = 0
        #stats = sort(stats, decreasing = TRUE)
        #tt[[i]] = stats
        tt[[i]] = res
        
    }
    names(tt) <- levels(cty)
    return(tt)
}


doROC <- function(seuObj, cellTypes){
    cty <- droplevels(as.factor(cellTypes))
    
    tt <- list()
    for (i in 1:nlevels(cty)) {
        res = Seurat::FindMarkers(seuObj, 
                                  ident.1 = levels(cty)[[i]], 
                                  test.use = "roc", 
                                  min.pct = 0, 
                                  min.cells.feature = 0,
                                  min.cells.group = 0,
                                  logfc.threshold = 0)
        tt[[i]] = res
        
    }
    names(tt) <- levels(cty)
    return(tt)
}


doBimod <- function(seuObj, cellTypes){
    cty <- droplevels(as.factor(cellTypes))
    
    tt <- list()
    for (i in 1:nlevels(cty)) {
        res = Seurat::FindMarkers(seuObj, 
                                  ident.1 = levels(cty)[[i]], 
                                  test.use = "bimod", 
                                  min.pct = 0, 
                                  min.cells.feature = 0,
                                  min.cells.group = 0,
                                  logfc.threshold = 0)
        tt[[i]] = res
        
    }
    names(tt) <- levels(cty)
    return(tt)
}


doLRPairwise <- function(seuObj, cellTypes){
    
    message("Logistic Regression")
    cty <- droplevels(as.factor(cellTypes))
    
    ttAve <- list()
    for (i in 1:nlevels(cty)) {
        
        message(paste0(levels(cty)[[i]], "...."))
        
        tt <- list()
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        for (j in seq_along(other_celltype)) {
            
            res = Seurat::FindMarkers(seuObj, ident.1 = anchor_celltype, 
                                      ident.2 = other_celltype[[j]], 
                                      test.use = "LR",  
                                      min.cells.feature = 0,
                                      min.cells.group = 0,
                                      logfc.threshold = 0,
                                      min.pct = 0)
            stats = res$p_val_adj
            stats = sapply(1:length(res$avg_logFC), function(x) {
                if (res$avg_logFC[[x]] > 0) { stats[[x]]} else { 1 - stats[[x]] }
            })
            names(stats) = rownames(res)
            stats = -log10(stats)
            stats[is.infinite(stats)] = 0
            tt[[j]] = stats
            
        }
        idx = names(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) { 
            stat = x
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}


