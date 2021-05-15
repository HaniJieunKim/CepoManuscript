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
        
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        tt = furrr::future_map(
            .x = other_celltype,
            .f = ~ compute_LRPairwise_once(other = .x,
                                           seuObj = seuObj,
                                           anchor = anchor_celltype), 
            .progress = TRUE)
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

compute_LRPairwise_once = function(seuObj, anchor, other) {
    
    res = Seurat::FindMarkers(seuObj, 
                              ident.1 = anchor, 
                              ident.2 = other, 
                              test.use = "LR",  
                              min.cells.feature = 0,
                              min.diff.pct = -Inf,
                              max.cells.per.ident = Inf,
                              min.cells.group = 0,
                              logfc.threshold = 0,
                              min.pct = 0)
    stats = res$p_val
    logFC = grep("FC", colnames(res), value = TRUE)
    stats = sapply(1:length(res[,logFC]), function(x) {
        if (res[,logFC][[x]] > 0) { stats[[x]]} else { 1 - stats[[x]] }
    })
    names(stats) = rownames(res)
    stats = -log10(stats)
    stats[is.infinite(stats)] = 0
    
    return(stats)
    
}


doROCPairwise <- function(seuObj, cellTypes){
    
    cty <- droplevels(as.factor(cellTypes))
    names(cty) <- colnames(seuObj)
    
    ttAve <- list()
    for (i in 1:nlevels(cty)) {
        
        message(paste0(levels(cty)[[i]], "...."))
        
        tt <- list()
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        for (j in seq_along(other_celltype)) {

            res = Seurat::FindMarkers(seuObj, 
                                      ident.1 = anchor_celltype, 
                                      ident.2 = other_celltype, 
                                      test.use = "roc", 
                                      min.pct = 0, 
                                      min.cells.feature = 0,
                                      min.cells.group = 0,
                                      logfc.threshold = 0)
            tt[[j]] = res
        }
        
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) { 
            stat = x$power
            
            logFC = grep("FC", colnames(x), value = TRUE)
            stat = sapply(1:length(x[,logFC]), function(y) {
                if (x[,logFC][[y]] > 0) { stat[[y]]} else { -(1 - stat[[y]]) }
            })
            names(stat) = rownames(x)
            stat[is.na(stat)] <- 0
            stat[is.infinite(stat)] <- 0
            
            stat <- sort(stat, decreasing=TRUE, na.last=TRUE)
            
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
    
}
