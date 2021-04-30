#' doDD
#' Function to find differentially distributed genes
#' 
#' @importFrom stats ks.test p.adjust
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @examples 
#' set.seed(1234)
#' n = 1000 ## genes, rows
#' p = 100 ## cells, cols
#' exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
#' rownames(exprsMat) = paste0("gene", 1:n)
#' colnames(exprsMat) = paste0("cell", 1:p)
#' cellTypes = sample(letters[1:3], size = p, replace = TRUE)
# system.time({doDD(exprsMat = exprsMat, cellTypes = cellTypes)})
#' ####### an example where `sce`` is a SingleCellExperiment object
#' mat <- logcounts(sce)
#' labels <- sce$lineage
#' 
#' dd_res <- doDD(mat, labels, exprs_pct=0.05, filter=FALSE)
doDD <- function(exprsMat, cellTypes){
    cellTypes <- droplevels(as.factor(cellTypes))
    tt <- list()
    for (i in 1:nlevels(cellTypes)) {
        tmp_celltype <- ifelse(cellTypes == levels(cellTypes)[i], 1, 0)
        tt[[i]] <- t(apply(exprsMat, 1, function(x) {
            x1 <- x[tmp_celltype == 0]
            x2 <- x[tmp_celltype == 1]
            ks <- stats::ks.test(x1, x2, alternative = "greater")
            return(c(stats=ks$statistic, 
                     pvalue=ks$p.value))
        }))
        tt[[i]] <- as.data.frame(tt[[i]])
        tt[[i]]$adj.pvalue <- stats::p.adjust(tt[[i]]$pvalue, method = "BH")
    }
    names(tt) <- levels(cellTypes)
    return(tt)
}

doDDPairwise <- function(exprsMat, cellTypes){
    message("DD")
    cty <- droplevels(as.factor(cellTypes))
    names(cty) <- colnames(exprsMat)

    ttAve <- list()
    for (i in 1:nlevels(cty)) {
        
        message(paste0(levels(cty)[[i]], "...."))
        
        tt <- list()
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        for (j in seq_along(other_celltype)) {
            tmp_exprsMat <- exprsMat[,cty %in% c(anchor_celltype, other_celltype[[j]])]
            tmp_celltype <- cty[cty %in% c(anchor_celltype, other_celltype[[j]])]
            tmp_celltype <- (ifelse(tmp_celltype == anchor_celltype, 1, 0))
            
            tt[[j]] <- t(apply(tmp_exprsMat, 1, function(x) {
                x1 <- x[tmp_celltype == 0]
                x2 <- x[tmp_celltype == 1]
                ks <- stats::ks.test(x1, x2, alternative = "greater")
                return(c(stats=ks$statistic, 
                         pvalue=ks$p.value))
            }))
            tt[[j]] <- as.data.frame(tt[[j]])
            tt[[j]]$adj.pvalue <- stats::p.adjust(tt[[j]]$pvalue, method = "BH")
        }
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) { 
            stat = x$`stats.D^+`
            names(stat) = rownames(x)
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}
