suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(data.table))

library(limma)
#' @importFrom limma eBayes lmFit
#' @importFrom methods new
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @examples 
#' 
#' ####### an example where `sce`` is a SingleCellExperiment object
#' mat <- logcounts(sce)
#' labels <- sce$lineage
#' 
#' ds_res <- Cepo(mat, labels, exprs_pct=0.05, filter=FALSE)
#' 
#' ####### A quick simulation. Note that we need to have colnames and rownames on matrix ##########
#' set.seed(1234)
#' n = 1000 ## genes, rows
#' p = 100 ## cells, cols
#' exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
#' rownames(exprsMat) = paste0("gene", 1:n)
#' colnames(exprsMat) = paste0("cell", 1:p)
#' cellTypes = sample(letters[1:3], size = p, replace = TRUE)
#' doLimma(exprsMat = exprsMat, cellTypes = cellTypes)
#' doChiSquared(exprsMat = exprsMat, cellTypes = cellTypes)
#' doBI(exprsMat = exprsMat, cellTypes = cellTypes)
#' doChiSquared2(exprsMat = exprsMat, cellTypes = cellTypes)
#' microbenchmark::microbenchmark(
#' doLimma(exprsMat = exprsMat, cellTypes = cellTypes),
#' doChiSquared(exprsMat = exprsMat, cellTypes = cellTypes),
#' doBI(exprsMat = exprsMat, cellTypes = cellTypes),
#' doChiSquared2(exprsMat = exprsMat, cellTypes = cellTypes),
#' times = 10
#' )
#' identical(doChiSquared2(exprsMat = exprsMat, cellTypes = cellTypes), doChiSquared(exprsMat = exprsMat, cellTypes = cellTypes))
#' bi = doBI_KW(exprsMat, cellTypes)[[1]]
#' limma = doLimma(exprsMat = exprsMat, cellTypes = cellTypes)[[1]]
#' plot(bi, limma[names(bi),]$t)
doLimmaPairwise <- function(exprsMat, cellTypes) {
    # input must be normalised, log-transformed data
    message("Limma")
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
            
            design <- stats::model.matrix(~tmp_celltype)
            
            y <- methods::new("EList")
            y$E <- tmp_exprsMat
            fit <- limma::lmFit(y, design = design)
            fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
            tt[[j]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
            if (!is.null(tt[[j]]$ID)) {
                tt[[j]] <- tt[[j]][!duplicated(tt[[j]]$ID),]
                rownames(tt[[j]]) <- tt[[j]]$ID
            }
            
        }
        
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) { 
            stat = x$t
            names(stat) = rownames(x)
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}

doEdgeRPairwise <- function(countMat, cellTypes) {
    
    message("edgeRQLFDetRate")
    cty <- droplevels(as.factor(cellTypes))
    names(cty) <- colnames(countMat)
    
    ttAve <- list()
    for (i in 1:nlevels(cty)) {
        
        tt <- df <- list()
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        for (j in seq_along(other_celltype)) {
            
            tmp_countMat <- countMat[,cty %in% c(anchor_celltype, other_celltype[[j]])]
            tmp_celltype <- cty[cty %in% c(anchor_celltype, other_celltype[[j]])]
            grp <- (ifelse(tmp_celltype == anchor_celltype, 1, 0))
            names(grp) <- names(tmp_celltype)
            
            dge <- DGEList(tmp_countMat, group = tmp_celltype)
            dge <- calcNormFactors(dge)
            cdr <- scale(colMeans(tmp_countMat > 0))
            design <- model.matrix(~ cdr + grp)
            dge <- estimateDisp(dge, design = design)
            fit <- glmQLFit(dge, design = design)
            qlf <- glmQLFTest(fit)
            tt <- topTags(qlf, n = Inf)
            
            #plotBCV(dge)
            #plotQLDisp(fit)
            #hist(tt$table$PValue, 50)
            #hist(tt$table$FDR, 50)
            #limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
            #plotSmear(qlf)
            
            df[[j]] = tt$table
            
        }
        names(df) <- other_celltype
        
        idx = rownames(df[[1]])
        stats = rowMeans(do.call(cbind, lapply(df, function(x) {
            
            res = x$F
            names(res) = rownames(x)
            res <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {res[y]} else {-res[y]}
            })
            return(res[idx])
            
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
    
}


###### Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
doMASTPairwise <- function(exprsMat, cellTypes) {
    # input must be normalised, log-transformed data
    message("MAST")
    cty <- droplevels(as.factor(cellTypes))
    names(cty) <- colnames(exprsMat)
    
    ttAve <- list()
    for (i in 1:nlevels(cty)) {
        message(paste0(levels(cty)[[i]], "...."))
        
        df <- list()
        anchor_celltype = levels(cty)[[i]]
        other_celltype = setdiff(levels(cty), anchor_celltype)
        
        for (j in seq_along(other_celltype)) {
            
            tmp_exprsMat <- exprsMat[,cty %in% c(anchor_celltype, other_celltype[[j]])]
            tmp_celltype <- cty[cty %in% c(anchor_celltype, other_celltype[[j]])]
            grp <- (ifelse(tmp_celltype == anchor_celltype, 1, 0))
            
            names(grp) <- names(tmp_celltype)
            cdr <- scale(colMeans(tmp_exprsMat > 0))
            sca <- FromMatrix(exprsArray = tmp_exprsMat, 
                              cData = data.frame(wellKey = names(grp), 
                                                 grp = grp, cdr = cdr))
            zlmdata <- zlm(~cdr + grp, sca)
            mast <- lrTest(zlmdata, "grp")
            
            df[[j]] = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                                 lambda = mast[, "cont", "lambda"],
                                 row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
            df[[j]]$fdr <- stats::p.adjust(df[[j]]$pval, method="BH")
            b <- getLogFC(zlmdata)
            df[[j]] <- cbind(df[[j]], b[b$contrast=="grp",c("logFC", "varLogFC","z"), with=F])
            df[[j]] <- df[[j]][order(df[[j]]$pval),]
        }
        
        idx = rownames(df[[1]])
        stats = rowMeans(do.call(cbind, lapply(df, function(x) {
            
            mast <- x$pval
            names(mast) <- rownames(x)
            x$logFC[is.na(x$logFC)] <- 0
            stat <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {mast[y]} else {1-mast[y]}
            })
            stat <- -log10(stat)
            stat[is.na(stat)] <- 0
            stat[is.infinite(stat)] <- 0
            
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
    
}

###### Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(genefilter))
doTestPairwise <- function(exprsMat, cellTypes) {
    # input must be normalised, log-transformed data
    message("t-test")
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
                
                res <- stats::t.test(x2, y=x1)
                return(c(stats=res$statistic,
                         pvalue=res$p.value))
            }))
            tt[[j]] <- as.data.frame(tt[[j]])
            tt[[j]]$adj.pvalue <- stats::p.adjust(tt[[j]]$pvalue, method = "BH")
        }
        
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) {
            stat = x$stats.t
            names(stat) = rownames(x)
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}


###### Adapted from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R
suppressPackageStartupMessages(library(edgeR))
doWilcoxonPairwise <- function(exprsMat, cellTypes) {
    # input must be normalised, log-transformed data
    message("Wilcoxon")
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
                
                res <- stats::wilcox.test(x ~ tmp_celltype)
                c(stats=res$statistic,
                  pvalue=res$p.value)
            }))
            tt[[j]] <- as.data.frame(tt[[j]])
            tt[[j]]$adj.pvalue <- stats::p.adjust(tt[[j]]$pvalue, method = "BH")
        }
        
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) {
            stat = -x$stats.W
            names(stat) = rownames(x)
            return(stat[idx])
        })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        
    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}





###### Adapted from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_voomlimma.R
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
doVoomPairwise <- function(exprsMat, cellTypes) {
    # input must be normalised, log-transformed data
    message("Voom Limma")
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
            
            design <- stats::model.matrix(~tmp_celltype)
            
            y <- methods::new("EList")
            #y$E <- exprsMat[keep, ]
            y$E <- tmp_exprsMat
            vm <- limma::voom(y, design = design)
            fit <- limma::lmFit(vm, design = design)
            fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
            tt[[j]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
            if (!is.null(tt[[j]]$ID)) {
                tt[[j]] <- tt[[j]][!duplicated(tt[[j]]$ID),]
                rownames(tt[[j]]) <- tt[[j]]$ID
            }
            
        }
        
        idx = rownames(tt[[1]])
        stats = rowMeans(do.call(cbind, lapply(tt, function(x) { 
                stat = x$t
                names(stat) = rownames(x)
                return(stat[idx])
            })))
        ttAve[[i]] <- sort(stats, decreasing=T, na.last=T)
        

    }
    names(ttAve) <- levels(cty)
    return(ttAve)
}