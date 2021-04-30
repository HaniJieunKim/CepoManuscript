plotHeatmap <- function(sce, gene.list, celltype, slim=FALSE, verbose=F, scale=T) {

    mat.ordered <- do.call(cbind, sapply(split(seq_len(ncol(sce)), celltype), function(i) logcounts(sce)[, i]))
    label.ordered <- do.call(c, sapply(split(seq_len(ncol(sce)), celltype), function(i) celltype[i]))
    
    exp.sel <- mat.ordered[unlist(gene.list),]
    rownames(exp.sel) <- paste(unlist(gene.list), 1:length(unlist(gene.list)), sep="_")
    annot <- data.frame(location = label.ordered)
    rownames(annot) <- colnames(mat.ordered)
    
    annot_row <- data.frame(genes = rep(names(gene.list), times=sapply(gene.list, length)))
    rownames(annot_row) <- paste(unlist(gene.list), 1:length(unlist(gene.list)), sep="_")
    
    require(pheatmap)
    require(RColorBrewer)
    breaksList = seq(0, 20, by = 1)
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250)
    
    if (verbose==TRUE) {
            pheatmap(exp.sel,
                     cluster_cols = F, annotation_col = annot, annotation_row = annot_row,
                     cluster_rows = F, show_colnames = F, show_rownames = T, color = color)
    } else if (slim==FALSE) {
        pheatmap(exp.sel,
                 cluster_cols = F, annotation_col = annot, annotation_row = annot_row,
                 cluster_rows = F, show_colnames = F, show_rownames = F, color = color)
        } else {
        pheatmap(exp.sel,
                 cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F, 
                 legend=F, annotation_legend=F, annotation_names_col=F, annotation_names_row=F)
    }
}



make_dyno = function(sce, prop = 1L){
    ncells = ncol(sce)
    if(prop == 1){
        
    } else {
        sce = sce[,sample(ncells %>% seq_len, prop*ncells, replace = FALSE)]
    }
    wrap_expression(
        counts = counts(sce) %>% t,
        expression = logcounts(sce) %>% t,
        cell_info = colData(sce) %>% as.data.frame %>% tibble::rownames_to_column("cell_id"),
        group_info = colData(sce) %>% as.data.frame,
        feature_info = rowData(sce) %>% as.data.frame %>% tibble::rownames_to_column("feature_id")
    )
}

#adapted from http://genoweb.toulouse.inra.fr/~pmartin/pgpmartin/2018/11/14/nicer-scatterplot-in-gggally/
GGscatterPlot <- function(data, mapping, color_vec, ..., 
                          method = "spearman") {
    
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    cor <- cor(x, y, method = method)
    
    #Assemble data frame
    df <- data.frame(x = x, y = y)
    df$cols <- color_vec
    
    #Get 2D density for alpha
    dens2D <- MASS::kde2d(df$x, df$y)
    df$density <- fields::interp.surface(dens2D , 
                                         df[,c("x", "y")])
    
    if (any(df$density==0)) {
        mini2D = min(df$density[df$density!=0]) #smallest non zero value
        df$density[df$density==0] <- mini2D
    }
    
    #Prepare plot
    my_pal = colorRampPalette(viridis::magma(n=10)[9:3])(35)
    
    pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
        ggplot2::geom_point(shape=16, size=0.8, show.legend = T) +
        scale_color_gradientn(colors=my_pal) +
        ggplot2::scale_alpha(range = c(.05, .6)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    theme_classic()
    
    return(pp)
}


getStats <- function(res, method=c("Limma",
                                   "Voom",
                                   "MAST",
                                   "ttest",
                                   "EdgeR",
                                   "Wilcoxon",
                                   "DD","DP", "LR",
                                   "Bimod", "SCVI", "M3Drop", 
                                   "ROC")) {
    
    if (method=="DP") {
        
        result <- lapply(res, function(x) {
            
            stats <- x$nonzero.pvalue.adj
            names(stats) <- rownames(x)
            x$dir[is.na(x$dir)] <- FALSE
            
            stats <- sapply(1:nrow(x), function(y) {
                if(x$dir[y]==T) {stats[y]} else {1-stats[y]}
            })
            
            stats <- -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
            
        })
        return(result)
    }
    
    if (method=="DD") {
        
        result <- lapply(res, function(x) {
            stats <- x$`stats.D^+`
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method %in% c("Limma", "Voom")) {
        
        result <- lapply(res, function(x) {
            stats <- x$t
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="EdgeR") {
        
        result <- lapply(res, function(x) {
            stats <- x$F
            names(stats) <- rownames(x)
            stats <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {stats[y]} else {-stats[y]}
            })
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="ttest") {
        
        result <- lapply(res, function(x) {
            stats <- x$stats.t
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="Wilcoxon") {
        
        result <- lapply(res, function(x) {
            stats <- -x$stats.W
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="MAST") {
        
        result <- lapply(res, function(x) {
            mast <- x$pval
            names(mast) <- rownames(x)
            x$logFC[is.na(x$logFC)] <- 0
            
            stats <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {mast[y]} else {1-mast[y]}
            })
            stats <- -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            stats <- sort(stats, decreasing=T, na.last=T)
            
            return(stats)
        })
        return(result)
    }
    
    if (method %in% c("LR", "Bimod")) {
        
        result <- lapply(res, function(x) {
            stats <- x$p_val_adj
            
            logFC = grep("FC", colnames(x), value = TRUE)
            stats = sapply(1:length(x[,logFC]), function(y) {
                if (x[,logFC][[y]] > 0) { stats[[y]]} else { 1 - stats[[y]] }
            })
            names(stats) = rownames(x)
            stats <- -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            stats <- sort(stats, decreasing=TRUE, na.last=TRUE)
            
            return(stats)
        })
        return(result)
    }
    
    if (method=="ROC") {
        
        result <- lapply(res, function(x) {
            stats <- x$myAUC
            
            logFC = grep("FC", colnames(x), value = TRUE)
            stats = sapply(1:length(x[,logFC]), function(y) {
                if (x[,logFC][[y]] > 0) { stats[[y]]} else { -(1 - stats[[y]]) }
            })
            names(stats) = rownames(x)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            stats <- sort(stats, decreasing=TRUE, na.last=TRUE)
            
            return(stats)
        })
        return(result)
    }
    
    if (method=="SCVI") {
        
        result <- lapply(res, function(x) {
            stats <- log2(x$bayes_factor)
            logFC = grep("lfc_mean", colnames(x), value = TRUE)
            stats = sapply(1:length(x[,logFC]), function(y) {
                if (x[,logFC][[y]] > 0) { stats[[y]]} else { -stats[[y]] }
            })
            names(stats) = rownames(x)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            stats <- sort(stats, decreasing=TRUE, na.last=TRUE)
            
            return(stats)
        })
        return(result)
    }
    
    if (method=="M3Drop") {
        
        result <- lapply(levels(factor(res$Group)), function(x) {
            
            stats = res$adj_pval
            stats = sapply(1:nrow(res), function(y) {
                if (res[y,"Group"] ==  x) { stats[[y]]} else { 0 }
            })
            names(stats) = rownames(res)
            stats = -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            stats <- sort(stats, decreasing=TRUE, na.last=TRUE)
            return(stats)
            
        })
        names(result) = levels(factor(res$Group))
        return(result)
    }
    
}

FindDivergentGenes <- function(anchor.list, other.list, n_anchor=30, n_other=50) {
    
    other.list <- lapply(other.list, function(x) {
        res <- lapply(x, function(y) names(y)[1:n_other])
        return(res)
    })
    other.list <- lapply(1:length(other.list[[1]]), function(x) {
        unique(unlist(lapply(other.list, function(y) y[[x]])))
    })
    anchor.list <- lapply(anchor.list, function(x) {
        res <- lapply(x, function(y) names(y)[1:n_anchor])
        return(res)
    })
    anchor.unique.list <- lapply(1:length(other.list), function(x) {
        setdiff(anchor.list[["Cepo"]][[x]], other.list[[x]])
        
    })
    names(anchor.unique.list) <- names(anchor.list[["Cepo"]])
    return(anchor.unique.list)
    
}
