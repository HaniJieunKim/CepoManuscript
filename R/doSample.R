doSample <- function(sce, method=c("Cepo", "Voom", "MAST", "DD"), n=100, p=0.05, workers = 2L, coi) {
    mat <- as.matrix(logcounts(sce))
    cty <- as.factor(sce$celltype)
    lapply(coi, function(celltype) {
        
        dx.rep <- mclapply(1:n, function(x) {
            set.seed(x)
            rest <- mat[,!cty %in% celltype]
            cty.rest <- as.character(cty[!cty %in% celltype])
            
            sub <- mat[, cty %in% celltype]
            sub <- sub[, sample(1:ncol(sub), round(ncol(sub)*p))]
            
            mat.sub <- cbind(rest, sub)
            cty.all <- c(cty.rest, rep(celltype, ncol(sub)))
            
            if (method == "Cepo") {
                ds.res <- Cepo(mat.sub, cty.all) 
                return(ds.res)
            } else if (method == "DD") {
                system.time(DD <- doDD(mat.sub, cty.all))
                dd.res <- getStats(DD, method="DD")
                return(dd.res)
            } else if (method == "Voom") {
                Voom <- doVoom(mat.sub, cty.all)
                voom.res <- getStats(Voom, method="Voom")
                return(voom.res)
            } else {
                system.time(MAST <- doMAST(mat.sub, cty.all))
                mast.res <- getStats(MAST, method="MAST")
                return(mast.res)
            }
        }, mc.cores=workers)
        dx.ave <- lapply(1:n, function(rep) {
            dx.rep[[rep]][[celltype]]
        })
        idx <- Reduce(intersect, lapply(dx.ave, function(x) names(x)))
        dx.ave <- do.call(cbind, lapply(dx.ave, function(x) x[idx]))
        return(dx.ave)
    })
}