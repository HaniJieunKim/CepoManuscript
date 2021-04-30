deconvolution = function(T, C, P = P, method = "NNLS") {
    
    if (method=="DeconRNAseq") {
        # "DeconRNASeq"){ #nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog) datasets and reference matrix: signatures, need to be non-negative. "use.scale": whether the data should be centered or scaled, default = TRUE 
        unloadNamespace("Seurat") #needed for PCA step
        library(pcaMethods) #needed for DeconRNASeq to work
        RESULTS = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(T), 
                                             signatures = as.data.frame(C), 
                                             proportions = NULL, checksig = FALSE, 
                                             known.prop = FALSE, use.scale = FALSE, 
                                             fig = FALSE)$out.all)
        colnames(RESULTS) = colnames(T)
        require(Seurat)
        
    } else if (method=="OLS"){
        
        RESULTS = apply(T,2,function(x) lm(x ~ as.matrix(C))$coefficients[-1])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
        
    } else if (method=="NNLS"){
        
        require(nnls)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- colnames(C)
        print("done")
    } else if (method=="FARDEEP"){
        
        require(FARDEEP)
        RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        
    } else if (method=="RLR"){ #RLR = robust linear regression
        
        require(MASS)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
        
    } else if (method=="DCQ"){#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
        
        require(ComICS)
        RESULTS = t(ComICS::dcq(reference_data = C, mix_data = T, marker_set = as.data.frame(row.names(C)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        
    } else if (method=="EN"){#standardize = TRUE by default. lambda=NULL by default 
        
        require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        
    } else if (method=="RIDGE"){ #alpha=0
        
        require(glmnet)
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        
    } else if (method=="LASSO"){ #alpha=1; shrinking some coefficients to 0. 
        
        require(glmnet)
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        
        RESULTS[is.na(RESULTS)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
        
    }
    
    RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)),]
    RESULTS = data.table::melt(RESULTS)
    colnames(RESULTS) <-c("CT","tissue","observed_values")
    
    P = P[gtools::mixedsort(rownames(P)),]
    P$CT = rownames(P)
    P = data.table::melt(P, id.vars="CT")
    colnames(P) <-c("CT","tissue","expected_values")
    
    RESULTS = merge(RESULTS,P)
    RESULTS$expected_values <-round(RESULTS$expected_values,3)
    RESULTS$observed_values <-round(RESULTS$observed_values,3)
    
    return(RESULTS)
}