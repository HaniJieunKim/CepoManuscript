rm_cols_rows_zeroes = function(x){
    x[rowSums(x) != 0, colSums(x) != 0]
}
nz_mean = function(x){mean(x[x != 0])}
nz_sd = function(x){sd(x[x != 0])}
nz_length = function(x){length(x[x != 0])}
nz_zprop = function(x, zprop){sample(1:length(x[x != 0]), round(length(x[x != 0])*zprop/100))}
nz_zeros = function(x, zprop){sum(x != 0)/length(x)}


range_transform = function(x, new_max = 0.99, new_min = 0.01){
    old_min = min(x, na.rm = TRUE)
    old_max = max(x, na.rm = TRUE)
    old_range = old_max - old_min
    new_range = (new_max - new_min)
    new_x = (((x - old_min) * new_range) / old_range) + new_min
    return(new_x)
}

thres_prob = function(x, thres = 0.5){ifelse(x < thres, 0, x)}

#' @param x a vector 
#' @param labels celltype labels, a character vector 
#' @param target_celltype an element of \code{target_celltype}, which is the celltype that DS algorithm will be simulated for
#' @param mean_scale numeric, should be strictly less than 1 to scale down the mean expression of other celltypes
#' @param sd_scale numeric, should be strictly great than 1 to scale down up sd of expression of other celltypes
#' @example 
#' set.seed(123)
#' n = 100
#' p = 10
#' mat = matrix(rpois(n*p, lambda = 10), nrow = p, ncol = n)
#' rownames(mat) = paste0("gene_", 1:p)
#' colnames(mat) = paste0("cell_", 1:n)
#' labels = sample(letters[1:3], size = ncol(mat), replace = TRUE)
#' replace_ds_genes(x = mat[1,], labels = labels)
replace_ds_genes = function(x, labels, target_celltype = sort(unique(labels))[1], mean_scale = 0.01, sd_scale = 2){
    stopifnot(is.vector(x))
    index_a = which(labels == target_celltype) ## All index of targeted celltype
    index_b = which(labels != target_celltype)
    
    xa = x[index_a]
    xb = x[index_b]
    
    xb_nz_replace = rnorm(n = nz_length(xb),
                          mean = mean_scale*nz_mean(xb), sd = sd_scale*nz_sd(xb))
    xb_replace = xb
    xb_replace[xb_replace != 0] = xb_nz_replace
    xb_replace = ifelse(xb_replace <= 0, 0, xb_replace)
    xb_replace = ifelse(xb_replace >= nz_mean(xa), 0, xb_replace)
    
    result = x
    result[index_b] = xb_replace
    return(result) 
}

replace_ds_genes_v2 = function(x, labels, target_celltype = sort(unique(labels))[1], mean_scale = 1, sd_scale = 2, add_zero = 0){
    
    stopifnot(is.vector(x))
    
    index_a = which(labels == target_celltype) ## All index of targeted celltype
    index_b = which(labels != target_celltype)
    
    xa = x[index_a]
    xb = x[index_b]
    
    
    xb_nz_replace = rnorm(n = nz_length(xb), 
                          mean = mean_scale*nz_mean(xb), sd = sd_scale*nz_sd(xb))
    
    xb_replace = xb
    xb_replace[xb_replace != 0] = xb_nz_replace 
    xb_replace = ifelse(xb_replace <= 0, 0, xb_replace)

    if(add_zero > 0){
        size_non_z = length(xb_replace[xb_replace != 0])
        non_z_idx = which(xb_replace != 0)
        xb_rand_zeros = sample(non_z_idx, size = round(add_zero/100*size_non_z))
        xb_replace[xb_rand_zeros] = 0
    }
    
    result = x
    result[index_b] = xb_replace
    return(result) 
}



compute_rank_auc = function(rank_tbl){
    rank_tbl %>% 
        dplyr::select(gene_name, ds, contains("rank")) %>% 
        tidyr::pivot_longer(cols = contains("rank"), 
                            names_to = "rank_type", 
                            values_to = "rank_value") %>% 
        group_by(rank_type) %>% 
        roc_auc(truth = factor(ds), rank_value)
}

compute_rank_roc = function(rank_tbl){
    rank_tbl %>% 
        dplyr::select(gene_name, ds, contains("rank")) %>% 
        tidyr::pivot_longer(cols = contains("rank"), 
                            names_to = "rank_type", 
                            values_to = "rank_value") %>% 
        group_by(rank_type) %>% 
        roc_curve(truth = factor(ds), rank_value)
}