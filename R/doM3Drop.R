#' doM3Drop
#' Function to find differentially expressed genes using doM3Drop
#' 
#' @importFrom M3Drop M3DropConvertData M3DropGetMarkers
#' @param exprsMat
#' @param cellTypes vector of cell type labels
#' 
#' @examples 


doM3Drop <- function(exprsMat, cellTypes, pseudocount=1){
    
    require(ROCR)
    
    cty <- droplevels(as.factor(cellTypes))
    norm <- M3Drop::M3DropConvertData(exprsMat, 
                              is.log=TRUE, pseudocount)
    ######## Note it may internally filter undetected genes
    marker_genes <- M3Drop::M3DropGetMarkers(norm, cty)
    marker_genes$adj_pval = stats::p.adjust(marker_genes$pval, method = "BH")
    
    return(marker_genes)
}
