Generator <- function(sce, phenoData, Num.mixtures = 50, pool.size = 1000, min.percentage = 1, max.percentage = 99, seed = 24){ 
    
    CT = levels(droplevels(as.factor(phenoData$cellType)))
    ?stopifnot(length(CT) >= 2)
    
    set.seed(seed)
    require(dplyr)
    require(gtools)
    
    cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE)
    cell.distribution = cell.distribution[cell.distribution$Var1 %in% CT,]
    colnames(cell.distribution) = c("CT","max.n")
    
    Tissues = list()
    Proportions = list()
    
    for(y in 1:Num.mixtures){
        
        print(y)
        #Only allow feasible mixtures based on cell distribution
        while(!exists("P")){
            
            num.CT.mixture = sample(x = 2:length(CT),1)
            selected.CT = sample(CT, num.CT.mixture, replace = FALSE)
            
            P = runif(num.CT.mixture, min.percentage, max.percentage) 
            P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
            P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)
            
            missing.CT = CT[!CT %in% selected.CT]
            missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
            
            P = rbind.data.frame(P, missing.CT)
            potential.mix = merge(P, cell.distribution)
            potential.mix$size = potential.mix$expected * pool.size
            
            if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
                rm(list="P") 
            }
            
        }
        
        # Using info in P to build T simultaneously
        chosen_cells <- sapply(which(P$expected != 0), function(x){
            
            n.cells = P$expected[x] * pool.size
            chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
                            n.cells)
            
            chosen
        }) %>% unlist()
        
        
        T <- Matrix::rowSums(sce[,colnames(sce) %in% chosen_cells]) %>% as.data.frame()
        colnames(T) = paste("mix",y,sep="")
        
        P = P[,c("CT","expected")]
        P$mix = paste("mix",y,sep="")
        
        Tissues[[y]] <- T
        Proportions[[y]] <- P
        
        rm(list=c("T","P","chosen_cells","missing.CT"))
        
    }
    
    P = do.call(rbind.data.frame, Proportions)
    T = do.call(cbind.data.frame, Tissues)
    
    P = data.table::dcast(P, CT ~ mix, 
                          value.var = "expected",
                          fun.aggregate = sum) %>% data.frame(.,row.names = 1) 
    
    P = P[,gtools::mixedsort(colnames(P))]
    
    return(list(T = T, P = P))
    
} 
