se_flow <- function(se, batch_size, by = "column", generator, random = TRUE, ...) {
    next_batch <- 0
    Ind <- seq(ncol(se))
    function() {
        next_batch <<- next_batch + 1
        if(by == "column") {            
            nbatch <- ceiling(ncol(se) / batch_size)
            if (next_batch > nbatch) {
                next_batch <<- 1
            }
            ## random start
            if(random && next_batch == 1) {
                sind <- sample(Ind, 1)
                if(sind > 1) {
                    Ind <<- c(seq(sind, ncol(se)), seq(sind - 1))
                }
            }
            
            if(next_batch == nbatch){
                idx  <- seq((next_batch - 1) * batch_size + 1,
                            ncol(se))
            }else{
                idx  <- seq((next_batch - 1) * batch_size + 1,
                            next_batch * batch_size)
            }
            dat <- se[, Ind[idx]]
        }
        generator(dat, ...)
    }
}

tpm <- function(counts, len) {
    x <- (counts + 2) / len
    log2(t(t(x)*1e6/colSums(x)))
}

zscore <- function(x) (x - mean(x)) / sd(x)

Norm <- function(rse_gene, fgenes, count, version=FALSE) {
    if(count == "rpkm"){
        lgC <- edgeR::rpkm(assays(rse_gene)$counts,
                              gene.length = rowData(rse_gene)$bp_length,
                              log = TRUE)
        ## zs <- t(apply(lgCnt1, 2, zscore))
    }else if(count == "cpm"){
        lgC <- edgeR::cpm(assays(rse_gene)$counts,
                           log = TRUE)
    }else if(count == "tpm"){
        lgC <- tpm(as.matrix(assays(rse_gene)$counts), len = rowData(rse_gene)$bp_length)
    }
    zs <- t(apply(lgC, 2, zscore))
    if(!version){
        fgenes <- sub("\\..*", "", fgenes)
        colnames(zs) <- sub("\\..*", "", colnames(zs))
    }
    zs <- zs[,match(fgenes, colnames(zs))]
    zs[is.na(zs)] <- 0
    return(zs)
}

genFun <- function(se, fgenes, count = "rpkm") {
    x <- Norm(se, fgenes, count)
    ts <- colData(se)$SMTS
    y <- model.matrix(~ ts + 0)
    colnames(y) <- sub("^ts", "", colnames(y))
    list(x = x, y = y)
}

##
ggFun <- function(se) {
    norm_se <- genFun(se)
    lgCnt_z <- t(norm_se$x)
    
    gmlist <- vector("list", ncol(se))
    for(i in seq(ncol(se))){
        gpm1 <- data.frame(gg2, exp=lgCnt_z[match(gg2[,1], rownames(lgCnt_z)), i])
        gm1 <- gpm1 %>% spread(X2, exp, fill=0) %>% dplyr::select(-X1) %>% as.matrix
        ## resampling
        r <- raster(gm1)
        extent(r) <- extent(c(-180, 180, -90, 90))
        s <- raster(nrow=64, ncol=64)
        s <- resample(r,s)
        gm1 <- as.matrix(s)
        gmlist[[i]] <- array(gm1, dim=c(dim(gm1),1))
    }
    ggmatrix <- abind(gmlist, along=0)
    list(x = ggmatrix, y = norm_se$y)
}
