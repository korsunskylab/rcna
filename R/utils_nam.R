## NOTE: this function ignores distances and only uses unweighted connectivities
diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(as.matrix(res)) ## remove sparsity: do we need this? 
}

.batch_kurtosis <- function(NAM, batches_vec) {
    purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
        Matrix::colMeans(NAM[i, ])
    }) %>% 
        dplyr::bind_cols() %>% 
        apply(1, moments::kurtosis)
}


.qc_nam <- function(NAM, batches=NULL) {
    N <- nrow(NAM)
    ## NOTE: added NULL check 
    if (is.null(batches) | length(unique(batches)) == 1) {
        message('only one unique batch supplied to qc')
        keep <- rep(TRUE, ncol(NAM))
        return(list(NAM = NAM, keep = keep))
    }
    kurtoses <- .batch_kurtosis(NAM, batches)
    threshold <- max(6, 2*median(kurtoses))
    message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}'))
    keep <- which(kurtoses < threshold)
    message(glue::glue('keeping {length(keep)} neighborhoods'))    
    return(list(NAM = NAM[, keep, drop = FALSE], keep = keep)) 
}

## NOTE: assumes that covariates are not categorical? 
## @param covs matrix of covariate (fixed effects?)
## @param batches vector of batch covariates (random effects?)
.resid_nam <- function(NAM, covs_mat=NULL, batches=NULL, ridge=NULL) {
    N <- nrow(NAM)
    NAM_ <- scale(NAM, center = TRUE, scale = FALSE)

    ncols_C <- 0
    if (!is.null(covs_mat)) {
        covs_mat <- scale(covs_mat)
        ncols_C <- ncols_C + ncol(covs_mat)
    }
    if (is.null(batches) | length(unique(batches)) == 1) {
        message('only one unique batch supplied to prep')
        if (is.null(covs_mat)) {
            M <- Matrix::Diagonal(n = N)
        } else {
            M <- Matrix::Diagonal(n = N) - covs_mat %*% solve(t(covs_mat) %*% covs_mat, t(covs_mat)) 
        }
        NAM_ <- M %*% NAM_ 
    } else {
        B <- model.matrix(~0+batches)
        B <- scale(B) 
        ncols_C <- ncols_C + ncol(B)
        if (is.null(covs_mat)) {
            C <- B
        } else {
            C <- cbind(B, covs_mat)
        }
        if (is.null(ridge)) {
            ridges <- c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 0)            
        } else {
            ridges <- ridge
        }
        for (ridge in ridges) {
            L <- Matrix::Diagonal(x = c(rep(1, ncol(B)), rep(0, ncol(C) - ncol(B))))  
            M <- Matrix::Diagonal(n = N) - C %*% solve(t(C) %*% C + ridge * nrow(C) * L, t(C))
            ## Why do you want to update this iteratively? 
            NAM_ <- M %*% NAM_
            kurtoses <- .batch_kurtosis(NAM_, batches)
            print(glue::glue('\twith ridge {ridge}, median batch kurtosis = {median(kurtoses)}'))
            if (median(kurtoses) <= 6) {
                break                     
            }
        }
    }
    
    return(list(
        ## NOTE: scale function in R gives slightly different results
        ##       than does similar function in python
        NAM_=scale(NAM_, center=FALSE, scale=TRUE), 
        M=M, 
        r=ncols_C
    ))
}


## NOTE: added option to pass number of PCs, for potential speedups
## NOTE: svs are actually eigenvalues, not SVs. I squared them to be consistent with python code. 
## NOTE: original implementation actually computes dot product: why? 
## NOTE: .resid_nam scales NAM columns, so this is correlation, not covariance
.svd_nam <- function(NAM, npcs) {
    if (missing(npcs) | npcs > .5 * min(dim(NAM))) {
        svd_res <- svd(NAM)
    } else {
        svd_res <- RSpectra::svds(NAM, k = npcs)
    }

    U_df <- svd_res$u[, seq_len(npcs)]
    colnames(U_df) <- paste0('PC', seq_len(npcs))
    rownames(U_df) <- rownames(NAM)
    V_df <- svd_res$v[, seq_len(npcs)]
    colnames(V_df) <- paste0('PC', seq_len(npcs))
    rownames(V_df) <- colnames(NAM)
    return(list(U=U_df, svs=svd_res$d^2, V=V_df))
}


# .nam <- function(data, nsteps=NULL, maxnsteps=15L) {
#     s <- model.matrix(~0+SampleID, data$obs)
#     colnames(s) <- gsub('^SampleID(.*)', '\\1', colnames(s))
#     rownames(s) <- data$obs$CellID
#     s <- s[, data$samplem$SampleID] ## Necessary? 
    
#     prevmedkurt <- Inf
#     ## CHECK: number of iterations matches 
#     for (i in seq_len(maxnsteps)) {
#         s <- diffuse_step(data, s)
#         medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis))
        
#         if (is.null(nsteps)) {
#             prevmedkurt <- medkurt
#             if (prevmedkurt - medkurt < 3 & i > 3) {
#                 message(glue::glue('stopping after {i} steps'))
#                 break 
#             }            
#         } else if (i == nsteps) {
#             break
#         }
#     }  
#     snorm <- t(prop.table(s, 2))
#     rownames(snorm) <- data$samplem$SampleID
#     colnames(snorm) <- data$obs$CellID
#     return(snorm)
# }
.nam <- function(data, nsteps=NULL, maxnsteps=15L) {
    f <- as.formula(as.character(glue('~0+{data$samplem_key}')))    
    s <- model.matrix(f, data$obs)
    colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
    rownames(s) <- data$obs[[data$obs_key]]    
    s <- s[, data$samplem[[data$samplem_key]]] ## Necessary? 
    
    prevmedkurt <- Inf
    ## CHECK: number of iterations matches 
    for (i in seq_len(maxnsteps)) {
        s <- diffuse_step(data, s)
        medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis))
        
        if (is.null(nsteps)) {
            prevmedkurt <- medkurt
            if (prevmedkurt - medkurt < 3 & i > 3) {
                message(glue::glue('stopping after {i} steps'))
                break 
            }            
        } else if (i == nsteps) {
            break
        }
    }  
    snorm <- t(prop.table(s, 2))
    rownames(snorm) <- data$samplem[[data$samplem_key]]
    colnames(snorm) <- data$obs[[data$obs_key]]
    return(snorm)
}



#' Build and decompose neighborhood abundance matrix 
#' 
#' @param data list containing samplem (sample-level metadata), 
#'        obs (cell-level metadata), and connectivities (dgCMatrix)
#' @param batches string(s) to denote batch variables. 
#' @param covs string(s) to denote covariate variables. 
#' @param filter_samples TBD
#' @param nsteps TBD
#' @param max_frac_pcs TBD
#' @param suffix TBD
#' @param force_recompute TBD 
#' @param verbose TBD
#' @return TBD
#' 
#' @export 
nam <- function(
    data, batches=NULL, covs=NULL, filter_samples=NULL,
    nsteps=NULL, max_frac_pcs=0.15, suffix='',
    force_recompute=FALSE, verbose=FALSE
) {
    ## TODO: check inputs! 
    ## Only one batch variables allowed
    ## covs are numeric, batches are categorical? 
    
    res <- list()
    ## TODO: add sample filtering 
#     if filter_samples is None:
#         if covs is not None:
#             filter_samples = ~np.any(np.isnan(covs), axis=1)
#         else:
#             filter_samples = np.repeat(True, data.N)
#     else:
#         filter_samples = _df_to_array(data, filter_samples)
    
    if (is.null(batches)) {
        batches_vec <- NULL
    } else {        
        ## Assumes categorical 
        batches_vec <- as.character(as.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
    }
    if (is.null(covs)) {
        covs_mat <- NULL
    } else {
        ## Assumes numeric
        covs_mat <- as.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))
    }
    
    ## (2) Compute NAM 
    ## TODO: add option to optionally not recompute 
    if (verbose) message('Construct NAM')
    NAM <- .nam(data, nsteps=nsteps)
    if (verbose) message('QC NAM')
    .res_qc_nam <- .qc_nam(NAM, batches_vec) 
    
    res[[paste0('NAM.T', suffix)]] <- t(.res_qc_nam[[1]])
    res[[paste0('keptcells', suffix)]] <- .res_qc_nam[[2]]
    res[[paste0('_batches', suffix)]] <- batches_vec

    ## (3) Decompose NAM 
    ## TODO: check if double brackets appropriate for multiple covs and/or batches
    ## TODO: check with Y&L if covs should be numerical and batches categorical
    ## NOTE: don't really need a separate SVD function, since it can be done in one line
    if (verbose) message('Residualize NAM')
    .res_resid_nam <- .resid_nam(NAM, covs_mat, batches_vec)
    res[[paste0('_M', suffix)]] <- .res_resid_nam$M
    res[[paste0('_r', suffix)]] <- .res_resid_nam$r

    if (verbose) message('Decompose NAM')
    npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))
    npcs <- min(npcs, nrow(data$samplem) - 1) ## make sure you don't compute all SVs    
    .res_svd_nam <- .svd_nam(.res_resid_nam$NAM_, npcs)
    res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
    res[[paste0('NAM_svs', suffix)]] <- .res_svd_nam$svs
    res[[paste0('NAM_varexp', suffix)]] <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
    res[[paste0('NAM_nbhdXpc', suffix)]] <- .res_svd_nam$V

    ## TODO: save covs for later
#     du['_covs'+suffix] = (np.zeros(0) if covs is None else covs)  
#     res[[paste0('_covs', suffix)]]
    return(res)
}
