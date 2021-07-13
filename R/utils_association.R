.association <- function(
    NAMsvd, M, r, y, batches_vec, 
    ks=NULL, Nnull=1000, 
    force_permute_all=FALSE, local_test=TRUE, seed=NULL
) {
    if (is.null(seed)) {
        set.seed(sample(1e6, 1))
    }
    if (force_permute_all) {
        batches_vec <- rep(1L, length(y))
    }
    
    # prep data
    U <- NAMsvd[[1]]
    sv <- NAMsvd[[2]]
    V <- NAMsvd[[3]]
    y <- scale(y)
    n <- length(y)

    if (is.null(ks)) {
        incr <- max(round(0.02*n), 1)
        maxnpcs <- min(4*incr, round(n/5))
        ks <- seq(incr, maxnpcs+1, incr)
    }

    .reg <- function(q, k) {
        Xpc <- U[, 1:k]
        beta <- t(Xpc) %*% q
        qhat <- Xpc %*% beta
        return(list(qhat = qhat, beta = beta))
    }

    .stats <- function(yhat, ycond, k) {
        ssefull <- crossprod(yhat - ycond)
        ssered <- crossprod(ycond)
        deltasse <-  ssered - ssefull
        f <- (deltasse / k) / (ssefull/n)
        p <- -pf(f, k, n-(1+r+k), log.p = TRUE)    
        r2 <- 1 - ssefull/ssered
        return(list(p=p, r2=r2))
    }

    .minp_stats <- function(z) {
        zcond <- M %*% z
        zcond <- scale(zcond, center = FALSE, scale = TRUE)
        qhats <- purrr::map(ks, function(k) .reg(zcond, k)$qhat)
        .tmp <- purrr::map2(qhats, ks, function(qhat, k) .stats(qhat, zcond, k))
        ps <- purrr::map_dbl(.tmp, 'p')
        r2s <- purrr::map_dbl(.tmp, 'r2')
        k_ <- which.min(ps)
        return(list(k=ks[k_], p=ps[k_], r2=r2s[k_]))
    }


    # get non-null f-test p-value
    .tmp <- .minp_stats(y)
    k <- .tmp$k; p <- .tmp$p; r2 <- .tmp$r2
    if (k == max(ks)) {
        warning(glue::glue('data supported use of {k} NAM PCs, which is the maximum considered. Consider allowing more PCs by using the "ks" argument.'))        
    }


    # compute coefficients and r2 with chosen model
    ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
    .tmp <- .reg(ycond, k)
    yhat <- .tmp$qhat; beta <- .tmp$beta
    r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))**2

    # get neighborhood scores with chosen model
    ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta/n)    

    # compute final p-value using Nnull null f-test p-values
    y_ <- conditional_permutation(batches_vec, y, Nnull)
    .tmp <- apply(y_, 2, .minp_stats)
    nullminps <- purrr::map_dbl(.tmp, 'p')
    nullr2s <- purrr::map_dbl(.tmp, 'r2')

    pfinal <- (sum(nullminps <= p+1e-8) + 1)/(Nnull + 1)
    if (sum(nullminps <= p+1e-8) == 0) {
        warning('global association p-value attained minimal possible value. Consider increasing Nnull')
    }

    # get neighborhood fdrs if requested
    fdrs <- NULL
    fdr_5p_t <- NULL 
    fdr_10p_t <- NULL
    if (local_test) {    
        message('computing neighborhood-level FDRs')
        Nnull <- min(1000, Nnull)
        y_ <- y_[, 1:Nnull]
        ycond_ <- scale(M %*% y_, center = FALSE, scale = TRUE)
        gamma_ <- crossprod(U[, 1:k], ycond_)
        nullncorrs <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_ / n)))

        maxcorr <- max(abs(ncorrs))
        fdr_thresholds <- seq(maxcorr/4, maxcorr, maxcorr/400)
        fdr_vals <- empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
        fdrs = data.frame(
    #         threshold = fdr_thresholds
            threshold = head(fdr_thresholds, -1),
            fdr = fdr_vals, 
            num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t)) 
        )
        # find maximal FDR<5% and FDR<10% sets
        if (min(fdrs$fdr) > 0.05) {        
            fdr_5p_t <- NULL
        } else {
            fdr_5p_t <- min(subset(fdrs, fdr < 0.1)$threshold)        
        }
        if (min(fdrs$fdr) > 0.05) {        
            fdr_10p_t <- NULL
        } else {
            fdr_10p_t <- min(subset(fdrs, fdr < 0.1)$threshold)
        }
    }


    res <- list(
        p = pfinal, 
        nullminps=nullminps,
        k=k,
        ncorrs=ncorrs, 
        fdrs=fdrs,
        fdr_5p_t=fdr_5p_t, 
        fdr_10p_t=fdr_10p_t,
        yhat=yhat, 
        ycond=ycond,
        ks=ks, 
        beta=beta,
        r2=r2, 
        r2_perpc=r2_perpc,
        nullr2_mean=mean(nullr2s), 
        nullr2_std=sd(nullr2s)
    )

    return(res)
}

            
#' Main function to perform CNA association 
#' 
#' @param data list containing samplem (sample-level metadata), 
#'        obs (cell-level metadata), and connectivities (dgCMatrix)
#' @param y vector with contrast variable value for association 
#' @param batches string(s) to denote batch variables. 
#' @param covs string(s) to denote covariate variables. 
#' @param nsteps TBD
#' @param suffix TBD
#' @param force_recompute TBD 
#' @param verbose TBD
#' @return TBD
#' 
#' @export 
association <- function(
    data, y, batches=NULL, covs=NULL, nsteps=NULL, suffix='',
    force_recompute=FALSE, verbose=TRUE
) {

    
    # formatting and error checking
    ## For association, batches needs to be a numeric vector
    if (is.null(batches)) {
        batches_vec <- rep(1, data$N)
    } else {
        batches_vec <- as.integer(as.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
    }
    
    ## CHECK: need _df_to_array in R? 
#     covs = _df_to_array(data, covs)
#     batches = _df_to_array(data, batches)
#     y = _df_to_array(data, y)
    
    ## TODO: check lengths
#     if y.shape != (data.N,):
#         raise ValueError(
#             'y should be an array of length data.N; instead its shape is: '+str(y.shape))

    ## TODO: add sample filtering 
#     if covs is not None:
#         filter_samples = ~(np.isnan(y) | np.any(np.isnan(covs), axis=1))
#     else:
#         filter_samples = ~np.isnan(y)

    ## Here, data has all the du things 
#     du = data.uns
    if (verbose) message('Build NAM PCs')
    nam_res <- nam(data, batches=batches, covs=covs, filter_samples=filter_samples,
                    nsteps=nsteps, suffix=suffix,
                    force_recompute=force_recompute)

    if (verbose) message('Perform association testing')
    res <- .association(
        NAMsvd=list(
            nam_res$NAM_sampleXpc,
            nam_res$NAM_svs,
            nam_res$NAM_nbhdXpc
        ),
        nam_res[[paste0('_M', suffix)]],
        nam_res[[paste0('_r', suffix)]],
        ## TODO: add filter_samples to nam results
#         y[nam_res[[paste0('_filter_samples', suffix)]]],
#         batches[nam_res[paste0('_filter_samples', suffix)]] 
        y, 
        batches_vec
    )

    # TODO: add info about kept cells
#     vars(res)['kept'] = du['keptcells'+suffix]

    return(res)
}
                             