#' Fit y to the data
#' 
#' @param q A vector of values to fit. A column-scaled numeric vector of length n (samples)
#' @param k Number of components of U to fit, `1 <= k <= ncol(U)`.
#' @param U the left singular vector matrix of SVD (samples × PCs).
#' 
#' @returns a list containing fitted values (qhat) and regression coefficients (beta) 
#'
#' @keywords internal
regress <- function(q, k, U) {
    stopifnot(is.numeric(q) || isa(q, "Matrix"))
    stopifnot(is.numeric(k) && length(k) == 1 && k >= 1 && k <= ncol(U))
    stopifnot(is.matrix(U) || isa(U, "Matrix"))
    stopifnot(length(q) == nrow(U))
    
    Xpc <- U[, 1:k]
    beta <- crossprod(Xpc, q)
    qhat <- Xpc %*% beta
    res <- list(qhat = qhat, beta = beta)
    return(res)
}

#' Calculate the F-test p-value
#' 
#' @param yhat y fit to the data. Must be a numeric vectors of length n.
#' @param ycond y corrected for batch and covariates. Must be a numeric vectors 
#' of length n.
#' @param k Number of components of U to test.
#' @param n The number of samples.
#' @param r degrees of freedom consumed by batch correction in [residNAM()].
#' Default: FALSE
#' 
#' @returns A list containing p-value and r^2 for the F-test
#' 
#' @keywords internal
calcStats <- function(yhat, ycond, k, n, r) {
    stopifnot(is.numeric(k) && length(k) == 1 && k >= 1)
    stopifnot(is.numeric(n) && length(n) == 1 && n > 1 + r + k)
    stopifnot(is.numeric(r) && length(r) == 1 && r >= 0)
    
    ssefull <- crossprod(yhat - ycond)
    ssered <- crossprod(ycond)
    deltasse <-  ssered - ssefull
    f <- (deltasse / k) / (ssefull / n)
    f <- as.vector(f)
    p <- stats::pf(f, k, n-(1+r+k), lower.tail = FALSE)
    r2 <- as.vector(1 - ssefull/ssered)
    res <- list(p = p, r2 = r2)
    return(res)
}

#' Calculate the p-value for all ks and return the minimum
#' 
#' @param M The annihilator matrix (M returned by [residNAM])
#' @param y A scaled numeric vector containing the variable of interest
#' @param ks A vector of ks to test.
#' @param U a matrix whose columns contain the left singular vectors of x 
#' (U returned by [svdNAM])
#' @param r the degrees of freedom already consumed by batch correction 
#' (r returned by [residNAM])
#' Default: FALSE
#' 
#' @returns A list containing the k, p, and r^2 for k which produces the minimum
#' p-value.
#' 
#' @keywords internal
minpStats <- function(M, y, ks, U, r) {
    stopifnot(is.matrix(M) || isa(M, "Matrix"))
    stopifnot(is.numeric(y))
    stopifnot(nrow(M) == ncol(M) && nrow(M) == length(y))
    stopifnot(is.numeric(ks) && all(ks >= 1) && all(ks <= ncol(U)))
    
    n <- length(y)
    zcond <- M %*% y
    zcond <- Scale(zcond, center = FALSE, scale = TRUE, method = "sd")
    qhats <- lapply(ks, function(k) regress(zcond, k, U)$qhat)
    .tmp <- lapply(seq_along(ks), function(i) calcStats(qhats[[i]], zcond, ks[i], n, r))
    ps <- vapply(.tmp, \(x) x$p, numeric(1))
    r2s <- vapply(.tmp, \(x) x$r2, numeric(1))
    k_ <- which.min(ps)
    res <- list(k = ks[k_], 
                p = ps[k_], 
                r2 = r2s[k_])
    return(res)
}

#' Calculate global and (optionally) local association tests.
#' 
#' @param NAMres The list output from the [nam()] function.
#' @param ks A numeric scalar selecting the number of components of the SVD
#' to test for global association. If null 4 values for k are selected between
#' n/50 and n/5. \cr
#' Default: NULL
#' @param Nnull A numeric specifying the number of null permutations.\cr
#' Default: 1000
#' @param force_permute_all A logical controlling whether permutations of `y` 
#' should preserve batch information. \cr
#' Default: FALSE 
#' @param allow_duplicate_perms Logical. If `FALSE`, duplicate permutations
#' of `y` are removed from the null distribution. Useful for small n where
#' the number of distinct permutations is limited. Default is set to `TRUE` to
#' reproduce original behavior, but should probably be set to `FALSE` for numerical
#' accuracy. \cr
#' Default: TRUE.
#' @param local_test A logical, whether or not to perform local test of neighborhood
#' correlations. \cr
#' Default: TRUE
#' @param seed A numeric seed to set. Set if you want to repeat exact permutations
#' as a previous run. If NULL, a random seed is chosen. \cr
#' Default: NULL
#' @param verbose A logical controlling verbosity.\cr
#' Default: FALSE
#'
#' @returns A named list with the following elements:
#'   \item{p}{Empirical global association p-value from permutation test.}
#'   \item{minps_null}{Numeric vector of length \code{Nnull} containing the
#'     minimum permutation p-values from the null distribution.}
#'   \item{k}{Integer. The number of NAM PCs selected by the minimum p-value
#'     criterion.}
#'   \item{ncorrs}{Matrix of neighborhood-level association scores, rows
#'     corresponding to neighborhoods.}
#'   \item{fdrs}{Data frame of FDR results from [empirical_fdrs()],
#'     or \code{NULL} if \code{local_test = FALSE}. Columns: \code{threshold},
#'     \code{fdr}, \code{num_detected}.}
#'   \item{fdr_5p_t}{Numeric scalar. Minimum threshold achieving FDR < 5\%,
#'     or \code{NULL} if none exists.}
#'   \item{fdr_10p_t}{Numeric scalar. Minimum threshold achieving FDR < 10\%,
#'     or \code{NULL} if none exists.}
#'   \item{yhat}{Fitted values of \code{y} from the selected k-component model.}
#'   \item{ycond}{Batch/covariate-conditioned and scaled \code{y} vector.}
#'   \item{ks}{Numeric vector of k values that were tested.}
#'   \item{beta}{Regression coefficients for the selected k-component model.}
#'   \item{r2}{R-squared for the selected model.}
#'   \item{r2_perpc}{Per-PC contribution to R-squared, length \code{k}.}
#'   \item{nullr2_mean}{Mean R-squared across null permutations.}
#'   \item{nullr2_std}{Standard deviation of R-squared across null permutations.}
#'   \item{seed}{The integer seed used for permutations, for reproducibility.}
#' 
#' @keywords internal
innerAssociation <- function(NAMres, 
                             ks = NULL, 
                             Nnull = 1000, 
                             force_permute_all = FALSE, 
                             allow_duplicate_perms = TRUE,
                             local_test = TRUE, 
                             seed = NULL, 
                             verbose = FALSE) {
    stopifnot(is.list(NAMres))
    # N.B.> these are only required if supplied catch it internally?
    # , "batches", "donors"
    stopifnot(all(c("NAM_sampleXpc", "NAM_svs", "NAM_nbhdXpc", 
                    "M", "r", "y", "resid.NAM.T") %in% names(NAMres)))
    y <- Scale(matrix(NAMres[["y"]], ncol = 1), ddof = 0, method = "sd")

    stopifnot(is.numeric(y))
    stopifnot(is.numeric(Nnull) && length(Nnull) == 1 && Nnull >= 1)
    stopifnot(is.logical(force_permute_all) && length(force_permute_all) == 1)
    stopifnot(is.logical(allow_duplicate_perms) && length(allow_duplicate_perms) == 1)
    stopifnot(is.logical(local_test) && length(local_test) == 1)
    
    if (is.null(seed)) {
        seed <- sample(1e6, 1)
    } else {
        stopifnot(is.numeric(seed) && length(seed) == 1)
    }
    set.seed(seed)
    
    # prep data
    U <- NAMres[["NAM_sampleXpc"]]
    sv <- NAMres[["NAM_svs"]]
    V <- NAMres[["NAM_nbhdXpc"]]
    M <- NAMres[["M"]] 
    r <- NAMres[["r"]]
    n <- length(y)
    NAM_resid <- Matrix::t(NAMres[["resid.NAM.T"]])
    
    batches_vec <- NAMres[["batches"]]
    group_vec <- NAMres[["donors"]]
    
    if (force_permute_all) {
        batches_vec <- rep(1L, length(y)) |> factor()
    } 
    
    if (is.null(ks)) {
        # unique handles situations where n < 16
        # N.B.: now checks all k instead of just 4 evenly spaced ks
        # ks <- seq(n/50, n/5, by = 1) |> ceiling() |> unique()
        # TODO: ask Yakir for rationale for this
        incr <- max(as.integer(0.02 * n), 1)
        maxnpcs <- max(min(4 * incr, as.integer(n / 5)), 1)
        ks <- seq(from = incr, to = maxnpcs, by = incr)
    } else if(length(ks) == 1){
        if(ks < dim(U)[2]){
            ks <- seq(1, ks, 1)
        } else {
            warning("user defined ks is too large. setting it to maximum allowable: ", 
                    dim(U)[2] - 1)
            ks <- seq(1, dim(U)[2] - 1, 1)
        }
    } 
    
    # get non-null f-test p-value
    mp <- minpStats(M, y, ks, U, r)
    k <- mp$k 
    p <- mp$p 
    r2 <- mp$r2
    if (k == max(ks)) {
        warning("Data supported use of ", k, 
                " NAM PCs, which is the maximum considered.", 
                " Consider allowing more PCs by using the 'ks' argument.", 
                call. = FALSE)
    }
    
    # compute coefficients and r2 with chosen model
    ycond <- M %*% y
    ycond <- Scale(ycond, center = FALSE, scale = TRUE, method = "sd")
    rg <- regress(ycond, k, U)
    yhat <- rg$qhat
    beta <- rg$beta
    r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))^2
    r2_perpc <- as.vector(r2_perpc)
    
    # TODO: something downstream doesn't like the sparse matrix hence not a Matrix here
    ncorrs <- Matrix::colMeans(y[,1] * NAM_resid) |> matrix(ncol = 1)
    rownames(ncorrs) <- colnames(NAM_resid)
    
    # compute final p-value using Nnull null f-test p-values
    # TODO: ensure group level is working.
    y_null <- conditional_permutation(Y = y, 
                                      B = batches_vec, 
                                      G = group_vec,
                                      num = Nnull, 
                                      duplicates.ok = allow_duplicate_perms, 
                                      seed = seed)
    .tmp <- apply(y_null, 2, \(z) minpStats(M, z, ks, U, r))
    # return(list(ynull = y_null, M = M, ks = ks, U = U, r = r, mps = .tmp))
    minps_null <- vapply(.tmp, \(x) x$p, numeric(1))
    nullr2s <- vapply(.tmp, \(x) x$r2, numeric(1))
    # add sqrt(.Machine$double.eps) for floating point maths and add 1 to avoid pfinal = 0
    pfinal <- (sum(minps_null <= p + sqrt(.Machine$double.eps)) + 1)/(Nnull + 1)
    if (sum(minps_null <= p + sqrt(.Machine$double.eps)) == 0) {
        warning("global association p-value attained minimal possible value.", 
                "Consider increasing Nnull", immediate. = TRUE, call. = FALSE)
    }
    
    fdr_5p_t <- NULL 
    fdr_10p_t <- NULL
    # get neighborhood fdrs if requested
    if (local_test) {
        if(verbose) message('computing neighborhood-level FDRs')
        Nnull <- min(1000, ncol(y_null))
        y_null <- y_null[, 1:Nnull]
        ycond_null <- M %*% y_null
        ycond_null <- Scale(ycond_null, center = FALSE, scale = TRUE, method = "sd")
        # TODO: Something downstream doesn't like sparse matrices
        #   Error in max(c(abs(ncorrs), abs(ncorrs_null))) : 
        #           invalid 'type' (list) of argument
        ncorrs_null <- abs(Matrix::crossprod(NAM_resid, ycond_null) / n) |> 
            as.matrix()
        
        # use both ncorrs and ncorrs_null to make sure we observe full range
        # maxcorr <- max(max(abs(ncorrs)), 0.001)
        maxcorr <- max(c(abs(ncorrs), abs(ncorrs_null), 0.001)) + sqrt(.Machine$double.eps)
        n.bins <- 400
        fdr_thresholds <- seq(from = maxcorr/4, maxcorr, maxcorr/n.bins)
        fdrs <- empirical_fdrs(z = ncorrs, 
                               znull = ncorrs_null, 
                               thresholds = fdr_thresholds) 
        
        # find maximal FDR<5% and FDR<10% sets
        if (min(fdrs$fdr) > 0.05) {
            fdr_5p_t <- NULL
        } else {
            fdr_5p_t <- min(subset(fdrs, fdr <= 0.05)$threshold)
        }
        if (min(fdrs$fdr) > 0.10) {
            fdr_10p_t <- NULL
        } else {
            fdr_10p_t <- min(subset(fdrs, fdr <= 0.1)$threshold)
        }
    }
    
    res <- list(p = pfinal, 
                minps_null = minps_null,
                k = k,
                ncorrs = ncorrs, 
                ncorrs_null = ncorrs_null,
                fdrs = fdrs,
                fdr_5p_t = fdr_5p_t, 
                fdr_10p_t = fdr_10p_t,
                yhat = yhat, 
                ycond = ycond,
                ks = ks, 
                beta = beta,
                r2 = r2, 
                r2_perpc = r2_perpc,
                nullr2_mean = mean(nullr2s), 
                nullr2_std = stats::sd(nullr2s),
                seed = seed)

    return(res)
}

#' Main function to perform CNA association analysis
#' 
#' @description The core function of CNA, builds the NAM, and tests y for 
#' association with neighborhoods.
#' 
#' @param data a list
#' @param nam.result Optional. A pre-computed result list from [nam()]. If
#' supplied, the NAM construction step is skipped. If both `data` and
#' `nam.result` are supplied, `nam.result` is used with a warning.
#' @param y A character string specifying the column in `samplem` containing
#' the variable of interest.
#' @param batches A character string specifying the column in `samplem` 
#' containing the batch variable. Only a single batch variable is currently supported. \cr
#' Default: NULL
#' @param covs A character string or vector specifying the column(s) in `samplem` 
#' containing the covariate variables.\cr
#' Default: NULL
#' @param donors A character sting specyfing the column in `samplem` containing
#' a grouping variable. If batch is NULL, then Null permutations are permuted based
#' on donors. \cr
#' Default: NULL
#' @param n.steps Number of steps to take during the random walk. If specified then
#' exactly this many steps is taken on the random walk. \cr
#' Default: NULL
#' @param suffix A character string to be appended to the results. Useful
#'  when running multiple [association()] calls on the same data.\cr
#' Default: ''
#' @param return.nam Logical controlling whether or not to return the NAM. \cr
#' Default: FALSE
#' @param filter.samples STUB. Currently ignored.\cr
#' @param Ks One of NULL, a numeric scalar, or a numeric vector. If NULL, then 
#' `seq(n/50, n/5, length.out = 4)` are checked. If a numeric scalar then 1:Ks are
#' checked. If a numeric vector, then all supplied values are checked. \cr
#' Default: NULL
#' @param N.nulls A positive integer. Number of null permutations used to
#' construct the empirical null distribution for the global association test.\cr
#' Default: 1000
#' @param force.permute.all A logical controlling whether permutations of `y` 
#' should preserve batch information. \cr
#' Default: FALSE 
#' @param allow.duplicate.perms Logical. If FALSE, duplicate permutations
#' of `y` are removed from the null distribution. Useful for small n where
#' the number of distinct permutations is limited. Default is set to TRUE to
#' reproduce original behavior, but should probably be set to FALSE for numerical
#' accuracy. \cr
#' Default: TRUE.
#' @param local.test Logical. If TRUE, neighbourhood-level association
#' scores and empirical FDRs are computed in addition to the global test.\cr
#' Default: TRUE.
#' @param seed An integer seed for the permutation RNG, or NULL for a
#' randomly selected seed. In either case, the seed is returned in the result to 
#' allow reproducible results if needed.\cr
#' Default: NULL
#' @param verbose Logical controlling the verbosity of the function.\cr
#' Default: FALSE
#' @param ... Additional parameters passed to \code{\link{nam}}
#' 
#' @details The list object that cna expects is:
#'  * `samplem` A data.frame containing sample-level metadata. Must have a single 
#'  row for each sample in the data.#'      
#'  * `obs` A data.frame containing cell-level metadata. Must have a single row 
#'  for each cell in the data. Must have at least 2 columns, one for cell id and
#'  one for sample id. #'      
#'  * `connectivities` A symmetric weighted adjacency \code{\link[Matrix]{sparseMatrix}}. 
#'  Likely generated during UMAP calculation.#'   
#'  * `samplem_key` A character string indicating the column in \code{samplem} 
#'  uniquely identifying samples.
#'  * `obs_key` A character string indicating the column in obs uniquely identifying 
#'  cells.#'   
#'  * `N` A numeric scalar indicating the number of samples i.e. nrow(samplem).
#'  
#' @return A named list with the following elements:
#'   \item{p}{
#'      Empirical global association p-value from permutation test.
#'   }
#'   \item{minps_null}{
#'      Numeric vector of length \code{Nnull} containing the
#'      minimum permutation p-values from the null distribution.
#'   }
#'   \item{k}{
#'      Integer. The number of NAM PCs selected by the minimum p-value
#'      criterion.
#'   }
#'   \item{ncorrs}{
#'      Matrix of neighborhood-level association scores, rows
#'      corresponding to neighborhoods.
#'   }
#'   \item{fdrs}{
#'      Data frame of FDR results from \code{\link{empirical_fdrs}}, or NULL if 
#'      \code{local_test = FALSE}. Columns: \code{threshold}, \code{fdr}, \code{num_detected}.
#'   }
#'   \item{fdr_5p_t}{
#'      Numeric scalar. Minimum threshold achieving FDR < 5\%, or NULL if none exists.
#'   }
#'   \item{fdr_10p_t}{
#'      Numeric scalar. Minimum threshold achieving FDR < 10\%, or NULL if none exists.
#'   }
#'   \item{yhat}{
#'      Fitted values of \code{y} from the selected k-component model.
#'   }
#'   \item{ycond}{
#'      Batch/covariate-conditioned and scaled \code{y} vector.
#'   }
#'   \item{ks}{
#'      Numeric vector of k values that were tested.
#'   }
#'   \item{beta}{
#'      Regression coefficients for the selected k-component model.
#'   }
#'   \item{r2}{
#'      R-squared for the selected model.
#'   }
#'   \item{r2_perpc}{
#'      Per-PC contribution to R-squared, length \code{k}.
#'   }
#'   \item{nullr2_mean}{
#'      Mean R-squared across null permutations.
#'   }
#'   \item{nullr2_std}{
#'      Standard deviation of R-squared across null permutations.
#'   }
#'   \item{seed}{
#'      The integer seed used for permutations, for reproducibility.
#'   }
#'   \item{NAM_embeddings}{
#'      The Neighborhood x PCs matrix (The Left singular vectors returned by 
#'      \code{\link[base]{svd}}). If \code{return.nam = TRUE}
#'   }
#'   \item{NAM_loadings}{
#'      The Sample x PCs matrix (The Right singular vectors returned by 
#'      \code{\link[base]{svd}}). If \code{return.nam = TRUE}
#'   }
#'   \item{NAM_svs}{
#'      The _squared_ singular values returned by \code{\link[base]{svd}}. 
#'      If \code{return.nam = TRUE}
#'   }
#' 
#' @export 
association <- function(data, nam.result, y, 
                        batches = NULL, 
                        donors = NULL,
                        covs = NULL, 
                        n.steps = NULL, 
                        suffix = '',
                        return.nam = TRUE, 
                        filter.samples = NULL,
                        # passed to inner association
                        Ks = NULL, 
                        N.nulls = 1000, 
                        force.permute.all = FALSE, 
                        allow.duplicate.perms = TRUE,
                        local.test = TRUE, 
                        seed = NULL,
                        verbose = TRUE, 
                        ...) {
    if(missing(data) && missing(nam.result)){
        stop("One of 'data' or 'nam.result' must be specified.")
    } else if(!missing(data) && !missing(nam.result)){
        warning("both 'data' and 'nam.result' were supplied,", 
                " using nam.result for association test.", 
                immediate. = TRUE, call. = FALSE)
    }
    if(!missing(data) && missing(nam.result)){
        stopifnot(all(c("samplem", "obs", "connectivities", 
                        "samplem_key", "obs_key", "N") %in% names(data)))
        stopifnot(is.character(y) && length(y) == 1)
        if(!is.null(batches)){
            stopifnot(batches %in% colnames(data$samplem))
        }
        if(!is.null(donors)){
            stopifnot(donors %in% colnames(data$samplem))
        }
        if(!is.null(covs)){
            stopifnot(all(covs %in% colnames(data$samplem)))
        }
        # buildNAM needs this, it is implicit and needs to be made explicit in documentation
        stopifnot(data$samplem_key %in% colnames(data$obs))
        stopifnot(!is.numeric(data$obs[,data$samplem_key]))
    }
    
    stopifnot(is.null(batches) || (length(batches) == 1  && is.character(batches)))
    stopifnot(is.null(donors) || (length(donors) == 1  && is.character(donors)))
    
    if(!is.null(n.steps)){
        stopifnot(length(n.steps) == 1 && is.numeric(n.steps))
    }
    stopifnot(length(suffix) == 1 && is.character(suffix))
    stopifnot(length(return.nam) == 1 && is.logical(return.nam))
    stopifnot(length(N.nulls) == 1 && is.numeric(N.nulls))
    
    stopifnot(length(force.permute.all) == 1 && is.logical(force.permute.all))
    stopifnot(length(allow.duplicate.perms) == 1 && is.logical(allow.duplicate.perms))
    stopifnot(length(local.test) == 1 && is.logical(local.test))
    stopifnot(length(verbose) == 1 && is.logical(verbose))
    
    if(!is.null(Ks)){
        stopifnot(is.numeric(Ks))
    }
    
    if (verbose) message('Build NAM PCs')
    if(missing(nam.result)){
        nam_res <- nam(data, y = y, 
                       batches = batches, 
                       donors = donors,
                       covs = covs, 
                       filter.samples = filter.samples,
                       n.steps = n.steps, 
                       suffix = suffix, 
                       verbose = verbose,
                       ...)
    } else {
        nam_res <- nam.result
    }

    if (verbose) message('Perform association testing')
    res <- innerAssociation(NAMres = nam_res,
                            ks = Ks, 
                            Nnull = N.nulls, 
                            force_permute_all = force.permute.all, 
                            allow_duplicate_perms = allow.duplicate.perms, 
                            local_test = local.test, 
                            seed = seed,
                            verbose = verbose)
    if (return.nam) {
        res <- c(res, nam_res)
    }
    
    return(res)
}

