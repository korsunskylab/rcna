#' Perform a single diffusion mapping step.
#' 
#' @param a A sparseMatrix M by M (cells x cells) adjacency matrix (connectivities) 
#' whose m,m'-th entry indicates the similarity between cells m and m' in the graph
#' @param s A sparseMatrix M by N (cells x samples) indicator matrix assigning 
#' the m-th cell to the n-th sample from which it was drawn.
#' 
#' @note This function ignores distances and only uses unweighted connectivities
#' 
#' @returns A M by N sparseMatrix.
#' 
#' @keywords internal
diffuseStep <- function(a, s) {
    stopifnot((isa(a, "Matrix") || is.matrix(a)) && nrow(a) == ncol(a))
    stopifnot(isa(s, "Matrix") || is.matrix(s))
    stopifnot(nrow(s) == ncol(a))
    
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    rownames(res) <- rownames(s)
    colnames(res) <- colnames(s)
    return(res) 
}

#' Calculate kurtosis accounting for batch
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param batch A factor containing batch information
#' 
#' @returns a vector of kurtosis values for each column (cell) in the NAM
#' 
#' @keywords internal
batchKurtosis <- function(NAM, batch) {
    stopifnot(is.factor(batch))
    stopifnot(isa(NAM, "Matrix") || is.matrix(NAM))
    stopifnot(length(batch) == nrow(NAM))
    
    bI <- Matrix::fac2sparse(batch)
    # mean batch NAM
    bNAM <- Matrix::Diagonal(x = 1 / Matrix::rowSums(bI)) %*% bI %*% NAM
    # bNAM <- Matrix::t(bI %*% NAM) %*% Matrix::Diagonal(x = 1 / Matrix::rowSums(bI))
    res <- Kurtosis(bNAM, margin = 2)
    return(res)
}

#' Filter NAM based on batch level kurtosis
#'
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param batch A factor containing batch information or NULL. If NULL, no
#' batch QC is performed (this is for pass-through). \cr
#' Default: NULL
#' @param donors A factor containing donor information.\cr
#' Default NULL
#' @param min.threshold A numeric scalar defining the minimum value which the 
#' threshold should take.\cr
#' Default: 6
#' @param verbose Logical. Controls verbosity.\cr
#' Default: FALSE
#' 
#' @note min.threshold is 6 here, to keep functionality the same as initially 
#' written, but there is no mention of this in the paper. The rationale for the
#' choice is that 6 is "a lot" of kurtosis for a normally distributed variable.
#'
#' @returns A list containing the filtered NAM (NAM) and an index vector indicating
#' which cells were kept (keep).
#' 
#' @keywords internal
qcNAM <- function(NAM, 
                  batch = NULL, 
                  donors = NULL,
                  min.threshold = 6,
                  verbose = FALSE) {
    stopifnot(isa(NAM, "Matrix") || is.matrix(NAM))
    stopifnot(is.numeric(min.threshold) && length(min.threshold) == 1 && min.threshold >= 0)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    if (!is.null(batch)) {
        stopifnot(is.factor(batch) && length(batch) == nrow(NAM))
    }
    if (!is.null(donors)) {
        stopifnot(is.factor(donors) && length(donors) == nrow(NAM))
    }
    if (is.null(batch) && is.null(donors)) {
        if(verbose) message('Only one unique batch supplied to qc')
        keep <- rep(TRUE, ncol(NAM))
        res <- list(NAM = NAM, keep = keep)
        # TODO: donor currently ignored
    } else {
        kurtoses <- batchKurtosis(NAM, batch)
        if(verbose) message("Median batch kurtosis: ", round(stats::median(kurtoses), 3))
        threshold <- max(min.threshold, 2*stats::median(kurtoses))
        if(verbose) message("Throwing out neighborhoods with batch kurtosis >= ", 
                            round(threshold, 3))
        keep <- which(kurtoses < threshold)
        if(verbose) message("keeping ", length(keep), " neighborhoods")
        res <- list(NAM = NAM[, keep, drop = FALSE], keep = keep)
    }
    
    ### TODO: remove 0 variance columns (cells)
    # # filter samples and drop any columns that then have zero variance
    # NAM = NAM[filter_samples]
    # zero_variance_col_ix = np.where(NAM.std(axis=0) == 0)[0]
    # nz_ix = np.flatnonzero(kept)
    # kept[nz_ix[zero_variance_col_ix]] = False
    # NAM = NAM.drop(columns=NAM.columns[zero_variance_col_ix])
    
    return(res) 
}

#' Perform batch correction of a NAM using (partial) ridge regression
#'
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param cov.mat a matrix of covariates. \cr
#' Default: NULL
#' @param batch a factor containing the batch to correct. \cr
#' Default: NULL
#' @param ridges A vector of ridges to try or NULL. If NULL ridges are calculated
#' automatically centered on `svd(X)$d[1]^2`. \cr
#' Default: NULL
#' @param ridge.crit A character string specifying how the ridge regression should
#' select optimal ridges. One of c("gcv", "kurtosis"). Use "kurtosis" to match
#' original cna implementation.
#' # Default: "kurtosis"
#' @param partial A character string or NULL. If NULL, ridge regression is performed
#' on all variables. If a character, specify which set of variables ridge regression
#' should be used on, the other set does not get penalized. Must be one of "batches" or
#' "covariates". \cr 
#' Default: NULL
#' @param verbose Logical controlling verbosity.\cr
#' Default: FALSE
#'
#' @returns A list containing the corrected NAM, the Hat matrix (M), and the 
#' effective degree of freedom (`sum(diag(H))` from the ridge hat matrix) (r)
#' 
#' @keywords internal
ridgeNAM <- function(NAM, 
                     cov.mat = NULL, 
                     batch = NULL, 
                     ridges = NULL,
                     ridge.crit = c("kurtosis", "gcv"),
                     partial = NULL,
                     verbose = FALSE){
    stopifnot(isa(NAM, "Matrix") || is.matrix(NAM))
    if(!is.null(partial)){
        partial <- match.arg(partial, c("batches", "covariates"))
    }
    if(!is.null(partial) && any(c(is.null(cov.mat), is.null(batch)))){
        warning("'partial' is not null, but only 'cov.mat' or 'batch' supplied.",
                " Using standard ridge regression.")
        partial <- NULL
    }
    if(!is.null(ridges)){
        stopifnot(is.numeric(ridges) && all(ridges >= 0))
    }
    ridge.crit <- match.arg(ridge.crit, c("kurtosis", "gcv"))
    
    dfs <- 0
    N <- nrow(NAM)
    I <- Matrix::Diagonal(n = N)
    
    if(!is.null(batch) && !is.null(cov.mat)){
        stopifnot(is.factor(batch) && length(batch) == nrow(NAM))
        stopifnot((is.matrix(cov.mat) || isa(cov.mat, "Matrix")) && nrow(cov.mat) == nrow(NAM))
        
        B <- Matrix::sparse.model.matrix(~ 0 + batch)
        C <- Matrix::Matrix(cov.mat)
        X <- Scale(cbind(B, C))
        
        if(is.null(partial)){ 
            L <- Matrix::Diagonal(n = ncol(X))
        } else if(partial == "batches"){ # only do ridge on batches, fit OLS to covariates
            L <- Matrix::Diagonal(x = c(rep(1, ncol(B)), rep(0, ncol(X) - ncol(B))))
            dfs <- ncol(C)
        } else { # only do ridge on covariates, fit OLS to batches
            L <- Matrix::Diagonal(x = c(rep(0, ncol(B)), rep(1, ncol(X) - ncol(B))))
            dfs <- ncol(B)
        }
    } else if(!is.null(batch)){
        stopifnot(is.factor(batch) && length(batch) == nrow(NAM))
        X <- Matrix::sparse.model.matrix(~ 0 + batch) |> Scale()
        L <- Matrix::Diagonal(n = ncol(X))
    } else if (!is.null(cov.mat)) {
        stopifnot((is.matrix(cov.mat) || isa(cov.mat, "Matrix")) && nrow(cov.mat) == nrow(NAM))
        X <- Scale(Matrix::Matrix(cov.mat))
        L <- Matrix::Diagonal(n = ncol(X))
    } else {
        warning("no batch or cov.mat supplied.", immediate. = TRUE)
        res <- list(NAM = NAM, 
                    M = I,
                    r = 0)
        return(res)
    }
    
    if(is.null(ridges)){
        sv_max <- svd(X)$d[1]^2
        ridges <- sv_max * 10^seq(3, -3, length.out = 20)
        # ridges <- c(10^seq(-5, 5, length.out = 20), sqrt(.Machine$double.eps))
    } 
    if(ridge.crit == "gcv"){
        # find the best ridge penalty
        gcv <- vapply(ridges, function(l) {
            H <- X %*% Matrix::solve(Matrix::crossprod(X, X) + l * L, Matrix::t(X))
            resids <- (I - H) %*% NAM
            num <- (1/N) * sum(resids^2)
            den <- ((1/N) * (1 - sum(Matrix::diag(H))))^2
            num/den
        }, numeric(1)) |> round(3)
        
        lambda <- ridges[which.min(gcv)]
        H <- X %*% Matrix::solve(Matrix::crossprod(X, X) + lambda * L, Matrix::t(X))
        M <- I - H
        NAM_r <- M %*% NAM
        dfs <- dfs + sum(Matrix::diag(H))
        if(is.null(batch)){
            kurtoses <- Kurtosis(NAM_r, na.rm = T) |> stats::median()
        } else {
            kurtoses <- batchKurtosis(NAM_r, batch) |> stats::median()
        }
        if(verbose){ 
            message("with ridge ", round(lambda,3), " median batch kurtosis = ", 
                    round(kurtoses, 3))
        }
    } else { # original python implementation
        NAM_r <- NAM
        for (l in ridges) {
            H <- X %*% Matrix::solve(Matrix::crossprod(X, X) + l * L, Matrix::t(X))
            M <- I - H
            NAM_r <- M %*% NAM_r
            if(is.null(batch)){
                kurtoses <- Kurtosis(NAM_r, na.rm = T) |> stats::median()
            } else {
                kurtoses <- batchKurtosis(NAM_r, batch) |> stats::median()
            }
            if(verbose){ message("with ridge ", l, " median batch kurtosis = ", kurtoses) }
            if (kurtoses <= 6) {
                break                     
            }
        }
        dfs <- ncol(X)
    }
    
    colnames(NAM_r) <- colnames(NAM)
    rownames(NAM_r) <- rownames(NAM)
    
    return(list(NAM = NAM_r, 
                M = M,
                r = dfs))
}

#' remove batch effects (including continuous covariates)
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param covs.mat A matrix of covariates.\cr
#' Default NULL
#' @param batch A factor containing batch information.\cr
#' Default NULL
#' @param donors A factor containing donor information.\cr
#' Default NULL
#' @param partial.by A character string. One of NULL, "batches", or "covariates". 
#' Ridge penalties are only applied to those variables specified.
#' If NULL then all variables are penalized. Setting to "batches" will reproduce
#' original `rcna` behavior.\cr
#' Default: NULL
#' @param ridges A numeric vector of ridge penalties to try. If NULL, a default
#' set of ridges centered on `svd(X)$d[1]^2` is used.\cr
#' Default NULL
#' @param ridge.crit A character string specifying how the ridge regression should
#' select optimal ridges. One of c("gcv", "kurtosis"). Use "kurtosis" to match
#' original cna implementation.\cr
#' # Default: "kurtosis"
#' @param ddof An integer controlling whether the population or sample denominator
#' is used for scaling the NAM. Set to 1L to calculate R's default `n - 1`
#' denominator for the standard deviation. Default is 0L to use `n` in the denominator
#' matching pythons std() calculation. \cr
#' Default: 0L
#' 
#' @note assumes that covariates are continuous
#' @note If both batch and cov.mat are null, this is a passthrough for centering
#' and scaling the nam.
#' 
#' @returns A list containing the residuals NAM (NAM), the annihilator matrix 
#' used to remove unwanted effects (M), and the approximate degrees of freedom
#' consumed (r).
#' 
#' @keywords internal
residNAM <- function(NAM, 
                     covs.mat = NULL, 
                     batch = NULL, 
                     donors = NULL,
                     partial.by = c("batches", "covariates"),
                     ridges = NULL, 
                     ridge.crit = c("kurtosis", "gcv"),
                     ddof = 0L,
                     scale.method = "sd",
                     verbose = FALSE) {
    stopifnot(isa(NAM, "Matrix") || is.matrix(NAM))
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    if (!is.null(batch)) {
        stopifnot(is.factor(batch) && length(batch) == nrow(NAM))
    }
    if (!is.null(covs.mat)) {
        stopifnot((is.matrix(covs.mat) || isa(covs.mat, "Matrix")) && nrow(covs.mat) == nrow(NAM))
    }
    ridge.crit <- match.arg(ridge.crit, c("kurtosis", "gcv"))
    stopifnot(is.numeric(ddof) && length(ddof) == 1 && ddof %in% c(0L, 1L))
    
    NAM <- Scale(NAM, center = TRUE, scale = FALSE)
    if(!is.null(batch) || !is.null(covs.mat)){
        if(!is.null(partial.by) && (!is.null(batch) && !is.null(covs.mat))){
            partial.by <- match.arg(partial.by, choices = c("batches", "covariates"))
        } else {
            partial.by <- NULL
        }
        NAMres <- ridgeNAM(NAM = NAM, cov.mat = covs.mat, batch = batch, 
                           ridges = ridges, ridge.crit = ridge.crit,
                           partial = partial.by, verbose = verbose)
    } else {
        # pass through to generate proper results list. 
        I <- Matrix::Diagonal(n = nrow(NAM))
        colnames(I) <- rownames(I) <- rownames(NAM)
        NAMres <- list(NAM = NAM, 
                       M = I,
                       r = 0)
    }
    NAMres$NAM <- Scale(NAMres$NAM, center = FALSE, scale = TRUE, 
                        ddof = ddof, method = scale.method)
    return(NAMres)
}

#' Compute the singular-value decomposition of a NAM matrix.
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param n.pcs maximum number of components to use in the SVD. Passed to the `k`
#' parameter of [RSpectra::svds()]. If NULL or greater than `0.5 * min(dim(NAM))`
#' then the full svd is calculated. \cr
#' Default: NULL
#' 
#' @note svs are actually eigenvalues, not SVs. I squared them to be consistent 
#' with python code. 
#' @note svdNAM scales NAM columns, so this is correlation, not covariance
#' 
#' @returns A list containing the U, D, and V' components of the SVD.
#' 
#' @keywords internal
svdNAM <- function(NAM, n.pcs = NULL) {
    stopifnot(is.matrix(NAM) || isa(NAM, "Matrix"))
    
    if (is.null(n.pcs) || n.pcs > .5 * min(dim(NAM))) {
        svd_res <- svd(NAM)
        # N.B.> this is how the original python code handled it. It is the same 
        #       as just svd, but might improve performance for certain dimensions of NAM
        # svd_res <- svd(Matrix::tcrossprod(NAM))
        # svd_res$v <- Matrix::t(Matrix::t(Matrix::t(NAM) %*% svd_res$u) / sqrt(svd_res$d))
        n.pcs <- min(dim(NAM))
    } else {
        svd_res <- RSpectra::svds(NAM, k = n.pcs)
    }
    
    U <- svd_res$u[, seq_len(n.pcs)]
    colnames(U) <- paste0('PC', seq_len(n.pcs))
    rownames(U) <- rownames(NAM)
    V <- svd_res$v[, seq_len(n.pcs)]
    colnames(V) <- paste0('PC', seq_len(n.pcs))
    rownames(V) <- colnames(NAM)
    # N.B.> if using the svd of cov then don't square d
    res <- list(U = U, 
                svs = svd_res$d^2, 
                V = V)
    return(res)
}

#' Calculate Neighborhood Association Matrix (NAM) using diffusion mapping
#' 
#' @param data A rcna data object (list). See [nam()] or [association()] for details.
#' @param n.steps A numeric, the number of steps to take during the random walk. 
#' If NULL, then `kurt.delta` is checked at each step to decide whether or not 
#' to take the next step. \cr
#' Default: NULL
#' @param min.steps A numeric value controlling the minimum number of steps to take. 
#' Only used if n.steps is NULL. \cr
#' Default: 3
#' @param max.steps A numeric value controlling the maximum number of steps to take. 
#' Only used if n.steps is NULL. \cr
#' Default: 15
#' @param kurt.delta The minimum change in kurtosis between steps under which
#' the walk stops (if min.steps have already been made). Only used if n.steps is NULL.\cr
#' Default: 3
#' 
#' @returns A sparseMatrix containing the NAM (A cell by sample neighborhood
#' association matrix).
#' 
#' @keywords internal
buildNAM <- function(data, 
                     n.steps = NULL, 
                     min.steps = 3L,
                     max.steps = 15L,
                     kurt.delta = 3, 
                     # TODO: add back median kurtosis stopping parameter.
                     verbose = FALSE) {
    
    stopifnot(all(c("samplem", "obs", "connectivities", 
                    "samplem_key", "obs_key") %in% names(data)))
    stopifnot(isa(data$connectivities, "Matrix"))
    stopifnot(is.numeric(min.steps) && length(min.steps) == 1 && min.steps >= 1)
    stopifnot(is.numeric(max.steps) && length(max.steps) == 1 && max.steps >= min.steps)
    stopifnot(is.numeric(kurt.delta) && length(kurt.delta) == 1 && kurt.delta > 0)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    stopifnot(data$samplem_key %in% colnames(data$obs))
    
    # N.B.> An alternative approach is documented in Matrix::sparse.model.matrix
    #       colnames are already properly formatted as well, but requires `methods`
    #       for using `as()` to coerce
    # s <- methods::as(data$obs[,data$samplem_key], "sparseMatrix") |> Matrix::t()
    # rownames(s) <- data$obs[[data$obs_key]]
    f <- stats::as.formula(paste0("~ 0 + ", data$samplem_key))
    s <- Matrix::sparse.model.matrix(f, data$obs)
    colnames(s) <- gsub(data$samplem_key, '\\1', colnames(s))
    rownames(s) <- data$obs[[data$obs_key]]
    
    prevmedkurt <- Inf
    for (i in seq_len(max.steps)) {
        s <- diffuseStep(data$connectivities, s)
        medkurt <- propTable(s, margin = 2) |> 
            Kurtosis(margin = 1) |> 
            stats::median()
        if(is.null(n.steps)) {
            kd <- prevmedkurt - medkurt
            if(verbose) message("Median kurtosis = ", round(medkurt, 3), 
                                " Kurtosis change = ", round(kd, 3), " at step ", i)
            if (kd < kurt.delta & i > min.steps) {
                if(verbose) message("stopping after ", i, " steps")
                break 
            }
            prevmedkurt <- medkurt
        } else if (i == n.steps) {
            break
        }
    }  
    
    sfin <- propTable(s, margin = 2)
    medkurt <- stats::median(Kurtosis(sfin, margin = 1))
    if(verbose) {
        message("Final kurtosis: ", round(medkurt, 3))
    }
    snorm <- Matrix::t(sfin)
    return(snorm)
}

#' Build and decompose a neighborhood abundance matrix 
#' 
#' @param data a list containing: 
#' samplem = (sample-level metadata), 
#' obs = (cell-level metadata), 
#' connectivities  = (sparse symmetric weighted adjacency matrix),
#' samplem_key = (character string indicating the column in samplem uniquely identifying samples),
#' obs_key = (character string indicating the column in obs uniquely identifying cells),
#' N = nrow(samplem_df).
#' @param y A character string specifying the column in `data$samplem` where 
#' the variable of interest is stored.
#' @param batches A character string denoting the column in `data$samplem` where
#' batch information is stored. Only a single batch variable is allowed. \cr
#' Default: NULL 
#' @param donors A factor containing donor information.\cr
#' Default NULL
#' @param covs A character string or vector denoting the column(s) `data$samplem`
#' where continuous covariate information is stored.\cr
#' Default: NULL
#' @param n.steps Numeric scalar controlling the number of steps to take during 
#' the random walk.\cr
#' Default: NULL
#' @param min.steps Numeric scalar. Minimum number of steps to take on the random
#' walk.\cr
#' Default: 3L
#' @param max.steps Numeric scalar. Maximum number of steps to take on the random
#' walk.\cr
#' Default: 15L
#' @param kurtosis.delta Numeric scalar. Minimum kurtosis change below which the
#' random walk can stop (after travelling at least `min.steps`).\cr
#' Default: 3
#' @param min.batch.kurtosis Sets the minimum threshold below which cells are kept
#' during QC. This ensures that when kurtosis is not extreme between batches
#' then cells are not filtered. `threshold <- max(min.batch.kurtosis, 2*median(kurtoses))`.\cr
#' Default: 6
#' @param max.frac.pcs The number of PCs to calculate for SVD. The minimum is 10 
#' and the maximum is the number of samples minus 1. \cr
#' Default: 0.15
#' @param partial.by A character string. One of NULL, "batches", or "covariates". 
#' Ridge penalties are only applied to those variables specified.
#' If NULL then all variables are penalized. Setting to "batches" will reproduce
#' original cna behavior.\cr
#' Default: NULL
#' @param ridges A numeric vector of ridges to try during ridge regression. When
#' NULL ridges are calculated automatically. \cr
#' Default: NULL
#' @param ridge.crit A character string specifying how the ridge regression should
#' select optimal ridges. One of c("gcv", "kurtosis"). Use "kurtosis" to match
#' original cna implementation.\cr
#' Default: "kurtosis"
#' @param ddof An integer controlling whether the population or sample denominator
#' is used for scaling the NAM. Default is 0L which calculates python's default 
#' denominator for the standard deviation. Set to 1L to use R's default `n - 1` 
#' in the denominator. \cr
#' Default: 0L
#' @param scale.method Character string specifying the method to use for scaling
#' the data. Because [base::scale] uses root mean square, this option allows 
#' scaling by standard deviation "sd" as an alternative. Use "sd" to mimic python
#' implementation of CNA, and "rms" to use R defaults. See [Scale] for details.\cr
#' Default: "sd"
#' @param suffix A character scalar for optionally adding a suffix to the output
#' components. \cr
#' Default: ""
#' @param verbose Logical controlling the verbosity of the function.\cr
#' Default: FALSE
#' @param filter.samples STUB. Not used currently.
#' 
#' @return A named list. All keys are optionally suffixed by \code{suffix}.
#'   \item{y}{Numeric vector of the (possibly coerced) response variable.}
#'   \item{raw.NAM.T}{sparseMatrix. Transposed raw NAM (cells × samples) before
#'     QC or residualisation.}
#'   \item{qc.NAM.T}{sparseMatrix. Transposed NAM after batch-kurtosis QC.}
#'   \item{resid.NAM.T}{sparseMatrix. Transposed NAM after residualisation.}
#'   \item{M}{The annihilator matrix from [residNAM()].}
#'   \item{r}{Numeric. Degrees of freedom consumed by the correction.}
#'   \item{NAM_sampleXpc}{Matrix (samples × PCs): left singular vectors.}
#'   \item{NAM_svs}{Numeric vector of squared singular values.}
#'   \item{NAM_varexp}{Numeric vector of variance explained per PC.}
#'   \item{NAM_nbhdXpc}{Matrix (cells × PCs): right singular vectors.}
#'   \item{keptcells}{Integer index vector of retained cell columns after QC.}
#'   \item{batches}{Factor of batch labels, or \code{NULL}.}
#'   \item{donors}{Factor of donor labels, or \code{NULL}.}
#'   \item{covariates}{sparseMatrix of covariates, or \code{NA_real_} if none.}
#' 
#' @export 
nam <- function(data, y,
                batches = NULL, 
                donors = NULL,
                covs = NULL, 
                # random walk parameters passed to buildNAM
                n.steps = NULL, 
                min.steps = 3L,
                max.steps = 15L,
                kurtosis.delta = 3,
                # minimum batch kurtosis passed to qcNAM
                min.batch.kurtosis = 6,
                max.frac.pcs = 0.15, 
                # passed to residNAM
                partial.by = c("batches", "covariates"),
                ridges = NULL,
                ridge.crit = c("kurtosis", "gcv"),
                ddof = 0L,
                scale.method = c("rms", "sd"),
                # general
                suffix = '',
                filter.samples = NULL,
                verbose = FALSE) {
    if(missing(data) || missing(y)){
        stop("missing parameters with no defaults.")
    }
    stopifnot(all(c("samplem", "obs", "connectivities", 
                    "samplem_key", "obs_key", "N") %in% names(data)))
    stopifnot(is.null(batches) || (is.character(batches) && length(batches) == 1))
    stopifnot(is.null(covs) || (is.character(covs)))
    
    # if data$obs$Sample is not a factor, it will cause data ordering shenanigans
    lvls <- data$samplem[[data$samplem_key]] |> unique()
    data$obs[[data$samplem_key]] <- factor(data$obs[[data$samplem_key]], levels = lvls)
    data$samplem[[data$samplem_key]] <- factor(data$samplem[[data$samplem_key]], levels = lvls)
    
    # passed to buildNAM
    stopifnot(is.null(n.steps) || (is.numeric(n.steps) && length(n.steps) == 1))
    stopifnot(is.numeric(min.steps) && length(min.steps) == 1)
    stopifnot(is.numeric(max.steps) && length(max.steps) == 1)
    stopifnot(is.numeric(kurtosis.delta) && 
                  length(kurtosis.delta) == 1 && 
                  kurtosis.delta > 0)
    # buildNAM needs this, it is implicit and needs to be made explicit in documentation
    stopifnot(data$samplem_key %in% colnames(data$obs))
    stopifnot(!is.numeric(data$obs[,data$samplem_key]))
    
    # passed to qcNAM
    stopifnot(is.numeric(min.batch.kurtosis) && length(min.batch.kurtosis) == 1)
    stopifnot(is.numeric(max.frac.pcs) && length(max.frac.pcs) == 1)
    # passed to residNAM
    if(!is.null(partial.by)){
        partial.by <- match.arg(partial.by, c("batches", "covariates"))
    }
    stopifnot(is.null(ridges) || is.numeric(ridges))
    # general
    stopifnot(is.character(suffix) && length(suffix) == 1)
    
    if(!is.null(n.steps) && n.steps > max.steps){
        warning("'n.steps' is larger than 'max.steps', updating 'max.steps' to allow", 
                " at least 'n.steps'.", call. = FALSE, immediate. = TRUE)
        max.steps <- n.steps
    }
    ## TODO: Filtering should happen at top level, either during data object
    ##       build or as a helper to manipulate existing data object
    res <- list()
    
    if (is.null(batches)) {
        batch_vec <- NULL
    } else {
        if(!batches %in% colnames(data$samplem)){
            stop("Could not find batches: ", batches, " in `colnames(data$samplem)`")
        }
        # batches must be categorical or factor
        # TODO: handle multiple batches
        batch_vec <- dplyr::pull(data$samplem, dplyr::one_of(batches)) |> 
            as.factor()
        if(any(is.na(batch_vec))){
            stop("batches contains NAs")
        }
    }
    if(is.null(donors)){
        group_vec <- NULL
    } else {
        if(!donors %in% colnames(data$samplem)){
            stop("Could not find donors: ", donors, " in `colnames(data$samplem)`")
        }
        group_vec <- dplyr::pull(data$samplem, dplyr::one_of(donors)) |> 
            as.factor()
        if(any(is.na(group_vec))){
            stop("batches contains NAs")
        }
    }
    if (is.null(covs)) {
        cov_mat <- NULL
    } else {
        # covariates must be numeric or factors
        cov_mat <- data$samplem[, covs, drop = FALSE] 
        stopifnot(all(vapply(cov_mat, is.numeric, logical(1))))
        anycovnas <- vapply(cov_mat, \(x) { any(is.na(x)) }, logical(1))
        if(any(anycovnas)){
            stop("Some covariates contain NAs: ", covs[which(anycovnas)])
        }
        cov_mat <- Matrix::Matrix(as.matrix(cov_mat))
    }
    
    if (verbose) message('Constructing NAM')
    NAM <- buildNAM(data, 
                    n.steps = n.steps, 
                    min.steps = min.steps,
                    max.steps = max.steps, 
                    kurt.delta = kurtosis.delta, 
                    verbose = verbose)
    if (verbose) message('QC-ing NAM')
    res_qc_nam <- qcNAM(NAM, 
                        batch = batch_vec, 
                        donors = group_vec,
                        min.threshold = min.batch.kurtosis, 
                        verbose = verbose) 
    
    # TODO: handle multiple y?
    # y is returned in order to avoid recalculating NAM if not desired.
    y <- data$samplem[, y, drop = TRUE]
    if(is.character(y)){
        y <- as.numeric(as.factor(y))
    } else {
        y <- as.numeric(y)
    }
    
    ## (3) Decompose NAM 
    if (verbose) message('Residualizing NAM')
    res_resid_nam <- residNAM(res_qc_nam$NAM, 
                              covs.mat = cov_mat, 
                              batch = batch_vec, 
                              donors = group_vec,
                              partial.by = partial.by,
                              ridges = ridges,
                              ridge.crit = ridge.crit,
                              ddof = ddof,
                              scale.method = scale.method,
                              verbose = verbose)
    if (verbose) message('Decomposing NAM')
    n_pcs <- max(10, ceiling(max.frac.pcs * nrow(data$samplem)))
    n_pcs <- min(n_pcs, nrow(data$samplem) - 1) ## make sure you don't compute all SVs    
    res_svd_nam <- svdNAM(res_resid_nam$NAM, n.pcs = n_pcs)
    
    res[['y']] <- y
    res[['raw.NAM.T']] <- Matrix::t(NAM)
    res[['qc.NAM.T']] <- Matrix::t(res_qc_nam$NAM)
    res[['resid.NAM.T']] <- Matrix::t(res_resid_nam$NAM)
    res[['M']] <- res_resid_nam$M
    res[['r']] <- res_resid_nam$r
    res[['NAM_sampleXpc']] <- res_svd_nam$U
    res[['NAM_svs']] <- res_svd_nam$svs
    res[['NAM_varexp']] <- res_svd_nam$svs / nrow(res_svd_nam$U) / nrow(res_svd_nam$V)
    res[['NAM_nbhdXpc']] <- res_svd_nam$V
    res[['keptcells']] <- res_qc_nam[[2]]
    res[['batches']] <- batch_vec
    res[['donors']] <- group_vec
    res[['covariates']] <- cov_mat
    
    names(res) <- paste0(names(res), suffix)
    
    return(res)
}
