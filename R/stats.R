#' Permute Y conditioned on B or G
#' 
#' @param Y A numeric vector of length n. The variable to permute.
#' @param B A factor (or object coercible via \code{split()}) of length n
#' indicating batch membership. Permutations are performed within each level of `B`.
#' @param G A factor (or object coercible via \code{split()}) of length n
#' indicating group membership. Permutations are performed within each level of `G`
#' if `B` is not specified.
#' @param num Positive integer. Number of permutations to generate. \cr
#' Default: 1000
#' @param seed An integer random seed for reproducibility. \cr
#' Default: NULL
#' @param duplicates.ok Logical. If FALSE, duplicate permutation rows are
#' removed from the output, which may result in fewer than `num` columns.\cr
#' Default: TRUE
#' 
#' @returns A numeric matrix of dimensions n × m, where `m <= num`. 
#' Each column is one permutation of `Y` within batch levels.
#'   
#' @keywords internal
conditional_permutation <- function(Y, 
                                    B = NULL, 
                                    G = NULL, 
                                    num = 1000, 
                                    seed = NULL, 
                                    duplicates.ok = TRUE) {
    stopifnot(is.numeric(Y))
    stopifnot(is.null(B) || length(B) == length(Y))
    stopifnot(is.null(G) || length(G) == length(Y))
    stopifnot(is.numeric(num) && length(num) == 1 && num >= 1)
    stopifnot(is.logical(duplicates.ok) && length(duplicates.ok) == 1)
    
    if(is.null(seed)){
        seed <- sample(1e6, 1)
    }
    set.seed(seed)
    if(!is.null(B) && !is.null(G)){
        warning("Both batch (B) and group (G) specified. Using only group permutations.",
                immediate. = TRUE)
        B <- G
    } else if(is.null(B) && !is.null(G)){
        B <- G
    } else if(is.null(B) && is.null(G)){
        B <- rep(1L, length(Y)) |> factor()
    } 
    
    y_perm <- lapply(seq_len(num), function(i) {
        res <- split(seq_len(length(Y)), B) |> 
            lapply(\(idx) {
                data.frame(idx, val = sample(Y[idx]))
            }) |> do.call(what = rbind)
        res[order(res$idx), "val", drop = TRUE]
    })
    y_perm <- matrix(unlist(y_perm), nrow = length(y_perm), byrow = TRUE)
    dup_y <- duplicated(y_perm)
    if(sum(dup_y) > 0.2*num){
        warning("Found ", sum(dup_y), " duplicated permutations.", 
                immediate. = TRUE, call. = FALSE)
    }
    if(!duplicates.ok){
        y_perm <- y_perm[!dup_y,]
    }
    return(t(y_perm))
}

#' calculate the number of coefficients within each threshold
#'
#' @param breaks A vector of thresholds 
#' @param z A matrix of coefficients 
#' @param atol A small numeric value to adjust bins to enforce corr <= bin edge.
#' It is subtracted directly from the breaks
#' @param rtol A small numeric value to adjust bins to enforce corr <= bin edge.
#' It is first multiplied by the squared coefficients before subtracting from the 
#' breaks.
#' 
#' @returns A numeric matrix of dimensions (\code{length(breaks)} x \code{ncol(z)}), 
#' where entry [i, j] is the number of rows in column j of \code{z} whose squared 
#' value exceeds \code{breaks[i]^2}.
#' 
#' @keywords internal
tail_counts <- function(breaks, z, atol = sqrt(.Machine$double.eps), rtol = 1e-5) {
    
    if (!is.matrix(z) && !isa(z, "Matrix")) {
        z <- Matrix::Matrix(z, ncol = 1)
    }
    
    n <- length(breaks)
    z2 <- z^2
    breaks2  <- breaks^2
    ix  <- order(breaks2)   
    iix <- order(ix)   
    # match python implementation
    bins <- c(breaks2[ix] - atol - rtol * breaks2[ix], Inf)

    tails <- apply(z2, 2, \(zn2){
        col <- tabulate(findInterval(zn2, bins), nbins = n)
        rev(cumsum(rev(col)))
    }) |> matrix(nrow = n) 
    
    res <- t(tails)[, iix, drop = FALSE]
    return(res)
}

#' compute the empirical FDRs via FDP procedure
#'
#' @param z A correlation matrix with a single column containing the observed 
#' correlations
#' @param znull A correlation matrix with p-permutations columns containing the 
#' permuted null correlations
#' @param thresholds a vector of bins for the distribution
#' 
#' @note The raw per-permutation FDP estimates are averaged across permutations
#' to obtain the FDR, which is then monotonised via `cummin()` to ensure
#' that the FDR curve is non-increasing as the threshold increases.
#' 
#' @returns A data.frame containing thresholds, FDRs, and number of neighborhoods 
#' in each threshold bin.
#' 
#' @keywords internal
empirical_fdrs <- function(z, znull, thresholds) {
    n <- length(thresholds) - 1
    tails <- tail_counts(thresholds, znull)
    ranks <- tail_counts(thresholds, z)
    # compute FDPs
    fdp <- sweep(tails, 2, ranks, '/')
    # make sure we didn't over-shoot upper bins too much
    # NaN occurs when there are 0 counts in both tails and ranks
    valid_thresh <- setdiff(1:ncol(fdp), unique(which(is.nan(fdp), arr.ind = TRUE)[,2]))
    ret_thresh <- thresholds[valid_thresh]
    fdp <- fdp[,valid_thresh]
    fdp[fdp > 1] <- 1
    fdr <- Matrix::colMeans(fdp)
    # enforce monotonicity
    fdr_monotone <- cummin(fdr)
    num_detected <- vapply(ret_thresh, function(.t) sum(abs(z) > .t), numeric(1))
    fdrs <- data.frame(threshold = ret_thresh,
                       fdr = fdr_monotone, 
                       num_detected = num_detected) 
    return(fdrs)
}

