conditional_permutation <- function(B, Y, num) {
    purrr::map(seq_len(num), function(i) {
        split(seq_len(length(Y)), B) %>% purrr::map(function(idx) {
            data.frame(idx, val=sample(Y[idx]))
        }) %>% dplyr::bind_rows() %>% 
            dplyr::arrange(idx) %>% 
            with(val)    
    }) %>% 
        purrr::reduce(Matrix::cbind2)
}


## ~2x slower than python version 
tail_counts <- function(z, znull) {
    apply(znull, 2, function(znulli) {
        as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))        
    })
}

empirical_fdrs <- function(z, znull, thresholds) {    
    n <- length(thresholds) - 1
    tails <- t(tail_counts(thresholds, znull)[1:n, ])
    ranks <- t(tail_counts(thresholds, z)[1:n, ])

    # compute FDPs
    fdp <- sweep(tails, 2, ranks, '/')
    fdr <- Matrix::colMeans(fdp)

    return(fdr)
}


