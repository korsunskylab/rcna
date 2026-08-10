test_that("stats working as expected", {
    set.seed(42)
    y <- rnorm(8)
    set.seed(42)
    b <- sample(letters[1:4], 8, replace = TRUE) |> factor()
    set.seed(42)
    g <- sample(LETTERS[1:2], 8, replace = TRUE) |> factor()
    nnull <- 10
    # conditional_permutation
    conditional_permutation() |> expect_error()
    foo <- conditional_permutation(Y = y, num = nnull)
    expect_equal(dim(foo), c(length(y), nnull)) 
    
    foo <- conditional_permutation(Y = y, B = b, num = nnull)
    expect_equal(dim(foo), c(length(y), nnull))
    
    conditional_permutation(Y = y, B = b, G = g, num = nnull) |> 
        expect_warning()
    suppressWarnings( # Warning: Found 530 duplicated permutations.
        foo <- conditional_permutation(Y = y, G = g, 
                                       num = 1000, seed = 11,
                                       duplicates.ok = FALSE)
    )
    expect_equal(dim(foo), c(length(y), 470)) 
    
    # 10k cells to test
    ncells <- 10000
    set.seed(42)
    z_corrs <- matrix(runif(ncells, -1, 1), ncol = 1)
    z_nulls <- lapply(1:100, \(i) {runif(ncells, -0.9, 0.9)}) |> 
        unlist() |> matrix(ncol = 100)
    br <- seq(0, 1, length.out = 400)
    
    # each column is a break
    # each row is an "experiment"
    tcz <- tail_counts(breaks = br, z = z_corrs)
    expect_equal(dim(tcz), c(1, 400))
    tcn <- tail_counts(breaks = br, z = z_nulls)
    expect_equal(dim(tcn), c(100, 400))

    
    fdrs <- empirical_fdrs(z = z_corrs, 
                           znull = z_nulls, 
                           thresholds = br)
    expect_s3_class(fdrs, "data.frame")
    expect_equal(ncol(fdrs), 3)
    expect_equal(colnames(fdrs), c("threshold", "fdr", "num_detected"))
})

