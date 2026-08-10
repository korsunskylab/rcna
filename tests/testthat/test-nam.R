test_that("nam work as expected", {
    suppressPackageStartupMessages({
        library(SingleCellExperiment)
    })
    SCE <- mockSCE(add.graph = TRUE)
    
    foo <- createCNAListSCE(sce = SCE,
                            y = "Treatment", 
                            sample.key = "SampleID", 
                            sample.vars = c("Mutation_Status"),
                            knn.graph = "kNN")
    foo$samplem$Clinical <- sample(1:10, 4)
    
    nam(data = foo,
        y = "Treatment",
        verbose = FALSE) |>
        expect_type("list")
    
    # error caused by non-numeric covariates
    expect_error(
        nam(data = foo, 
            y = "Treatment", 
            covs = "Mutation_Status",
            verbose = FALSE)
    )
    
    nam(data = foo,
        y = "Treatment",
        batches = "Mutation_Status", 
        covs = "Clinical",
        verbose = FALSE) |>
        expect_type("list")
    
    # error caused by kurtosis.delta = 0
    expect_error(
        nam(data = foo, 
            y = "Treatment", 
            batches = "Mutation_Status", 
            min.steps = 1, 
            max.steps = 5, 
            kurtosis.delta = 0,
            verbose = FALSE)
    )
    
    nam(data = foo,
        y = "Treatment",
        batches = "Mutation_Status",
        min.steps = 1, 
        max.steps = 5,
        kurtosis.delta = 1,
        min.batch.kurtosis = 0,
        max.frac.pcs = 1,
        partial.by = NULL,
        ridges = 1000:1,
        verbose = FALSE) |>
        expect_type("list")
    
    suppressMessages(
        nam(data = foo, 
            y = "Treatment", 
            batches = "Mutation_Status", 
            covs = "Clinical",
            min.steps = 1, 
            max.steps = 5, 
            kurtosis.delta = 1, 
            min.batch.kurtosis = 0, 
            max.frac.pcs = 1, 
            partial.by = NULL, 
            ridges = 1:1000, 
            ridge.crit = "gcv", 
            ddof = 0, 
            scale.method = "rms", 
            suffix = "hvjgj", 
            verbose = TRUE) |> 
            expect_type("list")
    )
    
    
})