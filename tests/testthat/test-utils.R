test_that("utils working as expected", {
    suppressPackageStartupMessages({
        library(SingleCellExperiment)
        library(Matrix)
    })

    SCE <- mockSCE(add.graph = TRUE)
    expect_error(
        createCNAListSCE(
            sce = SCE,
            y = "group",
            sample.key = "SampleID"
        ),
        "Parameters missing"
    )
    
    bad_graph <- Matrix::Matrix(0, nrow = 3, ncol = 3, sparse = TRUE)
    expect_error(
        createCNAListSCE(
            sce = SCE,
            y = "group",
            sample.key = "SampleID",
            knn.graph = bad_graph
        ),
        "'knn.graph' must be a square sparseMatrix with dims equal to 'ncol\\(sce\\)'"
    )
    SCE$sample_id <- SCE$SampleID[order(SCE$SampleID)]
    expect_error(
        createCNAListSCE(
            sce = SCE,
            y = "Treatment",
            sample.key = "sample_id",
            knn.graph = "kNN"
        ),
        "must collapse to a unique value"
    )
    SCE$sample_id <- NULL
    
    res <- createCNAListSCE(
        sce = SCE,
        y = "Treatment",
        sample.key = "SampleID",
        knn.graph = "kNN"
    )
    expect_type(res, "list")
    expect_named(res, c("samplem", "obs", "connectivities", "samplem_key", "obs_key", "N"))
    expect_equal(res$samplem_key, "SampleID")
    expect_equal(res$obs_key, "CellID")
    expect_equal(nrow(res$obs), ncol(SCE))
    expect_equal(res$N, length(unique(SCE$SampleID)))
    
    expect_warning( 
        expect_warning(
            resSCE <- associationSCE(sce = SCE, 
                                     y = "Treatment", 
                                     sample.key = "SampleID", 
                                     graph = "kNN"),
            "Data supported use of 1 NAM PCs"
        ),
        "Found 994 duplicated permutations"
    )
    
    expect_contains(colnames(colData(resSCE)), 
                    c("cna_ncorrs", 
                      "cna_ncorrs_fdr05", 
                      "cna_ncorrs_fdr10"))
    expect_contains(reducedDimNames(resSCE), "CNA")
    expect_contains(names(metadata(resSCE)), "CNA")
    expect_contains(names(metadata(resSCE)$CNA), 
                    c("p", "minps_null", "k", "ncorrs_null", "fdrs", "fdr_5p_t", 
                      "fdr_10p_t", "yhat", "ycond", "ks", "beta", "r2", "r2_perpc", 
                      "nullr2_mean", "nullr2_std", "seed", "y", "raw.NAM.T", "qc.NAM.T", 
                      "resid.NAM.T", "M", "r", "NAM_varexp", "keptcells"))
    
})