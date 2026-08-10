test_that("plots work as expected", {
    suppressPackageStartupMessages({
        library(ggplot2)
    })
    expect_s3_class(
        dimplotCoords(coords = data.frame(foo = 1:10, bar = 10:1),
                      y = runif(10)),
        "ggplot"
    )
    
    SCE <- mockSCE(add.graph = TRUE)
    suppressWarnings(
        SCE <- associationSCE(sce = SCE, 
                              y = "Treatment", 
                              sample.key = "SampleID", 
                              graph = "kNN")
    )
    
    expect_s3_class(
        dimplotCorr(sce = SCE),
        "ggplot"
    )
    suppressWarnings({
        foo <- createCNAListSCE(sce = SCE,
                                y = "Treatment",
                                sample.key = "SampleID",
                                knn.graph = "kNN")
        res <- association(foo, y = "Treatment", verbose = FALSE)
    })
    
    expect_s3_class(
        dimplotCorr(res = res),
        "ggplot"
    )
    
    expect_error(
        dimplotCorr(sce = sce, res = res)
    )
    
})