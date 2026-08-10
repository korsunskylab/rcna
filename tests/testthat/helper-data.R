make_rcna_data <- function(n_samples = 20, n_cells = 200, seed = 42, y.sd = 1) {
    set.seed(seed)
    sample_ids <- paste0("S", seq_len(n_samples))
    cell_ids   <- paste0("C", seq_len(n_cells))
    
    b <- rep(c("A", "B"), each = floor(n_samples / 2))
    d <- rep(paste0("D", 1:2), each = floor(n_samples / 2))
    
    if(n_samples %% 2 == 1){
        b <- factor(c(b, "B"))
        d <- factor(c(d, "D1"))
    }
    
    samplem <- data.frame(
        SampleID  = sample_ids,
        y         = rnorm(n_samples, sd = y.sd),
        batch     = b,
        donor     = d,
        covariate = rnorm(n_samples),
        stringsAsFactors = FALSE
    )
    
    obs <- data.frame(
        CellID   = cell_ids,
        SampleID = sample(sample_ids, n_cells, replace = TRUE),
        stringsAsFactors = FALSE
    )
    
    # sparse symmetric adjacency (connectivities)
    adj <- Matrix::rsparsematrix(n_cells, n_cells, density = 0.05, symmetric = TRUE)
    adj@x <- abs(adj@x)  # non-negative weights
    diag(adj) <- 0
    
    res <- list(
        samplem        = samplem,
        obs            = obs,
        connectivities = adj,
        samplem_key    = "SampleID",
        obs_key        = "CellID",
        N              = n_samples
    )
    return(res)
}

# make_rcna_sce <- function(n_samples = 20, n_cells = 200, seed = 42, y.sd = 1) {
#     set.seed(seed)
#     sample_ids <- paste0("S", seq_len(n_samples))
#     cell_ids   <- paste0("C", seq_len(n_cells))
#     
#     b <- rep(c("A", "B"), each = floor(n_samples / 2))
#     d <- rep(paste0("D", 1:2), each = floor(n_samples / 2))
#     
#     if(n_samples %% 2 == 1){
#         b <- factor(c(b, "B"))
#         d <- factor(c(d, "D1"))
#     }
#     
#     samplem <- data.frame(
#         SampleID  = sample_ids,
#         y         = rnorm(n_samples, sd = y.sd),
#         batch     = b,
#         donor     = d,
#         covariate = rnorm(n_samples),
#         stringsAsFactors = FALSE
#     )
#     
#     obs <- data.frame(
#         CellID   = cell_ids,
#         SampleID = sample(sample_ids, n_cells, replace = TRUE),
#         stringsAsFactors = FALSE
#     )
#     
#     # sparse symmetric adjacency (connectivities)
#     adj <- Matrix::rsparsematrix(n_cells, n_cells, density = 0.05, symmetric = TRUE)
#     adj@x <- abs(adj@x)  # non-negative weights
#     diag(adj) <- 0
#     
#     # res <- list(
#     #     samplem        = samplem,
#     #     obs            = obs,
#     #     connectivities = adj,
#     #     samplem_key    = "SampleID",
#     #     obs_key        = "CellID",
#     #     N              = n_samples
#     # )
#     return(res)
# }

mockSCE <- function(n.cells = 200,
                    n.genes = 2000,
                    use.sparse = TRUE,
                    transpose = FALSE, 
                    add.graph = FALSE, 
                    graph.name = "kNN") {
    # TODO: for certain dimensions force sparse matrix?
    if(n.genes > 900000){
        stop("cannot generate a SingleCellExperiment with > 900000 genes")
    }
    cell.means <- 2^stats::runif(n.genes, 2, 10)
    cell.disp <- 100 / cell.means + 0.5
    data <- stats::rnbinom(n.genes * n.cells,
                           mu = cell.means,
                           size = 1 / cell.disp)
    if (use.sparse) {
        cell.data <- Matrix::Matrix(data, 
                                    ncol = n.cells,
                                    sparse = TRUE)
    } else {
        cell.data <- matrix(data, ncol = n.cells)
    }
    rownames(cell.data) <- sprintf("Gene_%s", formatC(seq_len(n.genes), width = 5, flag = 0))
    colnames(cell.data) <- sprintf("Cell_%s", formatC(seq_len(n.cells), width = 3, flag = 0))
    cdat <- data.frame(
        SampleID = NA_character_, # just like having it as the first column
        Mutation_Status = sample(c("positive", "negative"), n.cells, replace = TRUE),
        Cell_Cycle = sample(c("S", "G0", "G1", "G2M"), n.cells, replace = TRUE),
        Treatment = sample(c("treat1", "treat2"), n.cells, replace = TRUE),
        Celltype = sample(c(
            "Macrophage", "B Cell", "Monocyte", "T Cell",
            "Neutrophil", "Dendritic Cell"
        ), n.cells, replace = TRUE)
    )
    cdat$SampleID <- paste0("sample", as.numeric(factor(paste0(cdat$Mutation_Status, "_", 
                                                               cdat$Treatment))))
    rdat <- data.frame(
        gene_id = paste0("ENSG", "00000", sample(100000:999999, n.genes, replace = FALSE)),
        gene_type = sample(c("snoRNA", "pseudogene", "miRNA", "lncRNA", "protein_coding"),
                           n.genes, replace = TRUE),
        chr = sample(paste0("chr", c(seq_len(22), "X", "Y", "MT")), n.genes, replace = TRUE)
    )
    if(!transpose){
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts = cell.data),
                                                          colData = cdat,
                                                          rowData = rdat)
    } else {
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Matrix::t(cell.data)),
                                                          colData = rdat,
                                                          rowData = cdat)
    }
    if(add.graph){
        hits <- S4Vectors::SelfHits(
            sample(ncol(sce), floor(ncol(sce)*0.1)),
            sample(ncol(sce), floor(ncol(sce)*0.1)),
            nnode = ncol(sce)
        )
        S4Vectors::mcols(hits)$value <- runif(floor(ncol(sce)*0.1))
        SingleCellExperiment::colPair(sce, graph.name) <- hits
    }
    
    return(sce)
}
