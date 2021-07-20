#' Extract information from Seurat object to create CNA object
#' 
#' @param seurat_object Initialized Seurat object. Assumes that FindNeighborhoods has been run.
#' @param samplem_key String denoting the name of the sample-level identifier (e.g. DonorID). 
#' @param samplem_vars Which sample-level covariates to include. 
#' @param graph_use Which graph to use. By default, will use first graph in seurat_object. 
#' 
#' @return 
#' 
#' @export 
create_object.Seurat <- function(seurat_object, samplem_key, samplem_vars, graph_use=NULL) {    
    if (length(names(seurat_object@graphs)) == 0) {
        stop('Must precompute graph in Seurat with FindNeighbors()')
    }
    if (is.null(graph_use)) {
        graph_use <- names(seurat_object@graphs)[[1]]
        message('Graph not specified. Using graph {graph_use}')
    } else {
        if (!graph_use %in% names(seurat_object@graphs)) {
            stop('Graph {graph_use} not in seurat object')
        }
    }
    samplem_vars <- c(samplem_vars, samplem_key)
    samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, one_of(samplem_vars))))
    obs_df <- tibble::rownames_to_column(seurat_object@meta.data, 'CellID')
    if (nrow(samplem_df) == nrow(obs_df)) {
        stop(
            'Sample-level metadata is same length as cell-level metadata.       
             Please check that samplem_vars are sample-level covariates.'
        )
    }

    rcna_object <- list(
        samplem = samplem_df,
        obs = obs_df, 
        connectivities = seurat_object@graphs[[graph_use]],
        samplem_key = samplem_key,
        obs_key = 'CellID',
        N = nrow(samplem_df)
    )
    return(rcna_object)
}


#' Extract information from Seurat object to create CNA object
#' 
#' @param seurat_object Initialized Seurat object. Assumes that FindNeighborhoods has been run.
#' @param test_var Contrast variable to test for association. 
#' @param samplem_key String denoting the name of the sample-level identifier (e.g. DonorID). 
#' @param graph_use Which graph to use. By default, will use first graph in seurat_object. 
#' @param batches Name of batch variable. Currently only one categorical variable allowed. 
#' @param covs Name(s) of other (numerical) covariates to control for. 
#' @param nsteps TBD
#' @param verbose TBD
#' @param assay Which seurat assay to save results under. 
#' @param key Which key to use for cached NAM PC dimensions. 
#' 
#' @return 
#' 
#' @export 
association.Seurat <- function(
    seurat_object, test_var, samplem_key, graph_use, 
    batches = NULL, covs = NULL, nsteps = NULL, verbose=TRUE, 
    assay=NULL, key='NAMPC_'
    ## TODO: put back these params: 
#     suffix = '', force_recompute = FALSE, 
) {
    ## (1) format data 
    covs_keep <- test_var
    if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
    if (!is.null(covs)) covs_keep <- c(covs_keep, covs)
    rcna_data <- create_object.Seurat(seurat_object, samplem_key, covs_keep, graph_use)
    yvals <- rcna_data$samplem[[test_var]]
    if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer') ) {
        stop('test_var is of class {class(yvals)}. It must be numeric variable for association testing.')
    }
        
    ## (2) do association
    cna_res <- association(
        data = rcna_data, 
        y = yvals, 
        batches, covs, nsteps, '', TRUE, verbose
    ) 
    
    if (is.null(assay)) {
        ## Choose first assay 
        assay <- names(seurat_object@assays)[[1]]
    }
    
    ## (3) save results 
    seurat_object[['cna']] <- Seurat::CreateDimReducObject(
        embeddings = cna_res$NAM_embeddings, 
        loadings = cna_res$NAM_loadings, 
        stdev = cna_res$NAM_svs, 
        assay = assay, 
        key = key,
        misc = cna_res ## Association results 
    )
    seurat_object@meta.data$cna_ncorrs <- cna_res$ncorrs[colnames(seurat_object), , drop=TRUE]
    seurat_object@meta.data$cna_ncorrs_fdr05 <- dplyr::case_when(
        abs(seurat_object@meta.data$cna_ncorrs) < cna_res$fdr_5p_t ~ 0,
        TRUE ~ seurat_object@meta.data$cna_ncorrs
    )
    seurat_object@meta.data$cna_ncorrs_fdr10 <- dplyr::case_when(
        abs(seurat_object@meta.data$cna_ncorrs) < cna_res$fdr_10p_t ~ 0,
        TRUE ~ seurat_object@meta.data$cna_ncorrs
    )
    
    return(seurat_object)
}
