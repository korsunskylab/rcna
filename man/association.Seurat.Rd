% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_Seurat.R
\name{association.Seurat}
\alias{association.Seurat}
\title{Extract information from Seurat object to create CNA object}
\usage{
association.Seurat(
  seurat_object,
  test_var,
  samplem_key,
  graph_use,
  batches = NULL,
  covs = NULL,
  nsteps = NULL,
  verbose = TRUE,
  assay = NULL,
  key = "NAMPC_"
)
}
\arguments{
\item{seurat_object}{Initialized Seurat object. Assumes that FindNeighborhoods has been run.}

\item{test_var}{Contrast variable to test for association.}

\item{samplem_key}{String denoting the name of the sample-level identifier (e.g. DonorID).}

\item{graph_use}{Which graph to use. By default, will use first graph in seurat_object.}

\item{batches}{Name of batch variable. Currently only one categorical variable allowed.}

\item{covs}{Name(s) of other (numerical) covariates to control for.}

\item{nsteps}{TBD}

\item{verbose}{TBD}

\item{assay}{Which seurat assay to save results under.}

\item{key}{Which key to use for cached NAM PC dimensions.}
}
\value{

}
\description{
Extract information from Seurat object to create CNA object
}
