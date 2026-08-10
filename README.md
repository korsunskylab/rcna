This package provides a native R implementation of the (`cna` python package)[https://github.com/immunogenomics/cna] by Reshef and Rumker. 

To get started, install `rcna`

```{r}
library(devtools)
install_github('korsunskylab/rcna')
```
and follow the (pre-rendedred vignette)[docs/rcna_demo.html]. You can also re-render the vignette yourself by installing the required packages and running
```{r}
# install rcna requirements and vignette requirements
devtools::install(dependencies = TRUE) 
# render the vignette
rmarkdown::render("vignettes/rcna_demo.Rmd", output_dir="docs")
```

Next, 
- Read the (paper)[https://www.nature.com/articles/s41587-021-01066-4].
- Read the (python code)[https://github.com/immunogenomics/cna].

# Release Notes
## August 2026
Version `0.1.0` removes a small set of Seurat-specific convenience helpers to reduce maintenance complexity, but `rcna` remains fully compatible with Seurat workflows via conversion to `SingleCellExperiment` objects.
```{r}
# convert Seurat to SCE
sce <- SeuratObject::as.SingleCellExperiment(obj)
SingleCellExperiment::colPair(sce, "connectivities") <- obj@graphs[["RNA_nn"]]

# Call rcna's analysis function(s) — replace `rcna_analysis` with the function you use
results <- rcna::associationSCE(sce=sce, y=...)
```