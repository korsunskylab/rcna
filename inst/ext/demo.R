#
# # # 
#
rm(list = ls())
gc()
# ctrl shift f10
devtools::document()
library(tictoc)
# # # for hand running tests...skip if you are going to run tests automagically
library(testthat)
source("tests/testthat/helper-data.R")

ltf <- list.files("tests/testthat/", pattern = "*.R")
tic()
for(i in ltf){
    tic(paste0("Started testing: ", i))
    testthat::test_file(paste0("tests/testthat/", i))
    toc()
}
toc()

testthat::test_local()

# Run tests and generate a report 17min
tic()
# options("covr.fix_parallel_mcexit" = TRUE)
covr::report()
toc()


devtools::check()
devtools::build(manual = TRUE)


# 
# # # check that we can replicate the workflow from the anndata
#
library(anndataR)
library(SingleCellExperiment)
library(ggplot2)
sce_obj <- read_h5ad("~/Git/cna-rcna-test/test.h5ad", as = "SingleCellExperiment")

reducedDim(sce_obj, "X_umap") |> 
    cbind(colData(sce_obj)[, c("case_coef", "case_coef_fdr")]) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = case_coef)) +
    ggplot2::scale_color_gradient2(low = "darkblue", 
                                   high = "darkred") + 
    theme_classic()

reducedDim(sce_obj, "X_umap") |> 
    cbind(colData(sce_obj)[, c("male_coef", "male_coef_fdr")]) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = male_coef)) +
    ggplot2::scale_color_gradient2(low = "darkblue", 
                                   high = "darkred") + 
    theme_classic()

reducedDim(sce_obj, "X_umap") |> as.data.frame() |> 
    cbind(leiden = as.factor(sce_obj$leiden)) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = leiden)) +
    theme_classic()

# Not easy to do gradient fill for violins
reducedDim(sce_obj, "X_umap") |> as.data.frame() |> 
    cbind(leiden = as.factor(sce_obj$leiden), 
          male_coef = sce_obj$male_coef) |> 
    ggplot(aes(x = leiden, y = male_coef, color = leiden)) + 
    geom_violin() +
    theme_classic()

# 
# # # get the association object including nam and connectivities and save to an rds 
#
library(reticulate)
library(HSSRscripts)
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test.py")


connectivities <- d$obsp["connectivities"]
connectivities <- as(connectivities, "CsparseMatrix")

obs <- d$obs |> rownames_to_column("cell_id") |> 
    dplyr::mutate(id = as.factor(id), case = as.factor(case), 
                  male = as.factor(male), batch = as.factor(batch))

samplem <- obs |> dplyr::distinct(id, .keep_all = T) |> 
    dplyr::mutate(cell_id = NULL, case_coef = NULL, case_coef_fdr = NULL)

rownames(samplem) <- paste0("sample_", samplem$id)

demo <- list(
    connectivities = connectivities,
    samplem = samplem,
    samplem_key = "id",
    obs = obs,
    obs_key = "cell_id",
    N = nrow(samplem)
)

saveRDS(demo, "inst/ext/demo_data.rds")

p <- builtins$vars(p)
saveRDS(p, "inst/ext/demo_result_nocov_nobatch.rds")




# 
# # # now let's see if rcna recapitulates p from connectivities
#
devtools::document(".")
rcna_demo <- readRDS("inst/ext/demo_data.rds")

a <- nam(data = rcna_demo, y = "case")
b <- association(rcna_demo, nam.result = a, y = "case")

p <- readRDS("inst/ext/demo_result_nocov_nobatch.rds")

# only attributes differ
all.equal(as(as.matrix(p$M), "CsparseMatrix"), 
          as(as(a$`_M`, "CsparseMatrix"), "dsCMatrix"))

p$namresid |> as.matrix() |> as("sparseMatrix") |> dim()
Matrix::t(a$resid.NAM.T) |> dim()

all.equal(as(as.matrix(p$namresid), "CsparseMatrix"), 
          as(Matrix::t(a$resid.NAM.T), "CsparseMatrix"))
# [1] TRUE

#################################
# because things seem to diverge quickly, going through step by step :/
##############
library(reticulate)
library(HSSRscripts)
devtools::document(".")
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test.py")
rcna_demo <- readRDS("inst/ext/demo_data.rds")

#
### checking diffusion step...
#
diffuse_step1 <- d1 |> as.matrix() |> as("CsparseMatrix")


f <- as.formula(paste0("~ 0 + ", rcna_demo$samplem_key))
s <- Matrix::sparse.model.matrix(f, rcna_demo$obs)
colnames(s) <- gsub(rcna_demo$samplem_key, '\\1', colnames(s))
rownames(s) <- rcna_demo$obs[[rcna_demo$obs_key]]
diffuseStep <- function(a, s) {
    stopifnot((isa(a, "Matrix") || is.matrix(a)) && nrow(a) == ncol(a))
    stopifnot(isa(s, "Matrix") || is.matrix(s))
    stopifnot(nrow(s) == ncol(a))
    
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(res) 
}

rdiffuse_s1 <- s <- diffuseStep(rcna_demo$connectivities, s)

all.equal(rdiffuse_s1, diffuse_step1)

rdiffuse_s2 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step2 <- d2 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s2, diffuse_step2)

s <- diffuseStep(rcna_demo$connectivities, s)
rdiffuse_s4 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step4 <- d4 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s4, diffuse_step4)

s <- diffuseStep(rcna_demo$connectivities, s)
s <- diffuseStep(rcna_demo$connectivities, s)
s <- diffuseStep(rcna_demo$connectivities, s)
rdiffuse_s8 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step8 <- d8 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s8, diffuse_step8)

#### diffusion steps match exactly.

#
### checking nam step...both took 4 steps
#
pysnorm <- as.matrix(snorm) |> as("CsparseMatrix")

rsnorm <- buildNAM(rcna_demo)

all.equal(pysnorm, rsnorm)
# [1] TRUE

#### nam steps match exactly.

#
### checking qcnam step...
#

# NO BATCH
pyqc_nam <- qc_nam |> as.matrix() |> as("CsparseMatrix")
rqcnam <- buildNAM(rcna_demo) |> qcNAM()

all.equal(pyqc_nam, rqcnam$NAM)
# TRUE

# qc works for both batch and batch+covs

#
### checking svdnam step...
#

rnam <- buildNAM(rcna_demo)
resid_rnam <- residNAM(rnam)
rsvd_nam <- svdNAM(resid_rnam$NAM)
all.equal(abs(as.matrix(svd_nam[[1]])), abs(as.matrix(rsvd_nam[[1]])))
# TRUE (after much painful bug finding)

#### svd steps match


#########################
# NAM processing matches exactly
###

#
### checking association step...
#
# rnam <- nam(rcna_demo, y = "case")
# 
# 
# r1 <- sweep(rnam$resid.NAM.T, 2, scale(rnam$y), "*") |>  Matrix::colMeans()
# r2 <- Matrix::colMeans(rnam$resid.NAM.T) * scale(rnam$y)[,1]
# r3 <- Matrix::colMeans(scale(rnam$y) * Matrix::t(rnam$resid.NAM.T)) |> matrix(ncol = 1)
# all.equal(r1, r2)
# V <- rnam$NAM_nbhdXpc
# sv <- rnam$NAM_svs
# ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta / n)

rass <- association(rcna_demo, y = "case")

# global p is identical
identical(rass$p, p$p)
# TRUE

identical(rass$k, p$k)
# TRUE

all.equal(as.numeric(rass$yhat), as.numeric(p$yresid_hat))
# TRUE

all.equal(as.numeric(rass$ycond), as.numeric(p$yresid))
# TRUE

all.equal(as.numeric(rass$ncorrs), as.numeric(p$ncorrs))
# "Mean relative difference: 0.001598857"

all.equal(as.numeric(rass$r2_perpc), as.numeric(p$r2_perpc))
# TRUE

all.equal(p$fdr_5p_t, rass$fdr_5p_t)
# "Mean relative difference: 1.563773e-08"
# basically .Machine$double.eps difference?
# sqrt(.Machine$double.eps) = 1.490116e-08

data.frame(py.thresh = p$fdrs$threshold, r.thresh = rass$fdrs$threshold[100:399])

plot(p$fdrs$threshold, rass$fdrs$threshold[100:399])
abline(0,1)

plot(p$fdrs$fdr, rass$fdrs$fdr[100:399])
abline(0,1)

plot(p$fdrs$num_detected, rass$fdrs$num_detected[100:399])
abline(0,1)

test.iters <- parallel::mclapply(1:100, \(i) {
    suppressWarnings(
        association(nam.result = rnam, y = "case", verbose = FALSE)$fdrs$fdr
    )
}, mc.cores = 12)

test.mat <- do.call(cbind, test.iters)

all.equal(test.mat[,1], test.mat[,2])

cor.p <- apply(test.mat[100:399,], 2, \(x){ 
    data.frame(cor = cor(p$fdrs$fdr, x), 
               mae = mean(abs(p$fdrs$fdr - x)))
}) |> do.call(what = rbind)

plot(cor.p)


##################################################
# One more time with batch and covariates 
######
# library(reticulate)
# library(HSSRscripts)
# devtools::document(".")
# use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
# builtins <- import_builtins()
# # py_run_string("import vars")
# source_python("~/Git/cna-rcna-test/test.py")
# rcna_demo <- readRDS("inst/ext/demo_data.rds")

# # # # NAM with batch
rcna_demo$samplem$male <- as.numeric(rcna_demo$samplem$male)
rnam2 <- nam(data = rcna_demo, y = "case", 
             batches = "batch", 
             ridge.crit = "kurtosis",
             ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 0),
             verbose = F)

prn <- p2$namresid |> as.matrix() |> as("Matrix")
rrn <- Matrix::t(rnam2$resid.NAM.T)
all.equal(prn, rrn)
plot(prn[,100], rrn[,100])
abline(0,1)

plot(prn[1,], rrn[1,], pch = 20, col = "red")
abline(0,1)

# # # # ASSOCIATION with batch
rass2 <- association(rcna_demo, y = "case", batches = "batch", seed = 12,
                     ridge.crit = "kurtosis", use.logp = F,
                     ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 
                                1e-1, 1e-2, 1e-3, 1e-4, 0))

plot(rass2$fdrs$fdr[100:399], p2$fdrs$fdr)
abline(0,1)

plot(p2$fdrs$threshold, p2$fdrs$fdr)
points(rass2$fdrs$threshold[100:399], rass2$fdrs$fdr[100:399], pch = 20, col ="red")


# # # # NAM with batch and covs
rnam3 <- nam(data = rcna_demo, y = "case", 
             batches = "batch", covs = "male",
             ridge.crit = "kurtosis",
             ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 0),
             verbose = F)

prn <- p3$namresid |> as.matrix() |> as("Matrix")
rrn <- Matrix::t(rnam3$resid.NAM.T)
all.equal(prn, rrn)
plot(prn[,1], rrn[,1])
abline(0,1)


# # # # ASSOCIATION with batch and covs
rass3 <- association(rcna_demo, y = "case", 
                     batches = "batch", covs = "male",
                     seed = 12, use.logp = F,
                     ridge.crit = "kurtosis", 
                     ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 
                                1e-1, 1e-2, 1e-3, 1e-4, 0))

plot(p2$fdrs$threshold, p2$fdrs$fdr)
points(rass3$fdrs$threshold[100:399], rass3$fdrs$fdr[100:399], pch = 20, col ="red")


plot(p2$ncorrs, rass2$ncorrs[,1], pch = 20, col = "red")
abline(0,1)

#
### Alternative demo file?
#
library(SingleCellExperiment)
library(bluster)
library(igraph)
library(flotsam)
library(HSSRscripts)
library(rnndescent)
library(BiocNeighbors)
library(anndataR)
library(uwot)
devtools::document(".")

mysce <- qs2::qs_read("~/Data/Melanie/Donlin-MS-18649/Processed/allSCObject_Processed_DimRed_Clustered_Markers_Annotated_SubClustered_DE_SLI_WBC_10subclust.qs2")

foo <- uwot::similarity_graph(reducedDim(mysce, "HARMONY"), n_neighbors = 20)

colPair(mysce, "connectivities") <- foo

colData(mysce) <- colData(mysce) |> 
    as.data.frame() |> 
    dplyr::select(Sample, SLI, WBCSF, Batch, label, custom_annotation) |> 
    DataFrame()

rowData(mysce)$subsets <- NULL
rowData(mysce)$subsets <- NULL

rowData(mysce) <- rowData(mysce) |> 
    as.data.frame() |> 
    dplyr::select(ID:subsets_Ribo, scranHVGs) |> 
    DataFrame()


anndataR::write_h5ad(mysce, 
                     "~/Data/Melanie/Donlin-MS-18649/Processed/allSCObject_Processed_DimRed_Clustered_Markers_Annotated_SubClustered_DE_SLI_WBC_10subclust.h5ad")
saveRDS(mysce, "inst/ext/real_data.rds")



library(reticulate)
library(SingleCellExperiment)
devtools::document(".")
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test_realdata.py")
mysce <- readRDS("inst/ext/real_data.rds")
foo <- uwot::similarity_graph(reducedDim(mysce, "HARMONY"), n_neighbors = 20)
colPair(mysce, "connectivities") <- foo
rcna_demo <- createCNAListSCE(mysce, y = "SLI", sample.key = "Sample", 
                              sample.vars = c("SLI", "WBCSF", "Batch"), 
                              graph = "connectivities")



# diffusion check on real
diffuse_step1 <- d1 |> as.matrix() |> as("CsparseMatrix")

f <- as.formula(paste0("~ 0 + ", rcna_demo$samplem_key))
s <- Matrix::sparse.model.matrix(f, rcna_demo$obs)
colnames(s) <- gsub(rcna_demo$samplem_key, '\\1', colnames(s))
rownames(s) <- rcna_demo$obs[[rcna_demo$obs_key]]
diffuseStep <- function(a, s) {
    stopifnot((isa(a, "Matrix") || is.matrix(a)) && nrow(a) == ncol(a))
    stopifnot(isa(s, "Matrix") || is.matrix(s))
    stopifnot(nrow(s) == ncol(a))
    
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    rownames(res) <- rownames(s)
    return(res) 
}

rdiffuse_s1 <- s <- diffuseStep(rcna_demo$connectivities, s)
all.equal(rdiffuse_s1, diffuse_step1)
rdiffuse_s2 <- s <- diffuseStep(rcna_demo$connectivities, s)
diffuse_step2 <- d2 |> as.matrix() |> as("CsparseMatrix")
all.equal(rdiffuse_s2, diffuse_step2)
s <- diffuseStep(rcna_demo$connectivities, s) # skip 3
rdiffuse_s4 <- s <- diffuseStep(rcna_demo$connectivities, s)
diffuse_step4 <- d4 |> as.matrix() |> as("CsparseMatrix")
all.equal(rdiffuse_s4, diffuse_step4)
s <- diffuseStep(rcna_demo$connectivities, s) # skip 5-7
s <- diffuseStep(rcna_demo$connectivities, s)
s <- diffuseStep(rcna_demo$connectivities, s)
rdiffuse_s8 <- s <- diffuseStep(rcna_demo$connectivities, s)
diffuse_step8 <- d8 |> as.matrix() |> as("CsparseMatrix")
all.equal(rdiffuse_s8, diffuse_step8)

# ALL # [1] TRUE


#
### checking nam step...both took 4 steps
#
devtools::document(".")
pysnorm <- as.matrix(snorm) |> as("CsparseMatrix")
rsnorm <- buildNAM(rcna_demo, max.steps = 3, min.steps = 3, verbose = T)
all.equal(pysnorm, rsnorm)
# [1] TRUE


pyqcnam <- as.matrix(qc_nam) |> as("CsparseMatrix")
rqcnam <- qcNAM(rsnorm)$NAM
all.equal(pyqcnam, rqcnam)
# [1] TRUE

rnam <- buildNAM(rcna_demo, max.steps = 3, min.steps = 3, verbose = T)
resid_rnam <- residNAM(rnam)
rsvd_nam <- svdNAM(resid_rnam$NAM)
all.equal(abs(as.matrix(svd_nam[[1]])), abs(as.matrix(rsvd_nam[[1]])))
# [1] TRUE

# testing API
# # starting with nam
a <- nam(data = rcna_demo, y = "SLI", max.steps = 3, min.steps = 3, verbose = T, 
         ridge.crit = "kurtosis", ddof = 1, scale.method = "sd")
b <- association(rcna_demo, nam.result = a, y = "SLI", seed = 20)

# # starting with association
b <- association(rcna_demo, y = "SLI", seed = 20,
                 max.steps = 3, min.steps = 3, verbose = T, 
                 ridge.crit = "kurtosis", ddof = 1, scale.method = "sd")

b <- associationSCE(sce = mysce, y = "SLI", sample.key = "Sample", 
                    graph = "connectivities", seed = 20, ridge.crit = "kurtosis",
                    max.steps = 3, min.steps = 3, ddof = 1, scale.method = "sd")


all.equal(b$p, p$p)
all.equal(S4Vectors::metadata(b)$CNA$p, p$p)




a2 <- nam(data = rcna_demo, y = "SLI", batches = "Batch", partial.by = "batches",
          max.steps = 3, min.steps = 3, verbose = T, ridge.crit = "kurtosis", 
          ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 0),
          ddof = 1, scale.method = "sd")
a3 <- nam(data = rcna_demo, y = "SLI", batches = "Batch", covs = 'WBCSF', 
          partial.by = "batches", max.steps = 3, min.steps = 3, verbose = T,
          ridge.crit = "kurtosis", ddof = 1, scale.method = "sd",
          ridges = c(1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 0))

b2 <- association(rcna_demo, nam.result = a2, y = "SLI", seed = 1)
b3 <- association(rcna_demo, nam.result = a3, y = "SLI", seed = 11)







all.equal(b$p, p$p)
all.equal(b2$p, p2$p)

# N>B> we find a seed that generates a p-value that matches python implementation
pseeds <- parallel::mclapply(1:30, \(i){
    association(rcna_demo, nam.result = a3, y = "SLI", return.nam = FALSE, seed = i)$p
}, mc.cores = 10) |> unlist()

which.min(abs(pseeds - p$p))
which.min(abs(pseeds - p2$p))
which.min(abs(pseeds - p3$p))


plot(b$fdrs$threshold, b$fdrs$fdr, col = "dodgerblue", pch = 20, 
     xlim = c(0,1), ylim = c(0,1), main = "CNA Comparison",
     xlab = "Corr Bins", ylab = "FDRs")
points(p$fdrs$threshold, p$fdrs$fdr, col = "coral", pch = 20)
legend(x = 0, y = 1, pch = 20, col = c("dodgerblue", "coral"), legend = c("rcna", "pycna"))

plot(b2$fdrs$threshold, b2$fdrs$fdr, col = "dodgerblue", pch = 20, 
     xlim = c(0,1), ylim = c(0,1), main = "CNA Comparison with Batches")
points(p2$fdrs$threshold, p2$fdrs$fdr, col = "coral", pch = 20)
legend(x = 0, y = 1, pch = 20, col = c("dodgerblue", "coral"), legend = c("rcna", "pycna"))

plot(b3$fdrs$threshold, b3$fdrs$fdr, col = "dodgerblue", pch = 20, 
     xlim = c(0,1), ylim = c(0,11), main = "CNA Comparison with Batches and Covs")
points(p3$fdrs$threshold, p3$fdrs$fdr, col = "coral", pch = 20)
legend(x = 0, y = 10, pch = 20, col = c("dodgerblue", "coral"), legend = c("rcna", "pycna"))






# only attributes differ...seems to be a general difference in data order between
# SCE and annData?
pM <- as(as.matrix(p$M), "CsparseMatrix")
aM <- as(as(a$M, "CsparseMatrix"), "dsCMatrix")
all.equal(pM, aM)
# [1] TRUE



spPnamresid <- as(as.matrix(p$namresid), "sparseMatrix") |> summary()
spAnamresid <- as(Matrix::t(a$resid.NAM.T), "CsparseMatrix") |> summary()
all.equal(spPnamresid, spAnamresid)
# TRUE

spPnamresid <- as(as.matrix(p2$namresid), "sparseMatrix") |> summary()
spAnamresid <- as(Matrix::t(a2$resid.NAM.T), "CsparseMatrix") |> summary()
all.equal(spPnamresid, spAnamresid)
# "Component “x”: Mean relative difference: 7.270417e-05"

spPnamresid <- as(as.matrix(p3$namresid), "sparseMatrix") |> summary()
spAnamresid <- as(Matrix::t(a3$resid.NAM.T), "CsparseMatrix") |> summary()
all.equal(spPnamresid, spAnamresid)
# "Component “x”: Mean relative difference: 7.395043e-05"

# plot(spPnamresid$x, spAnamresid$x, pch = ".")
# abline(0, 1, col = "red")

all.equal(abs(a$NAM_sampleXpc), abs(as.matrix(p$namresid_sampleXpc)))
# TRUE

p$namresid_nbhdXpc <- as.matrix(p$namresid_nbhdXpc)

all.equal(abs(a$NAM_nbhdXpc), abs(p$namresid_nbhdXpc))
# [1] "Mean relative difference: 0.04627592"
plot(a$NAM_nbhdXpc[,16], p$namresid_nbhdXpc[,16], pch = 20)
abline(0,1, col = "red")

sapply(seq_len(ncol(a$NAM_nbhdXpc)), \(i){
    all.equal(abs(unname(a$NAM_nbhdXpc[,i])), abs(p$namresid_nbhdXpc[,i]))
})
# linalg.svd and r svd have found different orthoganal vectors for the final 
# component

all.equal(a$NAM_svs[1:length(p$namresid_svs)], as.vector(p$namresid_svs))
# TRUE ... length differs because svd in python is calculated on tcrossprod(NAM)
#          instead of the full NAM so only 
# TODO: check with yakir whether this is just leftover from October 19, 2023 NOTICE on repo

#
##
###
####
# ASSOCIATION TESTING
####
###
##
#
devtools::document(".")
a <- nam(data = rcna_demo, y = "SLI", max.steps = 3, min.steps = 3, verbose = T, 
         ridge.crit = "kurtosis", ddof = 1, scale.method = "sd")

# N>B> we find a seed that generates a p-value that matches python implementation
pseeds <- parallel::mclapply(1:100, \(i){
    association(rcna_demo, nam.result = a, y = "SLI", return.nam = FALSE, seed = i)$p
}, mc.cores = 10) |> unlist()

myseed <- which.min(abs(pseeds - p$p))
b <- association(rcna_demo, nam.result = a, y = "SLI", return.nam = TRUE, seed = myseed)

b$p
# p$p
# [1] 0.004995005

all.equal(as.vector(p$beta), as.vector(b$beta))
# TRUE

all.equal(b$fdr_5p_t, p$fdr_5p_t)
# [1] "Mean relative difference: 0.003484692"

all.equal(b$fdr_10p_t, p$fdr_10p_t)
# [1] "Mean relative difference: 0.002774763"

all.equal(p$k, b$k)
# TRUE

identical(all(p$kept), all(b$keptcells))
# TRUE

all.equal(as.vector(p$ks), b$ks)
# TRUE

all.equal(as.matrix(p$M), as.matrix(b$M))
# TRUE

all.equal(as.matrix(p$nam), as.matrix(Matrix::t(b$raw.NAM.T)))
# [1] TRUE

all.equal(as.matrix(p$namresid), as.matrix(Matrix::t(b$resid.NAM.T)))
# [1] TRUE

all.equal(abs(b$NAM_nbhdXpc), abs(as.matrix(p$namresid_nbhdXpc)))
sapply(1:ncol(b$NAM_nbhdXpc), \(i){
    all.equal(abs(b$NAM_nbhdXpc[,i]), abs(as.matrix(p$namresid_nbhdXpc)[,i]))
})
# all TRUE except last

all.equal(abs(b$NAM_sampleXpc), abs(as.matrix(p$namresid_sampleXpc)))
# TRUE

# `namresid_varexp` is missing from output
p$namresid_varexp

all.equal(as.vector(p$ncorrs), unname(b$ncorrs[,1]))
# [1] TRUE

bfdrs <- b$fdrs
pfdrs <- p$fdrs

plot(pfdrs$threshold, pfdrs$fdr, pch = 2, col = "dodgerblue", xlim = c(0,1), ylim = c(0,1)) # 
points(bfdrs$threshold, bfdrs$fdr, pch = 6, col = "purple4")

all.equal(as.vector(p$yresid_hat), as.vector(b$yhat))
# [1] TRUE

all.equal(as.vector(p$yresid), as.vector(b$ycond))
# [1] TRUE

# NOT EXPECTED TO BE TRUE   
all.equal(as.vector(p$nullminps), as.vector(b$minps_null))
all.equal(p$nullr2_mean, b$nullr2_mean)
all.equal(p$nullr2_std, b$nullr2_std)

library(ggplot2)
library(ggExtra)

pl <- data.frame(py = as.vector(p$nullminps), r = as.vector(b$minps_null)) |> 
    ggplot(aes(x = py, y = r)) + 
    geom_point() + 
    theme(legend.position = "none")
ggMarginal(pl, type = "density")    

hist(as.vector(pl$data[,1]))
hist(as.vector(pl$data[,2]), add = T, col = rgb(0,1,0,0.5))
# NOT EXPECTED TO BE TRUE  

# # # N.B. seed above seed selection comment!
all.equal(p$p, b$p)
# [1] TRUE

all.equal(p$r, b$r)
# [1] TRUE

all.equal(p$k, b$k)
# [1] TRUE

all.equal(p$r2, b$r2)
# [1] TRUE

all.equal(as.vector(p$r2_perpc), as.vector(b$r2_perpc))


plot(tail_counts(b$fdrs$threshold, b$ncorrs)[,1], 
     tail_counts2(b$fdrs$threshold, b$ncorrs)[1,1:300])
abline(0,1)
plot(tail_counts3(b$fdrs$threshold, b$ncorrs)[1,], 
     tail_counts2(b$fdrs$threshold, b$ncorrs)[1,])
abline(0,1)

z <- b$fdrs$threshold
znull <- b$ncorrs_null
n <- length(z)
znull2 <- znull^2
z2  <- z^2
ix  <- order(z2)   
iix <- order(ix)   
atol = 1e-8
rtol = 1e-5
bins <- c(z2[ix] - atol - rtol * z2[ix], Inf)

hist_mat <- apply(znull2, 2, \(zn2){
    tabulate(findInterval(zn2, bins), nbins = n)
})

tails <- apply(hist_mat, 2, \(col) { rev(cumsum(rev(col))) })
tails <- t(tails)[, iix, drop = FALSE]

tails2 <- apply(znull2, 2, \(zn2){
    col <- tabulate(findInterval(zn2, bins), nbins = n)
    rev(cumsum(rev(col)))
})

library(microbenchmark)
library(ggplot2)
mb <- microbenchmark(
    r_tc = {tail_counts(b$fdrs$threshold, b$ncorrs)},
    pytc = {tail_counts2(b$fdrs$threshold, b$ncorrs)}
)
autoplot(mb)





