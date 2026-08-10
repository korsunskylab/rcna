#' Scatterplot with color 
#' 
#' Generic 2D plot with coloring. Takes a set of coordinates and y, where y is used
#' to color each row in coords. 
#' 
#' @param coords A matrix or dataframe with at least 2 columns containing the 
#' x and y coordinates to plot.
#' @param y A vector of values to use for coloring each point in coords.
#' @param label (optional) A character string to name the y vector.\cr
#' Default: NULL
#' @param reorder A character string specifying how the points should be ordered.
#' One of "shuffle" or "signif" where "shuffle" randomly shuffles the points to 
#' produce an even distribution of values and "signif" puts the largest values of
#' y on top.
#' @param threshold An (optional) numerical value below which values of y are 
#' not colored. Note that abs(y) > threshold is checked rather than y > threshold. \cr
#' Default: NULL
#' @param color.low,color.high Low and high color for correlation values.\cr
#' Default: "blue4" and "red4", respectively
#' @param point.size,point.alpha Numeric values controlling point size and alpha.\cr
#' Defaults: 0.5 (both)
#' 
#' @returns A ggplot object 
#' 
#' @export
dimplotCoords <- function(coords, y, 
                          label = NULL, 
                          reorder = c("shuffle", "signif"),
                          threshold = NULL,
                          color.low = "blue3",
                          color.high = "red3",
                          point.size = 0.5,
                          point.alpha = 0.5) {
    stopifnot(ncol(coords) >= 2)
    stopifnot(nrow(coords) == length(y))
    reorder <- match.arg(reorder, c("shuffle", "signif"))
    
    coords <- coords[,1:2]
    colnames(coords) <- c("Dim1", "Dim2")
    plt_df <- cbind.data.frame(coords, y = y) 
    
    if(!is.null(threshold)){
        plt_df <- dplyr::mutate(plt_df, passed = abs(y) > threshold)
    } else {
        plt_df <- dplyr::mutate(plt_df, passed = TRUE)
    }
    
    if(reorder == "shuffle"){
        plt_df <- dplyr::slice_sample(plt_df, prop = 1, replace = FALSE) 
    } else {
        plt_df <- dplyr::arrange(plt_df, abs(y))
    }
    
    plt <- ggplot2::ggplot(plt_df, ggplot2::aes(x = Dim1, y = Dim2)) + 
        ggplot2::geom_point(data = subset(plt_df, !passed), 
                            size = point.size, 
                            color = 'gray60', 
                            alpha = point.alpha) + 
        ggplot2::geom_point(data = subset(plt_df, passed), 
                            ggplot2::aes(color = y),
                            size = point.size, 
                            alpha = point.alpha) + 
        ggplot2::scale_color_gradient2(low = color.low, 
                                       mid = 'gray60',
                                       high = color.high) + 
        ggplot2::theme_classic(base_size = 16)
    
    if (!is.null(label)) 
        plt <- plt + ggplot2::labs(color = label)
    return(plt)
}

#' Scatter plot of correlation values based on FDR threshold 
#' 
#' @param sce An sce object to which has the results of [association] via 
#' [associationSCE]. Mutually exclusive with `res`.
#' @param res A result output from [association] function. Mutually exclusive 
#' with `sce`.
#' @param dim.use A character string or NULL. \cr
#' Default: "CNA"
#' @param label (optional) A character string to name the y vector, and by 
#' extension, the legend.\cr
#' Default: "corr"
#' @param fdr.thresh A numeric FDR threshold (\eqn{\alpha}), values larger than 
#' this are colored gray.
#' @param reorder A characterstring specifying whether the data should be reordered
#' to place significant values on top ("signif") or randomly shuffle the data ("shuffle").\cr
#' Default: "signif"
#' @param color.low,color.high Low and high color for correlation values.\cr
#' Default: "blue3" and "red3", respectively
#' @param point.size,point.alpha Numeric values controlling point size and alpha.\cr
#' Defaults: 0.5 (both)
#' 
#' @returns A ggplot object 
#' 
#' @export 
dimplotCorr <- function(sce, res, 
                        dim.use = "CNA",
                        label = "corr", 
                        fdr.thresh = 0.05,
                        reorder = "signif",
                        color.low = "blue3",
                        color.high = "red3",
                        point.size = 0.5,
                        point.alpha = 0.5) {
    
    if(missing(sce) && missing(res)){
        stop("Missing required `sce` or `res` inputs")
    }
    if(!missing(sce) && !missing(res)){
        stop("`sce` and `res` are mutually exclusive")
    } 

    if(!missing(sce)){
        if(!dim.use %in% SingleCellExperiment::reducedDimNames(sce)){
            stop(dim.use, " is missing from `reducedDim(sce)`.")
        }
        coords <- SingleCellExperiment::reducedDim(sce, "CNA")
        y <- sce$cna_ncorrs
        label <- "correlation"
        thresh <- S4Vectors::metadata(sce)$CNA$fdrs |> 
            subset(fdr < fdr.thresh) |> 
            dplyr::slice_head(n = 1) |> 
            dplyr::pull(threshold) 
    } else {
        coords <- res$NAM_nbhdXpc
        y <- res$ncorrs
        label <- "correlation"
        thresh <- subset(res$fdrs, fdr < fdr.thresh) |> 
            dplyr::slice_head(n = 1) |> 
            dplyr::pull(threshold) 
    }  
    if(length(thresh) == 0){
        warning("Specified threshold not observed in the results. Likely no FDR ",
                "reached the specified threshold. Setting it to 1 for this plot.", 
                immediate. = TRUE)
        thresh <- 1
    }
    
    plt <- dimplotCoords(coords = coords, 
                         y = y, 
                         label = label, 
                         reorder = reorder,
                         threshold = thresh,
                         color.low = color.low,
                         color.high = color.high, 
                         point.size = point.size, 
                         point.alpha = point.alpha)
    return(plt)
}

