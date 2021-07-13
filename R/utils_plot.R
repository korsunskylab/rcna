#' Scatterplot with color 
#' 
#' @param coords Coordinate table 
#' @param y Vector for colors 
#' @param label (optional) name for y vector 
#' @return ggplot object 
#' 
#' @export 
dimplot.generic <- function(coords, y, label=NULL) {
    plt <- cbind(coords, y) %>% 
        dplyr::sample_frac(1L, TRUE) %>% 
        data.frame() %>% 
        ggplot2::ggplot(ggplot2::aes_string(colnames(coords)[1], colnames(coords)[2], color = 'y')) + 
            ggplot2::geom_point(size = 0.5) + 
            ggthemes::scale_color_gradient2_tableau() + 
            ggplot2::theme_classic(base_size = 16)    
    if (!is.null(label)) 
        plt <- plt + ggplot2::labs(color = label)
    return(plt)
}

#' Scatterplot with ncorr and FDR  
#' 
#' @param coords Coordinate table 
#' @param res result from association function
#' @param fdr_thresh FDR threshold 
#' @return ggplot object 
#' 
#' @export 
#' @importFrom rlang .data
dimplot.ncorr <- function(coords, res, fdr_thresh=0.05) {
    .thresh <- subset(res$fdrs, fdr < fdr_thresh)$threshold[1]    
    plt_df <- cbind(coords, correlation = res$ncorrs) %>% 
        dplyr::mutate(passed = abs(correlation) > .thresh) %>% 
        dplyr::sample_frac(1L, TRUE) %>% 
        data.frame() 
    ggplot2::ggplot(plt_df, ggplot2::aes_string(colnames(coords)[1], colnames(coords)[2])) + 
        ggplot2::geom_point(data = subset(plt_df, !passed), size = .5, color = 'grey', alpha = .05) + 
        ggplot2::geom_point(data = subset(plt_df, passed), size = .5, ggplot2::aes(color = correlation), alpha = .5) + 
#         ggplot2::geom_point(data = subset(.data, passed), size = .5, color = 'grey', alpha = .05) + 
#         ggplot2::geom_point(data = subset(.data, passed), size = .5, aes(color = correlation), alpha = .5) + 
        ggthemes::scale_color_gradient2_tableau() + 
        ggplot2::theme_classic(base_size = 16)    
}

