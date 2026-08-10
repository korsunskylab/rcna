#' Calculate kurtosis of the rows or columns of a sparseMatrix 
#' 
#' @param x A sparseMatrix for which the row or column kurtosis should be calculated
#' 
#' @param margin If not missing, then an integer specifying the margin over which 
#' to calculate kurtosis 1 for rows and 2 for columns. If missing, the kurtosis of 
#' the whole matrix is calculated. 
#' 
#' @param na.rm A logical controlling whether NAs should be removed from the calculation.\cr
#' Default: FALSE
#' 
#' @return the row- or column-wise kurtosis of x
#' 
#' @keywords internal
Kurtosis <- function(x, margin, na.rm = FALSE){
    if(!missing(margin)){
        stopifnot(margin %in% c(1L, 2L) && length(margin) == 1)
    }
    stopifnot(is.matrix(x) || isa(x, "Matrix"))
    if(missing(margin)){
        mu <- Matrix::mean(x, na.rm = na.rm)
        m2 <- Matrix::mean(x^2, na.rm = na.rm)
        m3 <- Matrix::mean(x^3, na.rm = na.rm)
        m4 <- Matrix::mean(x^4, na.rm = na.rm)
    } else if(margin == 1){
        mu <- Matrix::rowMeans(x, na.rm = na.rm)
        m2 <- Matrix::rowMeans(x^2, na.rm = na.rm)
        m3 <- Matrix::rowMeans(x^3, na.rm = na.rm)
        m4 <- Matrix::rowMeans(x^4, na.rm = na.rm)
    } else if(margin == 2){
        mu <- Matrix::colMeans(x, na.rm = na.rm)
        m2 <- Matrix::colMeans(x^2, na.rm = na.rm)
        m3 <- Matrix::colMeans(x^3, na.rm = na.rm)
        m4 <- Matrix::colMeans(x^4, na.rm = na.rm)
    } else {
        stop("'margin' only supports none or a single dimension of a matrix (1 or 2).")
    }
    
    # Expanded 4th central moment: E[(X-mu)^4]
    central4 <- m4 - 4*mu*m3 + 6*mu^2*m2 - 3*mu^4
    sigma4 <- (m2 - mu^2)^2
    
    kurt <- central4 / sigma4
    return(kurt)
}

#' A prop.table alternative for sparseMatrices
#' 
#' @param x A sparseMatrix
#' @param margin A vector giving the margins to split by. E.g., for a matrix 1 
#' indicates rows, 2 indicates columns. 
#' 
#' @note Unlike [base::prop.table()], joint normalization (`margin = c(1, 2)`)
#' is not supported. Rows or columns that sum to zero will produce `NaN`.
#' 
#' @returns A proportion table for x over the given margin in sparseMatrix format.
#' 
#' @keywords internal
propTable <- function(x, margin){
    stopifnot(margin %in% c(1L, 2L) && length(margin) == 1)
    stopifnot(is.matrix(x) || isa(x, "Matrix"))
    if(margin == 1){
        res <- Matrix::t(Matrix::crossprod(x, Matrix::Diagonal(x = 1 /  Matrix::rowSums(x))))
    } else if(margin == 2){
        res <- x %*% Matrix::Diagonal(x = 1 /  Matrix::colSums(x)) 
    } else {
        stop("Only margins 1 or 2 are supported.")
    }
    dimnames(res) <- dimnames(x)
    return(res)
}

#' A scale function for sparseMatrix 
#' 
#' @param x A sparseMatrix
#' @param center A boolean, should x be centered
#' Default: TRUE
#' @param scale A boolean, should x be scaled
#' Default: TRUE
#' @param ddof An integer, either 0 or 1 specifying the denominator degrees of
#' freedom. Use 1L for the n - 1 denominator and 0L for n denominator.
#' Default: 1L
#' @param method A character string specifying the method for calculating the 
#' scaling vector. "rms" uses root-mean-square and "sd" uses standard deviation.
#' If centering then these are equivalent.
#' 
#' @note NAs are omitted as in [base::scale()]
#' 
#' @returns A dense sparseMatrix ("dgeMatrix")
#'
#' @keywords internal
Scale <- function(x, 
                  center = TRUE, 
                  scale = TRUE, 
                  ddof = 1L, 
                  method = c("rms", "sd")){
    method <- match.arg(method, c("rms", "sd"))
    stopifnot(isa(x, "Matrix") || is.matrix(x))
    stopifnot(any(ddof == c(0L, 1L)))
    if(center){
        centers <- Matrix::colMeans(x, na.rm = TRUE)
        x <- Matrix::t(Matrix::t(x) - centers)
    }
    if(scale){
        sden <- nrow(x) - ddof
        if(method == "rms"){
            scales <- sqrt(Matrix::colSums(x^2, na.rm = TRUE) / sden)
        } else {
            xbar <- Matrix::colMeans(x, na.rm = TRUE)
            scales <- sqrt(Matrix::rowSums((Matrix::t(x) - xbar)^2) / sden)
        }
        x <- Matrix::t(Matrix::t(x) / scales)
    }
    return(x)
}