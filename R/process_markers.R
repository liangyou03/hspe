#' Determines number of markers \code{n_markers} and marker list \code{mrkrs}.
#' @inheritParams hspe
process_markers <- function(Y, pure_samples, n_markers, markers, marker_method) {
    
    K <- length(pure_samples)
    
    if (is.null(markers)) {
        markers <- find_markers(Y = Y, pure_samples = pure_samples, marker_method = marker_method)
        if (is.null(n_markers)) {
            n_markers <- sapply(floor(0.1 * lengths(markers$L)), min, ncol(Y)/K)
        }
    }
    markers <- get_marker_list(markers)
    
    if (is.null(n_markers)) {
        n_markers <- lengths(markers)
    } else {
        if (length(n_markers) == 1) 
            n_markers <- rep(n_markers, K)
        
        wq_markers <- which(n_markers < 1)
        
        n_markers[wq_markers] <- floor(n_markers[wq_markers] * lengths(markers)[wq_markers])
        
    }
    
    n_markers <- sapply(n_markers, max, 1)
    
    mrkrs <- lapply(1:K, function(i) {
        markers[[i]][1:n_markers[i]]
    })
    names(mrkrs) <- names(pure_samples)
    
    return(list(n_markers = n_markers, mrkrs = mrkrs))
}
