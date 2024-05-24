run_ddc <- function(X, Y = NULL) {
    ddc_input <- X
    if (!is.null(Y)) {
        ddc_input <- cbind(ddc_input, Y)
    }
    ddc_res <- cellWise::DDC(ddc_input, DDCpars = list(fastDDC = TRUE, silent = TRUE))
    ddc_imp <- ddc_res$Ximp

    p <- ncol(X)
    X_imp <- ddc_imp[, 1:p]
    colnames(X_imp) <- colnames(X)
    Y_imp <- NA
    if (!is.null(Y)) {
        Y_imp <- as.matrix(ddc_imp[, (p + 1)], ncol = 1)
        colnames(Y_imp) <- names(Y)
    }

    return(list(X = X_imp, Y = Y_imp))
}
