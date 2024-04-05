scout1something <- function(x, y, p2, lam1s, lam2s, rescale, cov_method) {
    if (ncol(x) > 500) print("You are running scout with p1=1 and ncol(x) > 500. This will be slow. You may want to re-start and use p1=2, which is much faster.")
    if (min(lam1s) == 0 && min(lam2s) == 0 && ncol(x) >= nrow(x)) stop("don't run w/lam1=0 and lam2=0 when p>=n")
    if (sum(order(lam2s) == (1:length(lam2s))) != length(lam2s)) {
        stop("Error!!!! lam2s must be ordered!!!")
    }
    betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))
    if (sum(lam1s > 0 & lam1s < 1e-4) > 0 && ncol(x) >= nrow(x)) {
        warning("Non-zero lam1s that were smaller than 1e-4 were increased to 1e-4 to avoid problems with graphical lasso that can occur when p>n.")
        lam1s[lam1s < 1e-4 & lam1s != 0] <- 1e-4
    }
    for (i in 1:length(lam1s)) {
        g.out <- NULL
        if (i == 1 || is.null(g.out$w) || is.null(g.out$wi)) {
            if (lam1s[i] != 0) g.out <- glasso::glasso(est_cov(x, method = cov_method), rho = lam1s[i])
            if (lam1s[i] == 0) g.out <- list(w = est_cov(x, method = cov_method), wi = NULL)
        } else if (i != 1 && !is.null(g.out$w) && !is.null(g.out$wi)) {
            g.out <- glasso::glasso(est_cov(x, method = cov_method), rho = lam1s[i], start = "warm", w.init = g.out$w, wi.init = g.out$wi)
        }
        for (j in 1:length(lam2s)) {
            if (p2 == 0 || lam2s[j] == 0) {
                beta <- g.out$wi %*% est_cov(x, y, method = cov_method)
            } else if (p2 == 1 && lam2s[j] != 0) {
                if (j == 1) beta <- scout::crossProdLasso(g.out$w, est_cov(x, y, method = cov_method), rho = lam2s[j])$beta
                if (j != 1) {
                    if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
                        beta <- scout::crossProdLasso(g.out$w, est_cov(x, y, method = cov_method), rho = lam2s[j], beta.init = beta)$beta
                        # if got zero for smaller value of lambda 2,
                        # then no need to keep computing!!!
                    }
                }
            }
            if (rescale && sum(abs(beta)) != 0) beta <- beta * lsfit(x %*% beta, y, intercept = FALSE)$coef
            betamat[i, j, ] <- beta
        }
    }
    return(betamat)
}
