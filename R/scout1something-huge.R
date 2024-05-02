# Same as scout1something.R but using huge::huge.glasso instead
scout1something.huge <- function(x, y, p2, lam1s, lam2s, rescale, cov_method, scaleFun) {
    # if (ncol(x) > 500) {
    #     print("You are running scout with p1=1 and ncol(x) > 500. This will be slow. You may want to re-start and use p1=2, which is much faster.")
    # }

    if (min(lam1s) == 0 && min(lam2s) == 0 && ncol(x) >= nrow(x)) {
        stop("don't run w/lam1=0 and lam2=0 when p>=n")
    }

    betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))

    if (sum(lam1s > 0 & lam1s < 1e-4) > 0 && ncol(x) >= nrow(x)) {
        warning("Non-zero lam1s that were smaller than 1e-4 were increased to 1e-4 to avoid problems with graphical lasso that can occur when p>n.")
        lam1s[lam1s < 1e-4 & lam1s != 0] <- 1e-4
    }

    # If p2==0, no cov.output is needed
    cov.output <- !(p2 == 0)
    cov.est.xx <- est_cov(x, method = cov_method)
    cov.est.xy <- est_cov(x, y, method = cov_method)

    g.out.all <- huge::huge(cov.est.xx, method = "glasso", lambda = lam1s, cov.output = cov.output, verbose = FALSE)

    for (i in 1:length(lam1s)) {
        for (j in 1:length(lam2s)) {
            beta <- if (p2 == 0) {
                g.out.all$icov[[i]] %*% cov.est.xy
            } else {
                scout::crossProdLasso(g.out.all$cov[[i]], cov.est.xy, rho = lam2s[j])$beta
            }
            # Scaling coefficient, scout step 4
            # For simple linear regression w/o intercept, the coefficient of regression
            # of X*beta onto Y is the same as:
            #   cor(x%*%beta, y) / sd(x%*%beta) * sd(y)
            c.beta <- if (rescale & sum(abs(beta)) != 0) {
                as.numeric(est_cov(x %*% beta, y, method = cov_method, correlation = TRUE) * scaleFun(y) / scaleFun(x %*% beta))
            } else {
                1
            }

            betamat[i, j, ] <- beta * c.beta
        }
    }

    # return(list(betamat = betamat, c.betamat = c.betamat, x.betamat = x.betamat, y.mat = y.mat))
    return(betamat)
}
