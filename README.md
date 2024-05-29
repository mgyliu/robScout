# robscout

Implements RobScout, a robust covariance-regularized regression method.

Most of the core code is copied from https://github.com/cran/scout, which is the package that implements the Scout method for covariance-regularized regression.

## Installation

```r
devtools::install_github("mgyliu/robscout", dependencies = TRUE, build_vignettes = TRUE)
vignette("simulated_example")
```

# Example

Suppose you have training data `X_train` and `Y_train`. 

To run robScout(1,1) with DDC as the preprocessing step:
```r
library(robscout)

K <- 5 # Number of folds to use in cross-validation
scout_1something_stepwise(
    X_train, Y_train,
    K = K, p2 = 1, nlambda1 = 50, nlambda2 = 50,
    ddc_first = TRUE, ddc_with_response = TRUE,
    cov_method = "default"
)
```

To run robScout(2,1) with DDC as the preprocessing step:
```r
K <- 5 # Number of folds to use in cross-validation
cv.scout_alternating_lasso(
    X_train, Y_train,
    K = K, p1 = 2, nlambda1 = 50, nlambda2 = 50,
    ddc_first = TRUE, ddc_with_response = TRUE,
    cov_method = "default"
)
```

# References

Witten, Daniela M., and Robert Tibshirani. "Covariance-regularized regression and classification for high dimensional problems." Journal of the Royal Statistical Society Series B: Statistical Methodology 71.3 (2009): 615-636.
