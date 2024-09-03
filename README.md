
<!-- README.md is generated from README.Rmd. Please edit that file -->

\#fquantdr

<!-- badges: start -->
<!-- [![Build Status](https://travis-ci.com/elianachristou/fquantdr.svg?branch=master)](https://travis-ci.com/elianachristou/fquantdr) -->
<!-- badges: end -->

fquantdr: Dimension Reduction Techniques for Conditional Quantiles for
Functional Data
====================================================================

The R package `fquantdr` performs dimension reduction techniques for
conditional quantiles of a scalar response given the functional
predictors. Specifically, the method aims at replacing the
infinite-dimensional functional predictors with a few finite predictors
without losing important information on the conditional quantiles while
maintaining a flexible nonparametric model. For details of the
methodology, see Christou et al. (2024+).
<!-- [Christou, E., Solea, E., Wang, S., and Song, J. (2024+) Sufficient Dimension Reduction for Conditional Quantiles for Functional Data. *Journal*, volume, pages](link) -->

The main function of the package is `fcqs`, which estimates the
directions of the functional central quantile subspace. However, the
package includes more functions that are helpful to run `fcqs`.
Specifically, `mfsir` performs functional sliced inverse regression
(FSIR) of [Ferré and Yao
(2003)](https://doi.org/10.1080/0233188031000112845) and `sonf` performs
scalar-on-function linear regression. Moreover, `fundata` generates
functional data and `mcorr` computes the multiple correlation between
two matrices.

## Installation

You can install the released version of `fquantdr` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fquantdr")
```

and then use it during any `R` session by issuing the command

``` r
library(fquantdr)
```

The development version from [GitHub](https://github.com/) can be
installed using:

``` r
install.packages("devtools")
devtools::install_github("elianachristou/fquantdr")
```

## Introducing `fcqs`

This function computes the directions that span the $\tau$th functional
central quantile subspace. These directions represent functions that can
be linearly applied via the inner product to replace the
infinite-dimensional functional predictors with a few finite predictors
without losing important information on the conditional quantiles while
maintaining a flexible nonparametric model.

The function requires several inputs:

- `x` a 3-dimensional array ($n \times nt \times p$) representing the
  functional predictors, where $n$ is the number of observations, $nt$
  is the number of time points, and $p$ is the number of predictors.
- `y` a numeric vector of length $n$ representing the scalar response.
- `time` a numeric vector of length $nt$ representing the time points at
  which the functional data is evaluated.
- `nbasis` the number of basis functions for smoothing the functional
  predictors.
- `tau` the quantile level
- `dtau` the number of directions the user wants to extract. If not
  provided, the function will return $p$ directions.

The function then returns:

- `betacoef` the functional parameters that span the functional central
  quantile subspace
- `betax` the resulting sufficient predictors, calculated as the inner
  product between `betacoef` and `x`.

#### Example 1

This is a basic example that shows how to apply the function. First,
let’s define the basic parameters, such as the sample size, the number
of predictors, the number of times points, the vector of time points,
the quantile level, and the number of basis functions that we wish to
smooth the functional predictors.

``` r
library(fquantdr)

# define the parameters
set.seed(1234)
n <- 100
p <- 5
nt <- 101
time <- seq(0, 1, length = nt)
tau <- 0.5
nbasis <- 4
```

Then, we need to generate the functional predictors and the scalar
response. For the functional predictors, we can use the `fundata`
function that is available in the package. An additional argument,
`eta`, is required for the function and it represents the matrix of
coefficients with dimensions $n$ by $p \times nbasis$. For the scalar
response, we use the model
$$Y = 3 \langle \beta_1, X \rangle + \epsilon,$$ where $\epsilon$ is the
error term and $\beta_1$ represent an eigenfunction. In this case, the
$\tau$-th functional central quantile subspace is generated by
$\{\beta_1\}$.

``` r
# Generate the functional predictors
library(mvtnorm)
eta <- rmvnorm(n, mean = rep(0, p * nbasis))
data.output <- fundata(n, p, nbasis, time, eta)
xc <- data.output$xc

# Generate the scalar response
P <- eigen(cov(eta))$vectors
mfpca.scores <- eta %*% P
error <- rnorm(n)
y <- 3 * mfpca.scores[, 1] + error
```

Before moving on and for illustration purposes, we plot the first
functional predictor.

``` r
matplot(time, t(xc[, , 1]), type = "l", lty = 1, col = 1:n, xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

The purpose of the `fcqs` function is to estimate $\beta_1$ and form the
new predictor $\langle \beta_1, X \rangle$. We now run the function and
specify $d_\tau = 1$, so we can obtain the first direction.

``` r
result <- fcqs(xc, y, time, nbasis, tau, dtau = 1)
```

The first sufficient predictor $\langle \widehat{\beta}_1, X \rangle$ is
given by

``` r
result$betax
```

Since the true relationship is
$Y = 3 \langle \beta_1, X \rangle + \epsilon$, then a plot of $Y$
against $\langle \widehat{\beta}_1, X \rangle$ should display a linear
relationship.

``` r
plot(result$betax, y, xlab = 'First Sufficient Predictor')
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Another way to evaluate the performance of the methodology is to
calculate the correlation between the true $\langle \beta_1, X \rangle$
and the estimated $\langle \widehat{\beta}_1, X \rangle$. For that, we
can use the `mcorr` function of the package, which calculates the
multiple correlation between two matrices. The output is a number
between 0 and $d$, the dimension of the input. For this case, since we
have a one-dimensional sufficient predictor, a number closer to 1
indicates better performance.

``` r
true.pred <- mfpca.scores[, 1]
est.pred <- result$betax
mcorr(true.pred, est.pred)
#> [1] 0.9478022
```

#### Example 2

Let’s consider another example, where the $\tau$-th functional central
quantile subspace is two-dimensional. Let
$Y = \arctan(\pi \langle \beta_1, X \rangle) + 0.5 \sin(\pi \langle \beta_2, X \rangle / 6) + 0.1 \epsilon$.
Then,

``` r
y <- mfpca.scores[, 1]^3 + exp(mfpca.scores[, 2]) + error
result2 <- fcqs(xc, y, time, nbasis, tau, dtau = 2)

true.pred2 <- mfpca.scores[, 1:2]
est.pred2 <- result2$betax
mcorr(true.pred2, est.pred2)
#> [1] 1.664628
```

## Applications

Potential applications of fquantdr include, but are not limited to, the
following fields:

- Biomedical data analysis
- Financial data modeling
- Environmental data analysis

For more detailed examples, please refer to the package documentation.

## Functions

**mfsir: Functional Sliced Inverse Regression**

Performs FSIR that can incorporate multivariate functional predictors.

**Arguments:**

- X: (n x t x p) array of functional elements (centered).
- y: Response variable.
- H: Number of slices.
- nbasis: Number of basis functions.

**Returns:**

- phi: Eigenvectors.
- betas: Sufficient predictors.
- eigvalues: Eigenvalues.
- xfd.coef: Coefficients of the functional data.
- gx: Gram matrix of the basis functions.

### fcqs: Functional Central Quantile Subspace

Fits a scalar-on-function regression for multivariate functional data.

**Arguments:**

- x: (n x t x p) array of functional elements.
- y: Response variable.
- t: Vector of time points.
- tau: Quantile level
- d_tau: Dimension of the FCQS.
- nbasis: Number of basis functions.

**Returns:**

- ffun: Functional parameters that span the FCQS..
- betax: Resulting sufficient predictor.

## Additional Utility Functions

These utility functions support various operations related to functional
data analysis:

- **discretize**: Discretizes the response variable y into yunit slices,
  adding a small perturbation for stability.
- **slprob**: Computes slice proportions for the response variable.
- **slav**: Computes slice averages for the predictor coefficients.
- **symmetry**: Symmetrizes a matrix.

For more detailed examples and documentation, please refer to the
package documentation and vignettes.
