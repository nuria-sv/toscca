
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toscca

<!-- badges: start -->
<!-- badges: end -->

The goal of toscca is to find sparse canonical vectors maximally correlating two data sets. The number of nonzero variables must be stated beforehand. The descriptive algorithm can be seen below. 

![alt text](https://github.com/nuria-sv/toscca/blob/master/algorithm.png)

## Installation

You can install the development version of toscca from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("nuria-sv/toscca")
```

This is a basic example which shows you how to reproduce the analysis
described in “A framework for interpretation and testing of sparse
canonical correlations”.

## Simulations

We simulate data for 3 underlying signals of different size on a
high-dimensional setting. That is, $$
\mathbf{X} \in \mathbb{R}^{n\times p} \quad \text{and} \quad \mathbf{Y}\in \mathbb{R}^{n\times q}, 
\\
\text{where } n = 100, \text{ } p = 2500 \text{ and } q = 500.
$$

### Canonical Correlation Analysis

We use the method described in the paper, Thresholded Ordered Sparse
Canonical Correlation Analysis (TOSCCA), to uncover the underlying
processes linking the data.

``` r
X = standardVar(X0)
Y = standardVar(Y0)
K = 4                                       # number of components to be estimated
nonz_x = rep(100, K)                        # number of nonzero variables for X
nonz_y = rep(100, K)                        # number of nonzero variables for Y
init   = "uniform"                          # type of initialisation
cca_toscca  = toscca::toscca(X, Y, nonz_x, nonz_y, K, init, combination = FALSE)
#> 
#> __________________________________________ 
#>  For component K =  1 : 
#> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
#> k-fold cv max. cancor 
#>             0.9999112 
#> 
#>  ........................................ 
#>  # nonzero A: 100
#>  # nonzero B: 100
#>  ........................................
```

<img src="man/figures/README-toscca-1.png" width="100%" />

    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 8 
    #> k-fold cv max. cancor 
    #>             0.9805795 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 100
    #>  # nonzero B: 100
    #>  ........................................

<img src="man/figures/README-toscca-2.png" width="100%" />

    #> 
    #> __________________________________________ 
    #>  For component K =  3 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 10 
    #> k-fold cv max. cancor 
    #>             0.9325936 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 100
    #>  # nonzero B: 100
    #>  ........................................

<img src="man/figures/README-toscca-3.png" width="100%" />

    #> 
    #> __________________________________________ 
    #>  For component K =  4 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0.00049 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.1114875 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 100
    #>  # nonzero B: 100
    #>  ........................................

<img src="man/figures/README-toscca-4.png" width="100%" />

``` r
cpev_toscca = sapply(1:K, function(k) cpev.fun(X, X%*%cca_toscca$alpha[,1:k]))

perm_toscca = perm.toscca(X, Y, nonz_x, nonz_y, K = K, init, draws = 100, cancor = cca_toscca$cancor)
```

<img src="man/figures/README-toscca-5.png" width="100%" />

    #> Empirical p-values:
    #> 0
    #> 0
    #> 0
    #> 0.33
    #> NULL

We repeat the analysis using the Penalised Matrix Analysis approach from
the PMA R-package. I have modified the main function to allow for
different initialisations (such as random uniform and random normal) to
compare to our models’ performance.

``` r
pma_lambda = PMA::CCA.permute(X, Y, typex = "standard", typez = "standard")
cca_pma = CCA_pma_random(X, Y, K = 4, typex = "standard", typez = "standard", start = "uniform", 
                         penaltyx = pma_lambda$bestpenaltyx, pma_lambda$bestpenaltyz)
#> 123456789101112131415
#> 123456789101112131415
#> 123456789101112131415
#> 123456789101112131415
```

<img src="man/figures/README-plots-1.png" width="100%" />

Then we run the algorithm several times with different sparsity levels
for $\alpha$, keeping those for $\beta$ fixed. We observed the variables
making up the signal are repeatedly selected with a higher frequency
than the noise variables. The implication being that noise variable are
selected once the signal has been retrieved.

    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9998522 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.9452523 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>              0.999942 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 60
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 10 
    #> k-fold cv max. cancor 
    #>             0.9625605 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 60
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9999381 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 110
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 17 
    #> k-fold cv max. cancor 
    #>             0.9560377 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 110
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9998549 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 160
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 10 
    #> k-fold cv max. cancor 
    #>             0.9407342 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 160
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9997678 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 210
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 11 
    #> k-fold cv max. cancor 
    #>             0.9341117 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 210
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9997949 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 260
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0.00155 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.1725572 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 260
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9997038 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 310
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0.00232 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.2124609 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 310
    #>  # nonzero B: 10
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  1 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 4 
    #> k-fold cv max. cancor 
    #>             0.9996146 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 360
    #>  # nonzero B: 50
    #>  ........................................ 
    #> 
    #> __________________________________________ 
    #>  For component K =  2 : 
    #> | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  | 100 %          Common convergence error: 0 & Iterations: 12 
    #> k-fold cv max. cancor 
    #>             0.9458648 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 360
    #>  # nonzero B: 10
    #>  ........................................

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Cummulative percentage of explained variance (CPEV) and correlation between components.

AS discussed in the paper, a common measure to evaluated whether or not
to compute another component is the CPEV, as seen in the equation below.
However, if the estimated canonical vectors are correlated, then the
subsequent components are limited in the amount of new information they
are bringing by adding one more of them.

$$
\frac{tr(\mathbf{\gamma}_k^T\mathbf{\gamma}_k)}{\mathbf{X}^T\mathbf{X}}
$$

We therefore scale the CPEV according to correlation between previous
components. We refer to this as adjusted CPEV.

$$
 ACPEV(\mathbf{\gamma}_k) = CPEV(\mathbf{\gamma}_k)\cdot\prod_{i<k}(1 - |cor(\mathbf{\gamma}_i, \mathbf{\gamma}_k)|), \quad \text{for } k = 2, \dots, K
$$

``` r
cpev_pma = sapply(1:K, function(k) cpev.fun(X, X%*%cca_pma$u[,1:k]))

auto_cor = data.frame(first_second = c(cor(cca_pma$u[,1] , cca_pma$u[,2]), 
                                 cor(cca_toscca$alpha[,1], cca_toscca$alpha[,2])), 
                      first_third  = c(cor(cca_pma$u[,1] , cca_pma$u[,3]), 
                                 cor(cca_toscca$alpha[,1], cca_toscca$alpha[,3])), 
                      first_fourth = c(cor(cca_pma$u[,1] , cca_pma$u[,4]), 
                                 cor(cca_toscca$alpha[,1], cca_toscca$alpha[,4])),
                      second_third = c(cor(cca_pma$u[,2] , cca_pma$u[,3]), 
                                 cor(cca_toscca$alpha[,2], cca_toscca$alpha[,3])), 
                      third_fourth = c(cor(cca_pma$u[,3] , cca_pma$u[,4]), 
                                 cor(cca_toscca$alpha[,3], cca_toscca$alpha[,4]))
)
rownames(auto_cor) <- c("PMA", "TOSCCA")

adj_cpev_toscca = c(cpev_toscca[1], 
                    sapply(2:K, function(k) cpev_toscca[k]*prod(1-abs(auto_cor[2,k-1:k]))))
adj_cpev_pma    = c(cpev_pma[1],
                    sapply(2:K, function(k) cpev_pma[k]*prod(1-abs(auto_cor[1,k-1:k]))))

plot(t(auto_cor[2,]), type = "b", pch = 19, col = "deepskyblue", xlab = "Correlation to first component", ylab = "Autocorrelation", xaxt = "n", yaxt="n", ylim = c(-1,1), lty = 1)
axis(1, at=1:length(colnames(auto_cor)), lab=colnames(auto_cor), las=TRUE)
axis(2, at=pretty(c(-1,1)), lab=pretty(c(-1,1)), las=TRUE)
grid(nx = NULL, ny = NULL, lty = 2, col = alpha("lightgray", 0.4),lwd = 2) 
points(t(auto_cor[1,]), type = "b", pch = 15, col = "tomato", lty = 2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r

plot((adj_cpev_toscca), type = "b", pch = 19, col = "deepskyblue", xlab = "Component", ylab = "Adjusted CPEV", xaxt = "n", yaxt="n", ylim = c(0,max(max(adj_cpev_pma), max(adj_cpev_toscca))), lty = 1)
axis(1, at=1:K, lab=1:K, las=TRUE)
axis(2, at=pretty(c(0,max(max(adj_cpev_pma), max(adj_cpev_toscca)))), lab=pretty(c(0,max(max(adj_cpev_pma), max(adj_cpev_toscca)))), las=TRUE)
grid(nx = NULL, ny = NULL, lty = 2, col = alpha("lightgray", 0.4),lwd = 2) 
points((adj_cpev_pma), type = "b", pch = 15, col = "tomato", lty = 2)
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

``` r


kable(
  auto_cor,
  col.names = colnames(auto_cor),
  # row.names = row.names(auto_cor),
  # digits = 2,
  caption = "Autocorrelation between components"
  )
```

<table>
<caption>
Autocorrelation between components
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
first_second
</th>
<th style="text-align:right;">
first_third
</th>
<th style="text-align:right;">
first_fourth
</th>
<th style="text-align:right;">
second_third
</th>
<th style="text-align:right;">
third_fourth
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PMA
</td>
<td style="text-align:right;">
-0.0101010
</td>
<td style="text-align:right;">
0.6278225
</td>
<td style="text-align:right;">
-0.0101010
</td>
<td style="text-align:right;">
0.1227418
</td>
<td style="text-align:right;">
-0.0101010
</td>
</tr>
<tr>
<td style="text-align:left;">
TOSCCA
</td>
<td style="text-align:right;">
-0.0217817
</td>
<td style="text-align:right;">
-0.0096114
</td>
<td style="text-align:right;">
-0.0045191
</td>
<td style="text-align:right;">
-0.0040061
</td>
<td style="text-align:right;">
0.0097445
</td>
</tr>
</tbody>
</table>

``` r

adj_cpev = rbind(adj_cpev_pma, adj_cpev_toscca)
rownames(adj_cpev) <- c("PMA", "TOSCCA")
colnames(adj_cpev) <- paste0("K", 1:K)
kable(
  adj_cpev,
  col.names = colnames(adj_cpev_pma),
  # row.names = row.names(auto_cor),
  # digits = 2,
  caption = "Autocorrelation between components"
  )
```

<table>
<caption>
Autocorrelation between components
</caption>
<tbody>
<tr>
<td style="text-align:left;">
PMA
</td>
<td style="text-align:right;">
0.0020653
</td>
<td style="text-align:right;">
0.0081771
</td>
<td style="text-align:right;">
0.0068477
</td>
<td style="text-align:right;">
0.0120499
</td>
</tr>
<tr>
<td style="text-align:left;">
TOSCCA
</td>
<td style="text-align:right;">
0.0072966
</td>
<td style="text-align:right;">
0.0082782
</td>
<td style="text-align:right;">
0.0090524
</td>
<td style="text-align:right;">
0.0095376
</td>
</tr>
</tbody>
</table>
