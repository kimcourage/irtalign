# irtalign

<!-- badges: start -->
[![R-CMD-check](https://github.com/kimcourage/irtalign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kimcourage/irtalign/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Alignment optimization for multiple-group item response models.

irtalign estimates latent means and variances across many groups when item parameters are only approximately invariant, using the alignment method of Asparouhov and Muthén (2014). Dichotomous (2PL, 3PL) and polytomous (graded response, generalized partial credit) items can be mixed in one analysis. Each group is fitted with [mirt](https://cran.r-project.org/package=mirt), and the results include standard errors, the noninvariant parameters of each group, invariance R-squared, simplicity function contributions and factor mean comparisons.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kimcourage/irtalign")
```

## Example

Six countries answer four dichotomous and four four-category items. Item 2 is easier in country C3 and item 6 discriminates less in country C5.

```r
library(mirt)
library(irtalign)

set.seed(2026)
G <- 6
mu <- c(0, 0.3, -0.4, 0.6, 0.1, -0.2)
sigma <- sqrt(c(1, 1.3, 0.8, 1.1, 0.9, 1.2))
a <- matrix(c(1.2, 0.9, 1.5, 1.1, 1.4, 0.8, 1.0, 1.3))
d <- rbind(c(0.4, NA, NA), c(-0.6, NA, NA), c(1.0, NA, NA), c(0.2, NA, NA),
           c(1.5, 0.2, -1.0), c(2.0, 0.5, -0.8), c(1.0, -0.2, -1.5), c(1.8, 0.0, -1.2))
itemtype <- rep(c("2PL", "graded"), each = 4)
dat <- do.call(rbind, lapply(seq_len(G), function(g) {
  dg <- d
  if (g == 3) dg[2, 1] <- dg[2, 1] + 0.8   # item 2 is easier in country 3
  ag <- a
  if (g == 5) ag[6] <- 0.4                 # item 6 discriminates less in country 5
  resp <- simdata(ag, dg, N = 1500, itemtype = itemtype,
                  Theta = matrix(rnorm(1500, mu[g], sigma[g])))
  data.frame(country = paste0("C", g), resp)
}))

fit <- align_irt(dat, group = "country", nstarts = 30, seed = 1)
fit
```

```text
Alignment: FIXED, reference group C5 (automatic)
Groups: 6   Items: 8   Observations: 9,000
Simplicity function: SQRT, eps = 1e-04; metric: reference group variance 1
Simplicity function value: 98.2049 (reached by 15 of 31 starts)

 Group    N   Mean    SE Variance   SE 
    C1 1500 -0.103 0.108    1.002 0.254
    C2 1500  0.194 0.090    1.399 0.312
    C3 1500 -0.501 0.126    0.970 0.197
    C4 1500  0.575 0.133    1.290 0.369
    C5 1500  0.000 0.000    1.000 0.000
    C6 1500 -0.296 0.130    1.255 0.303
```

The reference group is the country whose mean is closest to zero in a free alignment. Means are in the metric of that group (mean 0, variance 1).

```r
summary(fit)
```

```text
Alignment: FIXED, reference group C5 (automatic)
Groups: 6   Items: 8   Observations: 9,000
Simplicity function: SQRT, eps = 1e-04; metric: reference group variance 1
Simplicity function value: 98.2049 (reached by 15 of 31 starts)

Approximate measurement (non)invariance for groups
 Slopes
   Item_1          C1 C2 C3 C4 C5 C6
   Item_2          C1 C2 C3 C4 C5 C6
   Item_3          C1 C2 C3 C4 C5 C6
   Item_4          C1 C2 C3 C4 C5 C6
   Item_5          C1 C2 C3 C4 C5 C6
   Item_6          C1 C2 C3 C4 (C5) C6
   Item_7          C1 C2 C3 C4 C5 C6
   Item_8          C1 C2 C3 C4 C5 C6
 Intercepts
   Item_1 c1       C1 C2 C3 C4 C5 C6
   Item_2 c1       C1 C2 (C3) C4 C5 C6
   Item_3 c1       C1 C2 C3 C4 C5 C6
   Item_4 c1       C1 C2 C3 C4 C5 C6
   Item_5 c1       C1 C2 C3 C4 C5 C6
   Item_5 c2       C1 C2 C3 C4 C5 C6
   Item_5 c3       C1 C2 C3 C4 C5 C6
   Item_6 c1       C1 C2 C3 C4 C5 C6
   Item_6 c2       C1 C2 C3 C4 C5 C6
   Item_6 c3       C1 C2 C3 C4 C5 C6
   Item_7 c1       C1 C2 C3 C4 C5 C6
   Item_7 c2       C1 C2 C3 C4 C5 C6
   Item_7 c3       C1 C2 C3 C4 C5 C6
   Item_8 c1       C1 C2 C3 C4 C5 C6
   Item_8 c2       C1 C2 C3 C4 C5 C6
   Item_8 c3       C1 C2 C3 C4 C5 C6

Noninvariant parameter-group combinations: 1.4%

Group means and variances
 Group    N   Mean    SE Variance   SE 
    C1 1500 -0.103 0.108    1.002 0.254
    C2 1500  0.194 0.090    1.399 0.312
    C3 1500 -0.501 0.126    0.970 0.197
    C4 1500  0.575 0.133    1.290 0.369
    C5 1500  0.000 0.000    1.000 0.000
    C6 1500 -0.296 0.130    1.255 0.303

Factor mean comparison at the 5% level in descending order
 Rank  Group      Mean   Groups with significantly smaller mean
    1  C4        0.575   C2 C5 C1 C6 C3
    2  C2        0.194   C5 C1 C6 C3
    3  C5        0.000   C6 C3
    4  C1       -0.103   C6 C3
    5  C6       -0.296   C3
    6  C3       -0.501   

Item parameters
   Item Parameter Noninvariant Invariant mean    R2 Contribution Variance
 Item_1         a          0/6          1.155 0.272        3.730    0.003
 Item_1        c1          0/6          0.497 0.986        3.467    0.003
 Item_2         a          0/6          0.961 0.752        4.143    0.004
 Item_2        c1          1/6         -0.546 0.000        7.117    0.110
 Item_3         a          0/6          1.434 0.585        5.998    0.022
 Item_3        c1          0/6          1.103 0.972        4.398    0.008
 Item_4         a          0/6          1.053 0.302        4.001    0.004
 Item_4        c1          0/6          0.328 0.986        3.344    0.002
 Item_5         a          0/6          1.330 0.791        4.388    0.006
 Item_5        c1          0/6          1.602 0.989        2.945    0.001
 Item_5        c2          0/6          0.294 0.989        3.316    0.002
 Item_5        c3          0/6         -0.929 0.975        3.980    0.007
 Item_6         a          1/6          0.756 0.259        6.063    0.038
 Item_6        c1          0/6          2.069 0.864        4.972    0.014
 Item_6        c2          0/6          0.542 0.976        3.504    0.003
 Item_6        c3          0/6         -0.743 0.967        2.873    0.001
 Item_7         a          0/6          0.967 0.719        2.930    0.001
 Item_7        c1          0/6          1.092 0.969        3.721    0.004
 Item_7        c2          0/6         -0.107 0.978        3.645    0.003
 Item_7        c3          0/6         -1.402 0.986        3.398    0.002
 Item_8         a          0/6          1.236 0.261        5.153    0.010
 Item_8        c1          0/6          1.893 0.970        4.229    0.006
 Item_8        c2          0/6          0.122 0.985        3.612    0.003
 Item_8        c3          0/6         -1.072 0.985        3.279    0.002
```

Groups in parentheses are noninvariant for that parameter. `Invariant mean` is the sample-size weighted mean of the aligned parameter over the invariant groups, `R2` is the share of the configural variation explained by the group means and variances, and `Contribution` is the parameter's share of the simplicity function.

```r
head(invariance_table(fit))
```

```text
    item param group configural  aligned        se invariant
1 Item_1     a    C1   1.157482 1.156143 0.1402028      TRUE
2 Item_1     a    C2   1.304257 1.102801 0.1455405      TRUE
3 Item_1     a    C3   1.139656 1.157123 0.1384519      TRUE
4 Item_1     a    C4   1.312302 1.155498 0.1448087      TRUE
5 Item_1     a    C5   1.261818 1.261818 0.1113588      TRUE
6 Item_1     a    C6   1.229224 1.097274 0.1510466      TRUE
```

## Functions

| Function | Output |
|---|---|
| `align_irt()` | Fits each group with `mirt::mirt()` and aligns the estimates |
| `align_fits()` | Aligns a list of single-group mirt fits |
| `align_pars()` | Aligns configural slopes and intercepts supplied as matrices |
| `summary()` | Noninvariance by parameter and group, means and variances, factor mean comparison and item parameter summary |
| `coef()` | Group means and variances with standard errors |
| `invariance_table()` | Configural and aligned parameters, standard errors and invariance status by group |
| `mean_comparison()` | Groups ranked by mean with the groups that have a significantly smaller mean |

## Model

Each item has one slope `a` and intercepts `c1, ..., cK` on the logit scale: `a * theta + c` for 2PL and 3PL items, the cumulative logits `a * theta + ck` for graded items, and the step intercepts `dk - d(k-1)` for generalized partial credit items. With configural estimates $a_{g0}, c_{g0}$ (latent mean 0 and variance 1 in every group), a group mean $\alpha_g$ and variance $\psi_g$ give the aligned parameters

$$a_g = a_{g0} / \sqrt{\psi_g}, \qquad c_g = c_{g0} - \alpha_g a_{g0} / \sqrt{\psi_g}.$$

The means and variances minimize

$$F = \sum_{\text{parameters}} \sum_{g < h} w_g w_h \, (\Delta_{gh}^2 + \epsilon)^{1/4}, \qquad w_g = \sqrt{N_g} / \overline{\sqrt{N}},$$

where $\Delta_{gh}$ is the difference between the aligned values of groups g and h, with $\epsilon = 10^{-4}$. `align = "fixed"` fixes the mean of the reference group at zero and `align = "free"` estimates all means. Standard errors follow from implicit differentiation of the optimality conditions and the delta method. The invariant set of each parameter starts from the largest set of groups connected by pairwise tests with p > 0.01 and is updated by comparing each group with the weighted mean of the set (p > 0.001).

## Supported models

Unidimensional single-group fits from `mirt::mirt()` with itemtypes `2PL`, `3PL`, `graded` and `gpcm`, and configural estimates of any model with one slope and one or more intercepts per item through `align_pars()`. The lower asymptote of 3PL items is reported in `$guess`.

## References

- Asparouhov, T., & Muthén, B. (2014). Multiple-group factor analysis alignment. *Structural Equation Modeling, 21*(4), 495–508.
- Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. *Journal of Statistical Software, 48*(6), 1–29.
- Muthén, B., & Asparouhov, T. (2014). IRT studies of many groups: The alignment method. *Frontiers in Psychology, 5*, 978.

## License

MIT
