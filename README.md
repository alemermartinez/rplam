Let first install package <code>rplam</code>.

``` r
library(devtools)
install_github("alemermartinez/rplam")
library(rplam)
```

The following is a real-data example of the implemention of the robust estimator based on B-splines under a partial linear additive model.

``` r
data(airquality)
x <- airquality
x <- x[ complete.cases(x), ]
x <- x[, c('Ozone', 'Solar.R', 'Wind', 'Temp','Month')]
y <- as.vector(x$Ozone)
X <- as.matrix(x[, c('Solar.R', 'Wind', 'Temp')])
Z <- as.matrix(x[, 'Month'])
```

As the is taken as a categorical variable, we convert it into a factor variable.

``` r
Z <- as.factor(Z)
```

The degree selected for the spline basis is 3 for the three additive functions.

``` r
degree.spline <- 3
```

The number of internal knots selected by the BIC criteria for the classical estimator is 9.

``` r
nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots
```

    ## [1] 9

And so is the number of internal knots selected using the RBIC criteria used for the robust proposal.

``` r
nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots
```

    ## [1] 9

Classical estimation of a partial linear additive model

``` r
fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)
```

and the robust proposal

``` r
fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)
```

When no number of internal knots is specified, select.knots.cl or select.knots.rob, respectively, is used.

The following three plots correspond to the classical (in red) and robust (in blue) fits for the additive functions with their respectively partial residuals.

    ## [1] 222  10

    ## 'data.frame':    222 obs. of  10 variables:
    ##  $ x1  : num  190 118 149 313 299 99 19 256 290 274 ...
    ##  $ x2  : num  7.4 8 12.6 11.5 8.6 13.8 20.1 9.7 9.2 10.9 ...
    ##  $ x3  : num  67 72 74 62 65 59 61 69 66 68 ...
    ##  $ re.1: num  1.462 0.827 -2.498 -5.139 -11.741 ...
    ##  $ re.2: num  12.15 -3.86 -21.69 -13.82 -11.35 ...
    ##  $ re.3: num  -8.2 -14.4 -24.1 -22.6 -23.2 ...
    ##  $ g1  : num  -5.358 2.9541 3.0449 0.0282 -1.696 ...
    ##  $ g2  : num  5.33 -1.74 -16.15 -8.66 -1.31 ...
    ##  $ g3  : num  -15 -12.3 -18.6 -17.4 -13.2 ...
    ##  $ Fit : Factor w/ 2 levels "Classical","Robust": 1 1 1 1 1 1 1 1 1 1 ...

![](README_files/figure-markdown_github/ggplot-1.png)![](README_files/figure-markdown_github/ggplot-2.png)![](README_files/figure-markdown_github/ggplot-3.png)
