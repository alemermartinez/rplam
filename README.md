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

The following three plots correspond to the classical (in red) and robust (in blue) fits for the additive functions with their respectively partial residuals. ![](README_files/figure-markdown_github/ggplot-1.png)![](README_files/figure-markdown_github/ggplot-2.png)![](README_files/figure-markdown_github/ggplot-3.png)

Since fits are not quite similar, we are going to study the residuals obtained by the robust estimator.

``` r
res.rob <- y-fit.rob$prediction
summary(res.rob)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -21.7551  -6.2581  -0.2832   2.6631   7.6366 103.6170

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

![](README_files/figure-markdown_github/residuals%20plots-1.png)

The residuals detected as outliers are:

``` r
p1.info$outliers[[1]]
```

    ## [1]  54.69827  70.38513  52.77027 103.61697

which correspond to observations 23, 34, 53 and 77. Highlighting the partial residuals of these four observations (in pink) we obtain the following plots:

    ## [1] 111  10

    ## 'data.frame':    111 obs. of  10 variables:
    ##  $ Temp   : num  190 118 149 313 299 99 19 256 290 274 ...
    ##  $ Wind   : num  7.4 8 12.6 11.5 8.6 13.8 20.1 9.7 9.2 10.9 ...
    ##  $ Solar.R: num  67 72 74 62 65 59 61 69 66 68 ...
    ##  $ re.1   : num  7.299 11.926 0.891 -3.972 -8.195 ...
    ##  $ re.2   : num  9.29 3.39 -18.36 -9.08 -10.12 ...
    ##  $ re.3   : num  -4.77 -8.03 -23.48 -21.53 -17.88 ...
    ##  $ g1     : num  0.5761 4.9682 6.051 1.0017 0.0892 ...
    ##  $ g2     : num  2.56 -3.57 -13.2 -4.1 -1.83 ...
    ##  $ g3     : num  -11.5 -15 -18.3 -16.6 -9.6 ...
    ##  $ Fit    : Factor w/ 1 level "Robust": 1 1 1 1 1 1 1 1 1 1 ...

    ## [1] 4 6

![](README_files/figure-markdown_github/ggplot%20highlighted-1.png)![](README_files/figure-markdown_github/ggplot%20highlighted-2.png)![](README_files/figure-markdown_github/ggplot%20highlighted-3.png)

Now, we remove these four observations from the original data set and re-calculate the classical estimator.
