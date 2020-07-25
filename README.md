The following real data example of the implementation of robust estimators based on B-splines under partially linear additive model is part of a work in progress done in collaboration with Prof. Dr. Graciela Boente.

Let first install package <code>rplam</code>.

``` r
library(devtools)
install_github("alemermartinez/rplam")
library(rplam)
```

We will use the Air Quality data set available in <code>R</code>.

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

Classical estimation of a partially linear additive model

``` r
fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)
```

and the robust proposal

``` r
fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)
```

When no number of internal knots is specified, functions select.knots.cl or select.knots.rob, respectively, is used.

The estimations obtained for each linear coefficient for both classical and robust approaches are shown it the following Table.

|                 |        col2|         col3|
|-----------------|-----------:|------------:|
| *λ*             |   49.230699|   42.6283885|
| *μ*             |  -13.088702|  -12.8228650|
| *π*             |   -4.768714|   -0.2128605|
| *β*             |   -2.936557|    1.5242527|
| *β*<sub>1</sub> |  -16.630411|   -9.2577636|

The following three plots correspond to the classical (in red and dashed line) and robust (in blue and solid line) fits for the additive functions with their respectively partial residuals. ![](README_files/figure-markdown_github/ggplot-1.png)![](README_files/figure-markdown_github/ggplot-2.png)![](README_files/figure-markdown_github/ggplot-3.png)

It seems curves obtained for each additive component are quite differents. For this reason, we are going to study the residuals obtained by the robust estimator.

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

which correspond to observations 23, 34, 53 and 77. Highlighting the partial residuals of these four observations (in pink) we obtain the following plots: ![](README_files/figure-markdown_github/ggplot%20highlighted-1.png)![](README_files/figure-markdown_github/ggplot%20highlighted-2.png)![](README_files/figure-markdown_github/ggplot%20highlighted-3.png)

Now, we remove these four observations from the original data set and re-calculate the classical estimator

``` r
out.pos <- c(23,34,53,77)
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- Z[-out.pos]

fit.del.cl <- plam.cl(y.del, Z.del, X.del)
```

and plot the new curves obtained with the classical fit (in red dashed line) using the data without the potential outliers identified by the robust fit together with the curves obtained by the robust fit on the original data set. ![](README_files/figure-markdown_github/ggplot%20final-1.png)![](README_files/figure-markdown_github/ggplot%20final-2.png)![](README_files/figure-markdown_github/ggplot%20final-3.png)
