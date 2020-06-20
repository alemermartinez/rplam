Let first install package <code>rplam</code>.

``` r
library(devtools)
install_github("alemermartinez/rplam")
library(rplam)
```

The following is a real-data example of the implemention of the robust estimator based on B-splines under a partial linear additive model.

``` r
attach(airquality)
aux <- !is.na(Solar.R)
y <- Ozone[aux]
x1 <- Temp[aux]
x2 <- Wind[aux]
x3 <- Solar.R[aux]
z <- Month[aux]

na.y <- !is.na(y)
y <- y[na.y]
x1 <- x1[na.y]
x2 <- x2[na.y]
x3 <- x3[na.y]
z <- z[na.y]

X <- cbind(x1,x2,x3)
Z <- as.matrix(z)
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
sal <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
sal$nknots
```

    ## [1] 9

And so is the number of internal knots selected using the RBIC criteria used for the robust proposal.

``` r
sal.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
sal.rob$nknots
```

    ## [1] 9
