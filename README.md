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

``` r
lim.cl <- lim.rob <- matrix(0, 2, 3)
x0 <- fit.full.cl$X
par(mfrow=c(1,3))
for(j in 1:3) {
  re <- fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-j])
  lim.cl[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5)
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full.cl$g.matrix[oo,j], lwd=5, col='red')
}
```

![](README_files/figure-markdown_github/plot%20cla-1.png)

``` r
par(mfrow=c(3,3))
for(j in 1:3) {
  re.cl <- fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-j])
  lim.cl[,j] <- c(min(re), max(re))
  re.rob <- fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-j])
  lim.rob[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5, ylim=lim.cl[,j])
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full.cl$g.matrix[oo,j], lwd=5, col='red')
  plot(re.rob ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5,ylim=lim.rob[,j])
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.rob$g.matrix[oo,j], lwd=5, col='blue')
  plot(x0[oo,j],fit.full.cl$g.matrix[oo,j], ylab="",type="l", col='red', ylim=c(min(lim.cl[1,j],lim.rob[1,j]),max(lim.cl[2,j],lim.rob[2,j])), lwd=5)
  lines(x0[oo,j],fit.rob$g.matrix[oo,j],col='blue',lwd=5)
}
```

![](README_files/figure-markdown_github/plot%20both-1.png)

``` r
DF <- as.data.frame(cbind(
  X,
  fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-1]),
  fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-2]),
  fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-3]),
  fit.full.cl$g.matrix,
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3]),
  fit.rob$g.matrix
))
names(DF) <- c("x1","x2","x3","re.cl.1", "re.cl.2", "re.cl.3", "cl.1", "cl.2", "cl.3", "re.rob.1", "re.rob.2", "re.rob.3", "rob.1", "rob.2", "rob.3")

library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.6.3

``` r
ggplot(DF,aes(x=x1,y=re.cl.1)) +
  geom_point(alpha=0.5,col='red')+
  geom_line(data=DF,aes(x=x1,y=cl.1),col='red',size=2)+
  geom_point(aes(x=x1,y=re.rob.1),col='blue',alpha=0.5)+
  geom_line(data=DF,aes(x=x1,y=rob.1),col='blue',size=2)+
  theme(
    axis.title.y=element_blank()
    )
```

![](README_files/figure-markdown_github/ggplot-1.png)

``` r
ggplot(DF,aes(x=x2,y=re.cl.2)) +
  geom_point(alpha=0.5,col='red')+
  geom_line(data=DF,aes(x=x2,y=cl.2),col='red',size=2)+
  geom_point(aes(x=x2,y=re.rob.2),col='blue',alpha=0.5)+
  geom_line(data=DF,aes(x=x2,y=rob.2),col='blue',size=2)+
  theme(
    axis.title.y=element_blank()
  )
```

![](README_files/figure-markdown_github/ggplot-2.png)

``` r
ggplot(DF,aes(x=x3,y=re.cl.3)) +
  geom_point(alpha=0.5,col='red')+
  geom_line(data=DF,aes(x=x3,y=cl.3),col='red',size=2)+
  geom_point(aes(x=x3,y=re.rob.3),col='blue',alpha=0.5)+
  geom_line(data=DF,aes(x=x3,y=rob.3),col='blue',size=2)+
  theme(
    axis.title.y=element_blank()
  )
```

![](README_files/figure-markdown_github/ggplot-3.png)
