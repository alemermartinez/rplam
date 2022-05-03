#Ejemplo de selección de variables#
source('rplam-fn.R')
mu <- 5
beta <- c(1,2,3,0,0,0)
q <- length(beta)
eta1 <- function(x){
  return( 2*sin(pi*x)-4/pi )
}
eta2 <- function(x){
  return( exp(x)-(exp(1)-1) )
}
p <- 4
n <- 200
set.seed(23)
Z <- matrix(runif(n*q),n,q)
X <- matrix(runif(n*p),n,p)

sigma <- 0.2
eps <- rnorm(n,0,sd=sigma)

y <- mu+Z%*%as.matrix(beta)+eta1(X[,1])+eta2(X[,2])+eps

#Primero estimo sin selección de variables
#-- De forma clásica
sal.cl <- plam.cl(y,Z,X)
sal.cl$coef.const
sal.cl$coef.lin
dim(sal.cl$g.matrix)
for(i in 1:p){
  p.res <- y-sal.cl$coef.const-Z%*%as.matrix(sal.cl$coef.lin)-rowSums(sal.cl$g.matrix[,-i])
  plot(X[,i],p.res)
  ord <- order(X[,i])
  lines(X[ord,i],sal.cl$g.matrix[ord,i],col=2,lwd=2)
}

#-- De forma robusta
sal.rob <- plam.rob(y,Z,X)
sal.rob$coef.const
sal.rob$coef.lin
dim(sal.rob$g.matrix)
for(i in 1:p){
  p.res <- y-sal.rob$coef.const-Z%*%as.matrix(sal.rob$coef.lin)-rowSums(sal.rob$g.matrix[,-i])
  plot(X[,i],p.res)
  ord <- order(X[,i])
  lines(X[ord,i],sal.rob$g.matrix[ord,i],col='blue',lwd=2)
}

#Ahora con selección de variables
#- Con k y lambdas fijos
lambdas1 <- rep(0.1,q)
lambdas2 <- rep(0.25,p)
nknots <- 0
degree.spline <- 3
maxit <- 100
MAXITER <- 100
bound.control <- 10^(-3)
sal.rob <- plam.rob.vs.nknots.lambdas(y, Z, X, lambdas1=lambdas1, lambdas2=lambdas2, nknots=nknots)
sal.rob$is.zero
sal.rob$coef.const
sal.rob$coef.lin
for(i in 1:p){
  p.res <- y-sal.rob$coef.const-Z%*%as.matrix(sal.rob$coef.lin)-rowSums(sal.rob$g.matrix[,-i])
  plot(X[,i],p.res)
  ord <- order(X[,i])
  lines(X[ord,i],sal.rob$g.matrix[ord,i],col='blue',lwd=2)
}



