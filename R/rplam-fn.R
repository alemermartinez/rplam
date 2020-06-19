library(MASS)
library(splines)

tukey.loss <- function(x,k=4.685){
  n <- length(x)
  salida <- rep(0,n)
  for(i in 1:n){
    if(abs(x[i])<=k){salida[i] <- (x[i])^6/(6*k^4)-(x[i])^4/(2*k^2)+(x[i])^2/2
    }else{salida[i] <- k^2/6}
  }
  return(salida)
}

#Knot selection

select.nknots.cl <- function(y=y,covariate=X,factor=Z,degree.spline=3){
  n <- length(y)
  d <- dim(X)[2]
  q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5

  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]

    Z.f <- as.factor(Z)
    lev.Z <- levels(Z.f)
    dummies <- matrix(0,n,nlevels(Z.f)-1)
    for(k in 1:(nlevels(Z.f)-1)){
      dummies[,k] <- as.numeric(Z.f == lev.Z[k+1])
    }


    #--- Classical estimator ---#

    #sal <- lm(y~Z.f+Mat.Slp.T1+Mat.Slp.T2+Mat.Slp.T3)
    sal <- lm(y~dummies+Xspline)
    betas <- as.vector(sal$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux - mean(aux)
    }

    regresion.hat <- predict(sal) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- nknots + degree.spline + 1
    BIC[nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+d)
  }
  posicion <- which.min(BIC)
  nknots <- posicion-1 #Decía "knots" en lugar de decir nknots... creo

  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots)
  return(salida)
}



select.nknots.rob <- function(y=y,covariate=X,factor=Z,degree.spline=3,seed=26){
  n <- length(y)
  d <- dim(X)[2]
  q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5

  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]

    Z.f <- as.factor(Z)
    lev.Z <- levels(Z.f)
    dummies <- matrix(0,n,nlevels(Z.f)-1)
    for(k in 1:(nlevels(Z.f)-1)){
      dummies[,k] <- as.numeric(Z.f == lev.Z[k+1])
    }

    #- Tukey MM estimator -#
    set.seed(seed)
    sal.r <- rlm(y~Z.f+Xspline ,method="MM",maxit=100)

    betas <- as.vector(sal.r$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux - mean(aux)
    }

    regresion.hat.r <- predict(sal.r) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- nknots + degree.spline + 1
    desvio.hat <- sal.r$s
    tuk <- tukey.loss( (y - regresion.hat.r)/desvio.hat )
    RBIC[nknots+1] <- log( (desvio.hat^2)*sum(tuk) )+ (log(n)/(2*n))*(nbasis+d)
  }
  posicion <- which.min(RBIC)
  nknots <- posicion-1

  salida <- list(nknots=nknots, RBIC=BIC, grid.nknots=grid.nknots)
  return(salida)
}

