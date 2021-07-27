library("fda")
library("robustbase")

tukey.loss <- function(x,k=4.685){
  n <- length(x)
  salida <- rep(0,n)
  for(i in 1:n){
    if(abs(x[i])<=k){salida[i] <- 1-(1-(x[i]/k)^2)^3 #(x[i])^6/(6*k^4)-(x[i])^4/(2*k^2)+(x[i])^2/2
    }else{
      salida[i] <- 1 #k^2/6
    }
  }
  return(salida)
}

my.norm.2 <- function(x){
  return( sqrt(sum(x^2)) )
}


#' Classical knot selection
# #' @importFrom splines bs
# #' @importFrom stats lm
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
select.nknots.cl <- function(y,Z,X,degree.spline=3){
  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))
  grilla.tes <- seq(0,1,length=n)

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
      base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]

    sal <- stats::lm(y~Z.aux+Xspline)
    betas <- as.vector(sal$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux #- mean(aux) No se necesita porque están centrados con la integral
    }

    regresion.hat <- stats::predict(sal) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    BIC[nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+q+1) #q+1 es la cantidad de lineales
  }
  posicion <- which.min(BIC)
  nknots <- posicion-1 #Dec?a "knots" en lugar de decir nknots... creo

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline
  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots, nbasis = nbasis, kj=kj)
  return(salida)
}



select.nknots.rob <- function(y, Z, X, degree.spline=3, method="MM", maxit=100){

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }


  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  RBIC <- rep(0,length(grid.nknots))
  grilla.tes <- seq(0,1,length=n)

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
      base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]



    #- Tukey MM estimator -#
    control <- lmrob.control(trace.level = 0,         # 0
                             nResample   =  500,      # 500 default
                             tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                             subsampling = 'simple',  #
                             rel.tol     = 1e-5,      # 1e-7
                             refine.tol  = 1e-5,      # 1e-7
                             k.max       = 2e3,       # 200
                             maxit.scale = maxit,       # 200 #2e3
                             max.it      = maxit)       # 50 #2e3
    sal.r  <- lmrob(y ~ Z.aux+Xspline, control = control)
    #sal.r <- MASS::rlm(y~Z.aux+Xspline, method=method, maxit=maxit)

    betas <- as.vector(sal.r$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]
    desvio.hat <- sal.r$scale #sal.r$s
    vv <- sum(Mpsi(sal.r$res / desvio.hat, cc  = control$tuning.psi,  psi = control$psi, deriv = -1))


    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux #- mean(aux) No se necesita porque están centrados con la integral
    }

    regresion.hat.r <- stats::predict(sal.r) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)

    tuk <- tukey.loss( (y - regresion.hat.r)/desvio.hat )
    RBIC[nknots+1] <- log( (desvio.hat^2)*sum(tuk) )+ (log(n)/(2*n))*(nbasis+q+1) #q+1 porque q de la parte lineal y 1 de la constante. O sea, q+1 es la cantidad de lineales.
  }
  posicion <- which.min(RBIC)
  nknots <- posicion-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, RBIC=RBIC, grid.nknots=grid.nknots, nbasis = nbasis, kj = kj, vv=vv)
  return(salida)

}

#' Classical Partial Linear Additive Model
# #' @importFrom splines bs
# #' @importFrom stats lm
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
plam.cl <- function(y, Z, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  if( is.null(nknots) ){
    AUX <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
    nbasis <- AUX$nbasis
    nknots <- AUX$nknots
    kj <- AUX$kj
  }else{
    nbasis <- d*(nknots + degree.spline)
    kj <- nknots + degree.spline
  }

  Mat.X <- as.list(rep(0,d))
  #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
  Xspline <- NULL

  grilla.tes <- seq(0,1,length=1000)
  for (ell in 1:d){
    nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
    base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                        norder = (degree.spline+1),
                                        breaks = nodos.spl)
    aux <- getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]
  #print(dim(Mat.X[[ell]]))

  sal <- stats::lm(y~Z.aux+Xspline)
  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    #correc[ell] <- mean(aux)
    #gs.hat[,ell] <- aux - mean(aux)
    gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- punto <- as.matrix(np.point)
      }else{
        prediccion <- punto <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- punto <- np.point
    }
    np <- dim(punto)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    #grilla.tes <- seq(0,1,length=n)
    for(ell in 1:d){
      nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
      base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- getbasismatrix(punto[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- getbasismatrix(grilla.tes, base.beta)
      spl.final.new <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final.new[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final.new[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }
    nMat.new <- dim(Mat.X.new[[ell]])[2]

    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat.new*(ell-1)+1):(nMat.new*ell)] %*% coef.spl[(nMat.new*(ell-1)+1):(nMat.new*ell)] )
        #prediccion[,ell] <- aux - correc[ell]
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat.new*(ell-1)+1):(nMat.new*ell)] %*% coef.spl[(nMat.new*(ell-1)+1):(nMat.new*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Robust Partial Linear Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
#' @export
plam.rob <- function(y, Z, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3, maxit=100, method="MM"){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  if( is.null(nknots) ){
    AUX <- select.nknots.rob(y, Z, X, degree.spline=degree.spline, method=method, maxit=maxit)
    nknots <- AUX$nknots
    nbasis <- AUX$nbasis
    kj <- AUX$kj
  }else{
    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline+1)
    kj <- (nknots + degree.spline) #(nknots + degree.spline + 1)
  }

  Mat.X <- as.list(rep(0,d))
  #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
  Xspline <- NULL
  grilla.tes <- seq(0,1,length=1000)
  for (ell in 1:d){
    nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
    base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                        norder = (degree.spline+1),
                                        breaks = nodos.spl)
    aux <- getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]
  control <- lmrob.control(trace.level = 0,         # 0
                           nResample   =  500,      # 500 default
                           tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                           subsampling = 'simple',  #
                           rel.tol     = 1e-5,      # 1e-7
                           refine.tol  = 1e-5,      # 1e-7
                           k.max       = 2e3,       # 200
                           maxit.scale = maxit,       # 200 #2e3
                           max.it      = maxit)       # 50 #2e3
  sal  <- lmrob(y ~ Z.aux+Xspline, control = control)
  #sal <- MASS::rlm(y~Z.aux+Xspline ,method=method,maxit=maxit)

  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]
  sigma.hat <- sal$s
  vv <- sum(Mpsi(sal$res / sigma.hat, cc  = control$tuning.psi,  psi = control$psi, deriv = -1))

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    #correc[ell] <- mean(aux)
    #gs.hat[,ell] <- aux - mean(aux)
    gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, alpha=alpha.hat, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj) #alpha.hat=alpha.hat+sum(correc)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- punto <- as.matrix(np.point)
      }else{
        prediccion <- punto <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- punto <- np.point
    }
    np <- dim(punto)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    #grilla.tes <- seq(0,1,length=n)
    for(ell in 1:d){

      nodos.spl   <- seq(min(X[,ell]), max(X[,ell]), length = (2+nknots))
      base.beta   <- create.bspline.basis(rangeval = c(0, 1),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- getbasismatrix(punto[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- getbasismatrix(grilla.tes, base.beta)
      spl.final.new <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final.new[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final.new[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }
    nMat.new <- dim(Mat.X.new[[ell]])[2]


    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #prediccion[,ell] <- aux - correc[ell]
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, alpha=alpha.hat, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion, vv=vv)
    return(salida)
  }
}

