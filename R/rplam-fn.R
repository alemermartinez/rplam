

#' Tukey Loss Function
#' @export
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

#' Tukey Loss Function
#' @export
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

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
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
      gs.hat[,ell] <- aux - mean(aux)
    }

    regresion.hat <- stats::predict(sal) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    BIC[nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+q+1) #q+1 es la cantidad de lineales
  }
  posicion <- which.min(BIC)
  nknots <- posicion-1 #Decía "knots" en lugar de decir nknots... creo

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline
  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots, nbasis = nbasis, kj=kj)
  return(salida)
}


#' Robust knot selection
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
# #' @importFrom MASS rlm
#' @export
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

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]
    #dim(Xspline)[2]/4



    #- Tukey MM estimator -#
    sal.r <- MASS::rlm(y~Z.aux+Xspline, method=method, maxit=maxit)

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

    regresion.hat.r <- stats::predict(sal.r) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    desvio.hat <- sal.r$s
    tuk <- tukey.loss( (y - regresion.hat.r)/desvio.hat )
    RBIC[nknots+1] <- log( (desvio.hat^2)*sum(tuk) )+ (log(n)/(2*n))*(nbasis+q+1) #q+1 porque q de la parte lineal y 1 de la constante. O sea, q+1 es la cantidad de lineales.
  }
  posicion <- which.min(RBIC)
  nknots <- posicion-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, RBIC=RBIC, grid.nknots=grid.nknots, nbasis = nbasis, kj = kj)
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
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
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
  correc <- rep(0,d)
  for(ell in 1:d){
    aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    correc[ell] <- mean(aux)
    gs.hat[,ell] <- aux - mean(aux)
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
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
    np <- dim(np.point)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X.new[[ell]] <- splines::bs( np.point[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }


    for(k in 1:np){
      for(ell in 1:d){
        aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        prediccion[,ell] <- aux - correc[ell]
      }
    }
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Robust Partial Linear Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
#' @export
plam.rob <- function(y, Z, X, nknots=NULL, knots=NULL, degree.spline=3, maxit=100, method="MM"){
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
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]

  sal <- MASS::rlm(y~Z.aux+Xspline ,method=method,maxit=maxit)
  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]
  sigma.hat <- sal$s

  gs.hat <- matrix(0,n,d)
  correc <- rep(0,d)
  for(ell in 1:d){
    aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    correc[ell] <- mean(aux)
    gs.hat[,ell] <- aux - mean(aux)
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl
  salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
  return(salida)
}

########################
#- Variable selection -#
########################

#' SCAD penalty function
#' @examples
#' x <- seq(-5, 5, length=100)
#' scad.p(x)
#' @export
scad.p <- function(x, lambda, a=3.7){
  x <- as.vector(x)
  nt <- length(x)
  sal <- rep(0,nt)
  for(i in 1:nt){
    if(abs(x[i])<= lambda){
      sal[i] <- lambda*abs(x[i])
    }
    if( (lambda<abs(x[i])) & (abs(x[i])<=a*lambda) ){
    sal[i] <- -(x^2-2*a*lambda*abs(x[i])+lambda^2)/(2*(a-1))
    }else{
      sal[i] <- (a+1)*lambda^2/2
    }
  }
  return(sal)
}

#' Derivative of the SCAD penality function
#' @export
scad.d <- function(t, lambda, a=3.7){
  t <- as.vector(t)
  t <- abs(t) #Esta condición no la ponen pero no está bien si no.
  nt <- length(t)
  sal <- rep(0,nt)
  for(i in 1:nt){
    if(t[i]<=lambda){
      sal[i] <- lambda
    }else{
      aux <- a*lambda-t[i]
      if(aux>0){
        sal[i] <-  aux/(a-1)
      }else{
        sal[i] <- 0
      }
    }
  }
  return(sal)
}

#' Variable selection in robust PLAM with fixed lambdas
#' @export
plam.rob.vs.lambdas <- function(y, Z, X, lambdas1, lambdas2, nknots=NULL, knots=NULL, degree.spline=3, maxit=100, method="MM", MAXITER=100, bound.control=10^(-3)){
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
    for (ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]
    sal <- MASS::rlm(y~Z.aux+Xspline ,method=method, maxit=maxit)
    xdesign <- sal$x
    sigma.hat <- sal$s
    #Construyo el beta 0
    beta.ini <- as.vector(sal$coefficients)
    #beta.hat <- beta.ini[-1]
    #coef.lin <- beta.ini[2:(q+1)]
    #coef.spl <- beta.ini[(q+2):(1+q+nMat*d)]
    nbetas <- length(beta.ini)

    normgammaj <- rep(0,d)
    corte <- 1
    iter <- 0
    while( (corte>bound.control) & (iter<MAXITER)){
      iter <- iter +1
      #print(iter)
      regresion.hat <- xdesign%*%beta.ini
      res <- y-regresion.hat
      an <- quantile(abs(res),2*(n^(-1/2)))
      W <- diag( as.vector(tukey.loss((res)/sigma.hat)/(an+(res)^2) ))

      Sigmalambda <- matrix(0,nbetas,nbetas)
      for(i in 1:(q+1)){
        Sigmalambda[i,i] <- scad.d(beta.ini[i],lambda=lambdas1[i])/abs(beta.ini[i])
      }
      for(i in 1:d){
        Hj <- Hj.matrix(X[,i], nknots, degree.spline)
        gammaj <- as.matrix(beta.ini[(q+2+nMat*(i-1)):(nMat*i+(q+1))])
        normgammaj[i] <- sqrt( t(gammaj)%*%Hj%*%gammaj )
        Sigmalambda[(q+2+nMat*(i-1)):(nMat*i+(q+1)),(q+2+nMat*(i-1)):(nMat*i+(q+1))] <- as.numeric(scad.d(normgammaj[i],lambda=lambdas2[i])/(normgammaj[i]))*Hj
      }

      #Sigmalambda <- diag( scad.d(beta.ini,lambda=lambda)/abs(beta.ini))

      #options(show.error.messages = TRUE)

      try.sal <- try(
        AUX <- solve(t(xdesign)%*%W%*%xdesign + 1/2*n*Sigmalambda)
        #aa <- solve(t(xdesign)%*%W%*%xdesign + 1/2*n*Sigmalambda)
        #print(det(aa))
        #dimaa <- dim(aa)[1]
        #AUX <- .C("inverse", as.double(aa), salida=as.double(0), as.integer(dimaa))$salida
      )

      if(class(try.sal)!= 'try-error'){
        beta1 <- as.vector(AUX%*%t(xdesign)%*%W%*%y)
        corte <- my.norm.2(beta.ini-beta1)/my.norm.2(beta.ini)
        beta.ini <- beta1
      }else{
        beta1 <- beta.ini
        iter <- MAXITER
      }
    }

    beta.hat <- beta1[-1]
    coef.lin <- beta1[2:(q+1)]
    coef.spl <- beta1[(q+2):(1+q+nMat*d)]
    alpha.hat <- beta1[1]

    gs.hat <- matrix(0,n,d)
    correc <- rep(0,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      correc[ell] <- mean(aux)
      gs.hat[,ell] <- aux - mean(aux)
    }

    is.zero <- c(abs(alpha.hat)<bound.control,abs(coef.lin)<bound.control,normgammaj<bound.control)
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, betas=beta1, coef.const=alpha.hat, coef.lin=coef.lin, coef.spl=coef.spl, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, xdesign=xdesign, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, normgammaj=normgammaj, is.zero=is.zero)
    return(salida)
  }

#' Variable selection in classical PLAM with fixed lambdas
#' @export
plam.cl.vs.lambdas <- function(y, Z, X, lambdas1, lambdas2, nknots=NULL, knots=NULL, degree.spline=3, MAXITER=100, bound.control=10^(-3)){
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
    AUX <- select.nknots.cl(y, Z, X, degree.spline=degree.spline)
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
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]
  sal <- stats::lm(y~Z.aux+Xspline)
  xdesign <- cbind(rep(1,n),Z.aux,Xspline)
  #Construyo el beta 0
  beta.ini <- as.vector(sal$coefficients)
  #beta.hat <- beta.ini[-1]
  #coef.lin <- beta.ini[2:(q+1)]
  #coef.spl <- beta.ini[(q+2):(1+q+nMat*d)]
  nbetas <- length(beta.ini)

  normgammaj <- rep(0,d)
  corte <- 1
  iter <- 0
  while( (corte>bound.control) & (iter<MAXITER)){
    iter <- iter +1
    #print(iter)
    regresion.hat <- xdesign%*%beta.ini
    res <- y-regresion.hat
    an <- quantile(abs(res),2*(n^(-1/2)))
    W <- diag( as.vector((1/2)*(res)^2/(an+(res)^2) ))

    Sigmalambda <- matrix(0,nbetas,nbetas)
    for(i in 1:(q+1)){
      Sigmalambda[i,i] <- scad.d(beta.ini[i],lambda=lambdas1[i])/abs(beta.ini[i])
    }
    for(i in 1:d){
      Hj <- Hj.matrix(X[,i], nknots, degree.spline)
      gammaj <- as.matrix(beta.ini[(q+2+nMat*(i-1)):(nMat*i+(q+1))])
      normgammaj[i] <- sqrt( t(gammaj)%*%Hj%*%gammaj )
      Sigmalambda[(q+2+nMat*(i-1)):(nMat*i+(q+1)),(q+2+nMat*(i-1)):(nMat*i+(q+1))] <- as.numeric(scad.d(normgammaj[i],lambda=lambdas2[i])/(normgammaj[i]))*Hj
    }

    #Sigmalambda <- diag( scad.d(beta.ini,lambda=lambda)/abs(beta.ini))

    #options(show.error.messages = TRUE)
    try.sal <- try(
      AUX <- solve(t(xdesign)%*%W%*%xdesign + 1/2*n*Sigmalambda)
    )

    if(class(try.sal)!= 'try-error'){
      beta1 <- as.vector(AUX%*%t(xdesign)%*%W%*%y)
      corte <- my.norm.2(beta.ini-beta1)/my.norm.2(beta.ini)
      beta.ini <- beta1
    }else{
      beta1 <- beta.ini
      iter <- MAXITER
    }
  }

  beta.hat <- beta1[-1]
  coef.lin <- beta1[2:(q+1)]
  coef.spl <- beta1[(q+2):(1+q+nMat*d)]
  alpha.hat <- beta1[1]

  gs.hat <- matrix(0,n,d)
  correc <- rep(0,d)
  for(ell in 1:d){
    aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    correc[ell] <- mean(aux)
    gs.hat[,ell] <- aux - mean(aux)
  }

  is.zero <- c(abs(alpha.hat)<bound.control,abs(coef.lin)<bound.control,normgammaj<bound.control)
  salida <- list(prediction=regresion.hat, betas=beta1, coef.const=alpha.hat, coef.lin=coef.lin, coef.spl=coef.spl, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, xdesign=xdesign, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, normgammaj=normgammaj, is.zero=is.zero)
  return(salida)
}


#' Hj matrix
#' @export
Hj.matrix <- function(Xj, nknots, degree.spline){
  if(nknots>0){
     knots <- stats::quantile(Xj,(1:nknots)/(nknots+1))
  }else{
    knots <- NULL
  }
  rango <- range(Xj)
  grilla <- seq(rango[1],rango[2],length=1000)
  Mat <- splines::bs(grilla, knots=knots, degree=degree.spline, intercept=FALSE)
  #dim(Mat)
  return( (t(Mat)%*%Mat)/length(grilla) )
}

#' Selection lambdas with robust BIC criteria
#' @export
select.rob.lambdas <- function(y, Z, X, grid.lambda1, grid.lambda2, nknots=NULL, knots=NULL, degree.spline=3, maxit=100, method="MM", MAXITER=100){
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
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]
  sal <- MASS::rlm(y~Z.aux+Xspline ,method=method, maxit=maxit)
  xdesign <- sal$x
  sigma.hat <- sal$s
  #Construyo el beta 0
  beta0 <- beta.ini <- as.vector(sal$coefficients)
  #beta.hat <- beta.ini[-1]
  #coef.lin <- beta.ini[2:(q+1)]
  #coef.spl <- beta.ini[(q+2):(1+q+nMat*d)]
  nbetas <- length(beta.ini)
  normgammaj0 <- rep(0,d)
  for(i in 1:d){
    Hj <- Hj.matrix(X[,i], nknots, degree.spline)
    gammaj <- as.matrix(beta.ini[(q+2+nMat*(i-1)):(nMat*i+(q+1))])
    normgammaj0[i] <- sqrt( t(gammaj)%*%Hj%*%gammaj )
  }

  grilla <- expand.grid(grid.lambda1,grid.lambda2)
  dim.grilla <- dim(grilla)[1]
  BIC <- rep(0,dim.grilla)
  for(i in 1:dim.grilla){
    #print(i)
    lambdas1 <- rep(grilla[i,1],q+1)
    lambdas2 <- rep(grilla[i,2],d)
    AUX2 <- plam.rob.vs.lambdas(y, Z, X, lambdas1, lambdas2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=maxit, method=method, MAXITER=MAXITER)
    betas <- AUX2$betas
    nbasis <- AUX2$nbasis
    xdesign <- AUX2$xdesign
    desvio.hat <- AUX2$sigma.hat
    normgammaj <- AUX2$normgammaj
    dfc <- sum( abs(betas[1:(q+1)])>10^(-3))
    dfn <- sum( normgammaj >10^(-3))
    regresion.hat <- xdesign%*%betas
    tuk <- tukey.loss( (y - regresion.hat)/desvio.hat )
    #tuk <- tukey.loss( (y - regresion.hat)/desvio.hat )*desvio.hat^2 #Este funciona peor
    BIC[i] <- mean(tuk) + dfn*(log(n/nbasis)/(n/nbasis)) + dfc*(log(n)/n)
  }

  position <- which.min(BIC)
  la1 <- grilla[position,1]
  la2 <- grilla[position,2]

  lambdas1 <- la1/abs(beta0[1:(q+1)])
  lambdas2 <- la2/normgammaj0

  salida <- list(la1=la1, la2=la2, lambdas1=lambdas1, lambdas2=lambdas2)
  return(salida)
}


#' Selection lambdas with classical BIC criteria
#' @export
select.cl.lambdas <- function(y, Z, X, grid.lambda1, grid.lambda2, nknots=NULL, knots=NULL, degree.spline=3, MAXITER=100){
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
    AUX <- select.nknots.cl(y, Z, X, degree.spline=degree.spline)
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
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]
  sal <- stats::lm(y~Z.aux+Xspline)
  xdesign <- cbind(rep(1,n),Z.aux,Xspline)
  #Construyo el beta 0
  beta0 <- beta.ini <- as.vector(sal$coefficients)
  #beta.hat <- beta.ini[-1]
  #coef.lin <- beta.ini[2:(q+1)]
  #coef.spl <- beta.ini[(q+2):(1+q+nMat*d)]
  nbetas <- length(beta.ini)
  normgammaj0 <- rep(0,d)
  for(i in 1:d){
    Hj <- Hj.matrix(X[,i], nknots, degree.spline)
    gammaj <- as.matrix(beta.ini[(q+2+nMat*(i-1)):(nMat*i+(q+1))])
    normgammaj0[i] <- sqrt( t(gammaj)%*%Hj%*%gammaj )
  }

  grilla <- expand.grid(grid.lambda1,grid.lambda2)
  dim.grilla <- dim(grilla)[1]
  BIC <- rep(0,dim.grilla)
  for(i in 1:dim.grilla){
    #print(i)
    lambdas1 <- rep(grilla[i,1],q+1)
    lambdas2 <- rep(grilla[i,2],d)
    AUX2 <- plam.cl.vs.lambdas(y, Z, X, lambdas1, lambdas2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER)
    betas <- AUX2$betas
    nbasis <- AUX2$nbasis
    xdesign <- AUX2$xdesign
    normgammaj <- AUX2$normgammaj
    dfc <- sum( abs(betas[1:(q+1)])>10^(-3))
    dfn <- sum( normgammaj >10^(-3))
    regresion.hat <- xdesign%*%betas
    tuk <- (1/2)*( (y - regresion.hat) )^2
    BIC[i] <- mean(tuk) + dfn*(log(n/nbasis)/(n/nbasis)) + dfc*(log(n)/n)
  }

  position <- which.min(BIC)
  la1 <- grilla[position,1]
  la2 <- grilla[position,2]

  lambdas1 <- la1/abs(beta0[1:(q+1)])
  lambdas2 <- la2/normgammaj0

  salida <- list(la1=la1, la2=la2, lambdas1=lambdas1, lambdas2=lambdas2)
  return(salida)
}


#' Robust Partial Linear Additive Model with Variable Selection
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
#' @export
plam.rob.vs <- function(y, Z, X, vs=TRUE, nknots=NULL, knots=NULL, degree.spline=3, rob.maxit=100, method="MM", MAXITER=100, bound.control=10^(-3)){
  if(vs=="TRUE"){
    grid.lambda1 <- seq(0,0.5,0.1)
    grid.lambda2 <- seq(0,0.5,0.1)
    sal <- select.rob.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method, MAXITER=MAXITER)
    lambda1 <- sal$la1
    lambda2 <- sal$la2
    if(lambda1==0.5){
     if(lambda2==0.5){
       grid.lambda1 <- seq(0.5,1,0.1)
       grid.lambda2 <- seq(0.5,1,0.1)
       sal <- select.rob.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method, MAXITER=MAXITER)
       lambda1 <- sal$la1
       lambda2 <- sal$la2
     }else{
       grid.lambda1 <- seq(0.5,1,0.1)
       grid.lambda2 <- seq(0,0.5,0.1)
       sal <- select.rob.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method, MAXITER=MAXITER)
       lambda1 <- sal$la1
       lambda2 <- sal$la2
     }
    }else{
       if(lambda2==0.5){
         grid.lambda1 <- seq(0,0.5,0.1)
         grid.lambda2 <- seq(0.5,1,0.1)
         sal <- select.rob.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method, MAXITER=MAXITER)
         lambda1 <- sal$la1
         lambda2 <- sal$la2
       }
    }
    #print(lambda1)
    #print(lambda2)
    lambdas1 <- sal$lambdas1
    lambdas2 <- sal$lambdas2
    sal <- plam.rob.vs.lambdas(y=y, Z=Z, X=X, lambdas1=lambdas1, lambdas2=lambdas2, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method, MAXITER=MAXITER, bound.control = bound.control)
    salida <- c(sal,lambda1=list(lambda1), lambda2=list(lambda2))
    return(salida)
  }else{
    sal <- plam.rob(y=y, Z=Z, X=X, nknots=nknots, knots=knots, degree.spline=degree.spline, maxit=rob.maxit, method=method)
    return(sal)
  }
}


#' Classical Partial Linear Additive Model with Variable Selection
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
#' @export
plam.cl.vs <- function(y, Z, X, vs=TRUE, nknots=NULL, knots=NULL, degree.spline=3, MAXITER=100, bound.control=10^(-3)){
  if(vs=="TRUE"){
    grid.lambda1 <- seq(0,0.5,0.1)
    grid.lambda2 <- seq(0,0.5,0.1)
    sal <- select.cl.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER)
    lambda1 <- sal$la1
    lambda2 <- sal$la2
    if(lambda1==0.5){
      if(lambda2==0.5){
        grid.lambda1 <- seq(0.5,1,0.1)
        grid.lambda2 <- seq(0.5,1,0.1)
        sal <- select.cl.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER)
        lambda1 <- sal$la1
        lambda2 <- sal$la2
      }else{
        grid.lambda1 <- seq(0.5,1,0.1)
        grid.lambda2 <- seq(0,0.5,0.1)
        sal <- select.cl.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER)
        lambda1 <- sal$la1
        lambda2 <- sal$la2
      }
    }else{
      if(lambda2==0.5){
        grid.lambda1 <- seq(0,0.5,0.1)
        grid.lambda2 <- seq(0.5,1,0.1)
        sal <- select.cl.lambdas(y=y, Z=Z, X=X, grid.lambda1=grid.lambda1, grid.lambda2=grid.lambda2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER)
        lambda1 <- sal$la1
        lambda2 <- sal$la2
      }
    }
    #print(lambda1)
    #print(lambda2)
    lambdas1 <- sal$lambdas1
    lambdas2 <- sal$lambdas2
    sal <- plam.cl.vs.lambdas(y=y, Z=Z, X=X, lambdas1=lambdas1, lambdas2=lambdas2, nknots=nknots, knots=knots, degree.spline=degree.spline, MAXITER=MAXITER, bound.control = bound.control)
    salida <- c(sal,lambda1=list(lambda1), lambda2=list(lambda2))
    return(salida)
  }else{
    sal <- plam.cl(y=y, Z=Z, X=X, nknots=nknots, knots=knots, degree.spline=degree.spline)
    return(sal)
  }
}


