#' Derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' @param r a vector of real numbers
#' @param k a positive tuning constant.
#'
#' @return A vector of the same length as \code{x}.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.tukey(r=x, k = 1.5)
#'
#' @export
#' @importFrom stats lm
#' @importFrom fda getbasismatrix
#' @importFrom fda create.bspline.basis
#' @importFrom robustbase lmrob
#' @importFrom robustbase lmrob.control
#' @importFrom splines bs
#' @useDynLib rplam, .registration = TRUE
psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w)
}

#' Derivative of the Tukey Loss Function
#' @export
psi.tukey.derivative <- function(r, k=4.685){
  u <- abs(r/k)
  w <- ((1-u)*(1+u))^2-4*u^2*((1-u)*(1+u))
  w[u>1] <- 0
  return(w)
}

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

# Tukey's weight function "Psi(r)/r"
psi.w <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}

#' Euclidean norm of a vector
#'
#' This function calculates the Euclidean norm of a vector.
#'
#' @param x A real vector.
#'
#' @return The Euclidean norm of the input vector.
#'
#' @author Alejandra Mercedes Martinez, \email{ammartinez@conicet.gov.ar}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' my.norm.2(x)
#'
#' @export
my.norm.2 <- function(x){
  return( sqrt(sum(x^2)) )
}


#' Selection of the number of knots for the least-squares estimator under a PLAM
#'
#' This function automatically selects the number of internal knots for the B-spline approximation. It uses the BIC criterion with the quadratic function.
#'
#' @param y a vector of real numbers.
#' @param Z a matrix of numbers corresponding to the covariates entering linearly into the model.
#' @param X a matrix of numbers corresponding to the covariates entering in the additive component of the model.
#' @param degree.spline spline degree. Defaults to \code{'3'}.
#'
#' @return A list with the following components:
#' \item{nknots}{Number of internal knots selected by the procedure.}
#' \item{grid.nknots}{Grid of internal knots used.}
#' \item{BIC}{Values of the BIC criterion for each element of the grid.}
#' \item{kj}{Number of elements of the B-spline basis used to approximate each additive function. It is calculated as \code{nknots + degree.spline}.}
#' \item{nbasis}{Total number of elements of the basis of B-splines. This corresponds to \code{d* kj} where \code{d} is the number of covariates entering in the additive part.
#'
#' @references
#' Boente G. and Martinez A. (2023). A robust spline approach in partially linear additive models. Computational Statistics and Data Analysis, 178, 107611.
#'
#' \author Alejandra Martinez, \email{ammartinez@conicet.gov.ar}
#'
#' @examples
#' set.seed(11)
#' n <- 100
#' z1 <- rnorm(n)
#' z2 <- rbinom(n, 4, 1/2)
#' x1 <- runif(n,-1,1)
#' x2 <- runif(n,-1,1)
#' err <- rnorm(n, 0, 0.1)
#' regre <- 2+3*z1-4*z2+x1^3+2*sin(pi*x2)
#' y <- regre + err
#' Z <- cbind(z1,z2)
#' X <- cbind(x1,x2)
#' sal <- select.nknots.cl(y, Z, X)
#'
#' @export
select.nknots.cl <- function(y, Z, X, degree.spline = 3){
  n <- length(y)
  d <- dim(X)[2]

  r <- degree.spline-1
  if(r<1){
    cat("No se cumple la hipótesis de: 1<=r")
  }

  #Corregir para que solo tire un warning
  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  lim.inf.kj <- ceiling(max(n^(1/(2*r+1))/2,degree.spline+1))
  lim.sup.kj <- floor(8+2*n^(1/(2*r+1)))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]

      #Centered with the integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])

    }
    nMat <- dim(Mat.X[[1]])[2]

    #Computing the least-squares regression estimator
    sal <- stats::lm(y~Z.aux+Xspline)

    betas <- as.vector(sal$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    }

    regresion.hat <- stats::predict(sal) #which is also alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline)

    #Computing the BIC criterion
    BIC[nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+q+1) #q+1 is the amount of linear components
  }

  posicion <- which.min(BIC)
  nknots <- posicion-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots, nbasis = nbasis, kj=kj)

  return(salida)
}


#' Selection of the number of knots for the least-squares estimator for additive models
#'
#' This function automatically selects the number of internal knots for the B-spline approximation. It uses the BIC criterion with the quadratic function.
#'
#' @param y a vector of real numbers.
#' @param X a matrix of numbers corresponding to the covariates entering in the additive component of the model.
#' @param degree.spline spline degree. Defaults to \code{'3'}.
#'
#' @return A list with the following components:
#' \item{nknots}{Number of internal knots selected by the procedure.}
#' \item{grid.nknots}{Grid of internal knots used.}
#' \item{BIC}{Values of the BIC criterion for each element of the grid.}
#' \item{kj}{Number of elements of the B-spline basis used to approximate each additive function. It is calculated as \code{nknots + degree.spline}.}
#' \item{nbasis}{Total number of elements of the basis of B-splines. This corresponds to \code{d* kj} where \code{d} is the number of covariates entering in the additive part.
#'
#' @references
#' Boente G. and Martinez A. (2023). A robust spline approach in partially linear additive models. Computational Statistics and Data Analysis, 178, 107611.
#'
#' \author Alejandra Martinez, \email{ammartinez@conicet.gov.ar}
#'
#' @examples
#' set.seed(11)
#' n <- 100
#' x1 <- runif(n,-1,1)
#' x2 <- runif(n,-1,1)
#' err <- rnorm(n, 0, 0.1)
#' regre <- 2+x1^3+2*sin(pi*x2)
#' y <- regre + err
#' X <- cbind(x1,x2)
#' sal <- select.nknots.cl.am(y, X)
#'
#' @export
select.nknots.cl.am <- function(y, X, degree.spline=3){
  q <- 0
  n <- length(y)
  d <- dim(X)[2]

  r <- degree.spline-1
  if(r<1){
    cat("No se cumple la hipótesis de: 1<=r")
  }

  lim.inf.kj <- ceiling(max(n^(1/(2*r+1))/2, degree.spline+1))
  lim.sup.kj <- floor(8+2*n^(1/(2*r+1)))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){

    Mat.X <- as.list(rep(0,d))
    Xspline <- NULL
    for (ell in 1:d){
      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]

      #Centered with the integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])

    }
    nMat <- dim(Mat.X[[1]])[2]

    sal <- stats::lm(y~Xspline)
    betas <- as.vector(sal$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    }

    regresion.hat <- stats::predict(sal)
    nbasis <- d*(nknots + degree.spline)

    BIC[nknots-lim.inf.nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+q+1)
  }
  posicion <- which.min(BIC)
  nknots <- posicion+lim.inf.nknots-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots, nbasis = nbasis, kj=kj)

  return(salida)
}



#' Robust knot selection for plam
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
select.nknots.rob <- function(y, Z, X, degree.spline=3, maxit=100){

  r <- degree.spline-1
  if(r<1){
    cat("No se cumple la hipótesis de: 1<=r<=ell-2")
  }
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


  lim.inf.kj <- ceiling(max(n^(1/(2*r+1))/2,degree.spline+1))
  lim.sup.kj <- floor(8+2*n^(1/(2*r+1)))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  RBIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])

    }
    nMat <- dim(Mat.X[[1]])[2] #Decía ell
    #dim(Xspline)[2]/4



    #- Tukey MM estimator -#
    control <- robustbase::lmrob.control(trace.level = 0,         # 0
                             nResample   =  500,      # 500 default
                             tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                             subsampling = 'simple',  #
                             rel.tol     = 1e-5,      # 1e-7
                             refine.tol  = 1e-5,      # 1e-7
                             k.max       = 2e3,       # 200
                             maxit.scale = maxit,       # 200 #2e3
                             max.it      = maxit)       # 50 #2e3
    tmp <- try(robustbase::lmrob(y ~ Z.aux+Xspline, control = control))
    #sal.r  <- robustbase::lmrob(y ~ Z.aux+Xspline, control = control)
    if( class(tmp)[1] != 'try-error'){ #Comentar esto si no funciona y sacar la def de tmp
      sal.r <- tmp
    betas <- as.vector(sal.r$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      #gs.hat[,ell] <- aux #- mean(aux) No necesita porque los splines están centrados
      gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    }

    regresion.hat.r <- stats::predict(sal.r) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    desvio.hat <- sal.r$s
    tuk <- tukey.loss( (y - regresion.hat.r)/desvio.hat )
    RBIC[nknots-lim.inf.nknots+1] <- log( (desvio.hat^2)*sum(tuk) )+ (log(n)/(2*n))*(nbasis+q+1) #q+1 porque q de la parte lineal y 1 de la constante. O sea, q+1 es la cantidad de lineales.
    }else{
      RBIC[nknots-lim.inf.nknots+1] <- NA
    }

  }
  posicion <- which.min(RBIC)
  nknots <- posicion+lim.inf.nknots-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, RBIC=RBIC, grid.nknots=grid.nknots, nbasis = nbasis, kj = kj)
  return(salida)


}


#' Robust knot selection for plam
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
select.nknots.rob.am <- function(y, X, degree.spline=3, maxit=100){

  n <- length(y)
  d <- dim(X)[2]
  q <- 0

  r <- degree.spline-1
  if(r<1){
    cat("No se cumple la hipótesis de: 1<=r<=ell-2")
  }

  lim.inf.kj <- ceiling(max(n^(1/(2*r+1))/2, degree.spline+1))
  lim.sup.kj <- floor(8+2*n^(1/(2*r+1)))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  RBIC <- rep(0,length(grid.nknots))


  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    Xspline <- NULL
    for (ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X[[ell]] <- spl.final[,-1]

      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[1]])[2] #Decía ell
    #dim(Xspline)[2]/4



    #- Tukey MM estimator -#
    control <- robustbase::lmrob.control(trace.level = 0,         # 0
                             nResample   =  500,      # 500 default
                             tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                             subsampling = 'simple',  #
                             rel.tol     = 1e-5,      # 1e-7
                             refine.tol  = 1e-5,      # 1e-7
                             k.max       = 2e3,       # 200
                             maxit.scale = maxit,       # 200 #2e3
                             max.it      = maxit)       # 50 #2e3
    tmp <- try(robustbase::lmrob(y ~ Xspline, control = control))
    #sal.r  <- robustbase::lmrob(y ~ Z.aux+Xspline, control = control)
    if( class(tmp)[1] != 'try-error'){ #Comentar esto si no funciona y sacar la def de tmp
      sal.r <- tmp
      betas <- as.vector(sal.r$coefficients)
      beta.hat <- betas[-1]
      coef.lin <- betas[2:(q+1)]
      coef.spl <- betas[(q+2):(1+q+nMat*d)]
      alpha.hat <- betas[1]

      gs.hat <- matrix(0,n,d)
      for(ell in 1:d){
        #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #gs.hat[,ell] <- aux #- mean(aux) No se necesita porque están centrados con la integral
        gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }

      regresion.hat.r <- stats::predict(sal.r) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

      nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
      desvio.hat <- sal.r$s
      tuk <- tukey.loss( (y - regresion.hat.r)/desvio.hat )
      RBIC[nknots-lim.inf.nknots+1] <- log( (desvio.hat^2)*sum(tuk) )+ (log(n)/(2*n))*(nbasis+q+1) #q+1 porque q de la parte lineal y 1 de la constante. O sea, q+1 es la cantidad de lineales.
    }else{
      RBIC[nknots-lim.inf.nknots+1] <- NA
    }

  }
  posicion <- which.min(RBIC)
  nknots <- posicion+lim.inf.nknots-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, RBIC=RBIC, grid.nknots=grid.nknots, nbasis = nbasis, kj = kj)
  return(salida)


}


#' Classical Partial Linear Additive Model
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

    grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

    if(nknots>0){
      aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
    }else{
      nodos.spl <- c(min(X[,ell]), max(X[,ell]))
    }

    #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                        norder = (degree.spline+1),
                                        breaks = nodos.spl)
    aux <- fda::getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])

  }
  nMat <- dim(Mat.X[[1]])[2] #Decía ell

  sal <- stats::lm(y~Z.aux+Xspline)
  betas <- as.vector(sal$coefficients)
  #Arreglo de NA's
  betas[is.na(betas)] <- rep(0,sum(is.na(betas)))

  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    #correc[ell] <- mean(aux) #Esto ya no lo necesito porque integran 0
    #gs.hat[,ell] <- aux - mean(aux)
    gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj)
      #list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- X.new <- as.matrix(np.point)
      }else{
        prediccion <- X.new <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- X.new <- np.point
    }

    np <- dim(X.new)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X.new[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])

    }


    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #prediccion[,ell] <- aux - correc[ell] #Esto ya no lo necesito
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
      #list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Classical Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
am.cl <- function(y, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  n <- length(y)
  d <- dim(X)[2]
  q <- 0

  if( is.null(nknots) ){
    AUX <- select.nknots.cl.am(y,X,degree.spline=degree.spline)
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
    grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

    if(nknots>0){
      aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
    }else{
      nodos.spl <- c(min(X[,ell]), max(X[,ell]))
    }

    #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                        norder = (degree.spline+1),
                                        breaks = nodos.spl)
    aux <- fda::getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[1]])[2] #Decía ell

  sal <- stats::lm(y~Xspline)
  betas <- as.vector(sal$coefficients)
  #Arreglo de NA's
  betas[is.na(betas)] <- rep(0,sum(is.na(betas)))

  beta.hat <- betas[-1]
  #coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    if(length((nMat*(ell-1)+1):(nMat*ell))!=1){
      #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      #correc[ell] <- mean(aux)
      #gs.hat[,ell] <- aux - mean(aux)
      gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    }else{
      #aux <- Xspline[,(nMat*(ell-1)+1):(nMat*ell)] * coef.spl[(nMat*(ell-1)+1):(nMat*ell)]
      #correc[ell] <- mean(aux)
      #gs.hat[,ell] <- aux - mean(aux) #No lo necesito porque integran 0
      gs.hat[,ell] <- Xspline[,(nMat*(ell-1)+1):(nMat*ell)] * coef.spl[(nMat*(ell-1)+1):(nMat*ell)]
    }
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj)
      #list(prediction=regresion.hat, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- X.new <- as.matrix(np.point)
      }else{
        prediccion <- X.new <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- X.new <- np.point
    }
    np <- dim(X.new)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X.new[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])

    }


    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #prediccion[,ell] <- aux - correc[ell] #Ya integran 0
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
      #list(prediction=regresion.hat, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Robust Partial Linear Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
plam.rob <- function(y, Z, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3, maxit=100){
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
    AUX <- select.nknots.rob(y, Z, X, degree.spline=degree.spline, maxit=maxit)
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
    grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

    if(nknots>0){
      aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
    }else{
      nodos.spl <- c(min(X[,ell]), max(X[,ell]))
    }

    #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                        norder = (degree.spline+1),
                                        breaks = nodos.spl)
    aux <- fda::getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])

  }
  nMat <- dim(Mat.X[[1]])[2] #Decía ell

  control <- robustbase::lmrob.control(trace.level = 0,         # 0
                           nResample   =  500,      # 500 default
                           tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                           subsampling = 'simple',  #
                           rel.tol     = 1e-5,      # 1e-7
                           refine.tol  = 1e-5,      # 1e-7
                           k.max       = 2e3,       # 200
                           maxit.scale = maxit,       # 200 #2e3
                           max.it      = maxit)       # 50 #2e3
  sal  <- robustbase::lmrob(y ~ Z.aux+Xspline, control = control)
  #MASS::rlm(y~Z.aux+Xspline ,method=method,maxit=maxit)
  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]
  sigma.hat <- sal$s

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    #correc[ell] <- mean(aux)
    #gs.hat[,ell] <- aux - mean(aux) #Esto ya no lo necesito porque integran 0
    gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj)
      #list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- X.new <- as.matrix(np.point)
      }else{
        prediccion <- X.new <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- X.new <- np.point
    }
    np <- dim(X.new)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){

      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X.new[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }


    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #prediccion[,ell] <- aux - correc[ell] #Esto ya no lo necesito
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
      #list(prediction=regresion.hat, sigma.hat=sigma.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Robust Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
am.rob <- function(y, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3, maxit=100){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  n <- length(y)
  d <- dim(X)[2]
  q <- 0

  if( is.null(nknots) ){
    AUX <- select.nknots.rob.am(y, X, degree.spline=degree.spline, maxit=maxit)
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
    grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

    if(nknots>0){
      aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
    }else{
      nodos.spl <- c(min(X[,ell]), max(X[,ell]))
    }

    #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                             norder = (degree.spline+1),
                                             breaks = nodos.spl)
    aux <- fda::getbasismatrix(X[,ell], base.beta)
    naux <- dim(aux)[2]
    #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

    #Centrado con la integral
    spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
    spl.final <- aux
    for (j in 1:naux){
      centroj=mean(spl.center[,j])
      spl.final[,j]=aux[,j]-centroj
    }
    Mat.X[[ell]] <- spl.final[,-1]

    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[1]])[2] #Decía ell

  control <- robustbase::lmrob.control(trace.level = 0,         # 0
                           nResample   =  500,      # 500 default
                           tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
                           subsampling = 'simple',  #
                           rel.tol     = 1e-5,      # 1e-7
                           refine.tol  = 1e-5,      # 1e-7
                           k.max       = 2e3,       # 200
                           maxit.scale = maxit,       # 200 #2e3
                           max.it      = maxit)       # 50 #2e3
  sal  <- robustbase::lmrob(y ~ Xspline, control = control)
  #MASS::rlm(y~Z.aux+Xspline ,method=method,maxit=maxit)
  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  #coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]
  sigma.hat <- sal$s

  gs.hat <- matrix(0,n,d)
  #correc <- rep(0,d)
  for(ell in 1:d){
    if(length((nMat*(ell-1)+1):(nMat*ell))!=1){
      #aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      #correc[ell] <- mean(aux)
      #gs.hat[,ell] <- aux - mean(aux) #Esto ya no lo necesito porque integran 0
      gs.hat[,ell] <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    }else{
      #aux <- Xspline[,(nMat*(ell-1)+1):(nMat*ell)] * coef.spl[(nMat*(ell-1)+1):(nMat*ell)]
      #correc[ell] <- mean(aux)
      gs.hat[,ell] <- Xspline[,(nMat*(ell-1)+1):(nMat*ell)] * coef.spl[(nMat*(ell-1)+1):(nMat*ell)]
    }
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj)
      #list(prediction=regresion.hat, sigma.hat=sigma.hat, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- X.new <- as.matrix(np.point)
      }else{
        prediccion <- X.new <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- X.new <- np.point
    }
    np <- dim(X.new)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){
      grilla.tes <- seq(min(X[,ell]),max(X[,ell]),length=n)

      if(nknots>0){
        aa <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
        nodos.spl <- c(min(X[,ell]), aa, max(X[,ell]))
      }else{
        nodos.spl <- c(min(X[,ell]), max(X[,ell]))
      }

      #Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      base.beta   <- fda::create.bspline.basis(rangeval = c(min(X[,ell]), max(X[,ell])),
                                          norder = (degree.spline+1),
                                          breaks = nodos.spl)
      aux <- fda::getbasismatrix(X.new[,ell], base.beta)
      naux <- dim(aux)[2]
      #Mat.X[[ell]] <- aux-t(matrix(colMeans(aux),naux,n))

      #Centrado con la integral
      spl.center   <- fda::getbasismatrix(grilla.tes, base.beta)
      spl.final <- aux
      for (j in 1:naux){
        centroj=mean(spl.center[,j])
        spl.final[,j]=aux[,j]-centroj
      }
      Mat.X.new[[ell]] <- spl.final[,-1]

      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])

    }


    for(k in 1:np){
      for(ell in 1:d){
        #aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        #prediccion[,ell] <- aux - correc[ell] #Ya no lo necesito porque integran 0
        prediccion[,ell] <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      }
    }
    salida <- list(prediction=regresion.hat, sigma.hat=sigma.hat, g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
      #list(prediction=regresion.hat, sigma.hat=sigma.hat, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


pos.est <- function(y, sigma.hat, typePhi, ini=NULL, epsilon=1e-6, iter.max=10){
  yp <- y[ tmp<-!is.na(y) ]
  if(is.null(ini)){
    ini <- median(yp)
  }
  corte <- 10
  iter <- 0
  n <- length(yp)
  prob <- rep(1,n)
  k.h <- 1.345
  k.t <- 4.685
  if(typePhi=='Huber'){
    beta <- .C("huber_pos", as.integer(n), as.double(yp), as.double(ini), as.double(epsilon),
               as.double(sigma.hat), as.double(prob), as.double(k.h), as.integer(iter.max), salida=as.double(0) )$salida
  }
  if(typePhi=='Tukey'){
    beta <- .C("tukey_pos", as.integer(n), as.double(yp), as.double(ini), as.double(epsilon),
               as.double(sigma.hat), as.double(prob), as.double(k.t), as.integer(iter.max), salida=as.double(0) )$salida
  }
  return(beta)
}
