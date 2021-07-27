source("rplam-simu-fn.R")



#library(splines)

huber.loss <- function(x,k=1.345){
  n <- length(x)
  salida <- rep(0,n)
  for(i in 1:n){
    if(abs(x[i])<=k){salida[i] <- (x[i])^2
    }else{salida[i] <- 2*k*abs(x[i])-k^2}
  }
  return(salida)
}

tukey.loss <- function(x,k=4.685){
  n <- length(x)
  salida <- rep(0,n)
  for(i in 1:n){
    if(abs(x[i])<=k){salida[i] <- (x[i])^6/(6*k^4)-(x[i])^4/(2*k^2)+(x[i])^2/2
    }else{salida[i] <- k^2/6}
  }
  return(salida)
}

tms <- function(a, alpha=.1) {
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}


simulacion<-function(ejemplo, NITER,anterior,tipo.cont,n,desvio.epsilon,epsilon,degree.spline){
  
  AUX <- 0
  
  seed <- 123
  
  coef.rob <- matrix(0,NITER+anterior,2)
  ecm.g1.rob <- rep(0,NITER+anterior)
  ecm.g2.rob <- rep(0,NITER+anterior)
  ecm.global.rob <- rep(0,NITER+anterior)
  ecm.g1.rob.t <- rep(0,NITER+anterior)
  ecm.g2.rob.t <- rep(0,NITER+anterior)
  ecm.global.rob.t <- rep(0,NITER+anterior)
  ecm.g1.rob.tt <- rep(0,NITER+anterior)
  ecm.g2.rob.tt <- rep(0,NITER+anterior)
  ecm.global.rob.tt <- rep(0,NITER+anterior)
  
  ecm.g1.rob.int <- rep(0,NITER+anterior)
  ecm.g2.rob.int <- rep(0,NITER+anterior)
  ecm.g1.rob.int.t <- rep(0,NITER+anterior)
  ecm.g2.rob.int.t <- rep(0,NITER+anterior)
  ecm.g1.rob.int.tt <- rep(0,NITER+anterior)
  ecm.g2.rob.int.tt <- rep(0,NITER+anterior)
  
  mseyback.rob <- matrix(0,NITER+anterior,6)
  
  inicio<-anterior+1
  final<-anterior+NITER
  
  nombreRob<-paste("PLAM_splines_new_Clasico_ej", ejemplo,"_C",tipo.cont,"desde_",inicio,"hasta_",final,"_newm3.txt",sep="")
  nombre.g1<-paste("PLAM_splines_new_Clasico_ej", ejemplo,"_C",tipo.cont,"desde_",inicio,"hasta_",final,"_g1_newm3.txt",sep="")
  nombre.g2<-paste("PLAM_splines_new_Clasico_ej", ejemplo,"_C",tipo.cont,"desde_",inicio,"hasta_",final,"_g2_newm3.txt",sep="")
  
  for(iter in inicio:final){
    
    print(c("iter",iter))
    
    #############################
    # FIJO LA SEMILLA
    #############################
    seed.new <- seed + 17*iter
    set.seed(seed.new)
    
    
    #############################
    # EJEMPLO 1
    #############################
    
    if( ejemplo==1 ){
      
     
      #Generamos las covariables
      mu <- c(0,0,1/2,1/2) #c(0,0,1/2,-1/2)    
      Sigma <- matrix( c(1,0, 1/(2*sqrt(3)), 1/(6*sqrt(3)), 0, 1, 1/(2*sqrt(3)), 1/(6*sqrt(3)) , 1/(2*sqrt(3)), 1/(2*sqrt(3)), 1/4, 0, 1/(6*sqrt(3)), 1/(6*sqrt(3)), 0, 1/6), 4, 4)
      
      #Si no tengo el rmvtnorm:
      U <- eigen(Sigma)$vectors
      Lambda <- matrix(0, 4,4)
      diag(Lambda) <- sqrt(eigen(Sigma)$values)
      
      RaizCuadradaSigma <- U%*%Lambda%*%(t(U))
      Z <- matrix(rnorm(4*n), 4,n)
      MediasMu <- t( matrix( rep(mu, n*4), 4, n))
      covariables <- t( RaizCuadradaSigma%*%Z ) + MediasMu
      
      #covariables <- as.matrix(rmvnorm(n, mu, Sigma))
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones  
      funcion.g1 <- function(t1){
        return( 2*sin(pi*(t1-1/2)) )
      }
      funcion.g2 <- function(t2){
        return(2*cos(pi*(t2-1/2-pi/2))) #return( 2*cos(4*pi*t2) )
      }

    }
    
    #############################
    # EJEMPLO 2
    #############################
    
    if( ejemplo==2 ){
      
      #Generamos las covariables
      #mu <- rep(0,4)    
      #Sigma <- matrix(0,4,4) 
      #diag(Sigma) <- rep(1,4)
      #covariables <- as.matrix(rmvnorm(n, mu, Sigma))
      covariables <- matrix( rnorm(n*4), n, 4)
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1))
      }
      funcion.g2 <- function(t2){
        return(1/2*t2)
      }
    }
    
    
    #############################
    # EJEMPLO 3
    #############################
    
    if( ejemplo==3 ){
      
      #Generamos las covariables
      covariables <- matrix(runif(n*4,-1,1),n,4)
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1))
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-1/2*(exp(1)-exp(-1)) )
      }
    }
    
    #############################
    # EJEMPLO 4
    #############################
    
    if( ejemplo==4 ){
      
      corr.unif <- function(rho, n) {
        sp <- 2 * sin( rho * pi / 6)
        si <- chol(matrix(c(1, sp, sp, 1), 2, 2))
        return( pnorm( matrix( rnorm(2*n), n, 2) %*% si ) )
      }
      
      tmp <- corr.unif(rho= 0.7, n=n)
      X1 <- as.vector(tmp[,1])
      T1 <- as.vector(tmp[,2])
      
      X2 <- runif(n, 0, 1)
      T2 <- runif(n, 0, 1)
      
      #Generamos las covariables
      covariables <- cbind(X1,X2,T1,T2)
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    #############################
    # EJEMPLO 5  # Modelo 3 en el paper
    #############################
    
    if( ejemplo==5 ){
      
      #Generamos las covariables
      XT <- matrix(runif(n*2), n, 2)
      T1 <- XT[,1]
      T2 <- XT[,2]
      u1 <- rnorm(n, mean=0, sd= 0.1)
      u2 <- rnorm(n, mean=0, sd= 0.1)
      X1 <- XT[,1]+(XT[,2])^2 +u1 #XT[,1]*(XT[,2])^2 +u1
      X2 <- 1/2*( exp(XT[,1])-1) +u2 #1/2*( exp(XT[,1]*XT[,2])-1) +u2
      
      #Generamos las covariables
      covariables <- cbind(X1,X2,T1,T2)
      X <- covariables[,1:2]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    
    #############################
    # EJEMPLO 6
    #############################
    
    if( ejemplo==6 ){
      
      
      covariables <- matrix( rnorm(n*4, 1/2, sd=1/2), n, 4)
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1))
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-exp(7/12) )
      }
    }
    
    #############################
    # EJEMPLO 7
    #############################
    
    if( ejemplo==7 ){
      
      #Generamos las covariables
      covariables <- matrix( runif(n*4), n, 4)
      X <- covariables[,1:2]
      XT <- covariables[,3:4]
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    
    
    #############################
    # EJEMPLO 8
    #############################
    
    if( ejemplo==8 ){
      
      #Generamos las covariables
      X1 <- (rbinom(n,3,1/2))/3
      X2 <- (rbinom(n,5,1/5))/5
      X <- cbind(X1,X2)
      XT <- matrix( runif(n*2), n, 2)
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    
    
    #############################
    # EJEMPLO 9
    #############################
    
    if( ejemplo==9 ){
      
      #Generamos las covariables
      Xaux <- t( (rmultinom(n,10,prob=c(1/4,1/2,1/4)))/10 )
      X <- Xaux[,1:2]
      XT <- matrix( runif(n*2), n, 2)
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    
    
    #############################
    # EJEMPLO 10
    #############################
    
    if( ejemplo==10 ){
      
      #Generamos las covariables
      XT <- matrix( runif(n*2), n, 2)
      
      X1 <- rbinom(n,5,1/4)/5
      W <- rbinom(n,1,1/2)
      X2 <- 1/2*( (XT[,1]<(2/3)) + W)
      
      X <- cbind(X1,X2) 
      
      #Generamos los errores
      eps<-rnorm(n,0,desvio.epsilon)
      
      #Funciones
      funcion.g1 <- function(t1){
        return(2*sin(pi*t1)-4/pi)
      }
      funcion.g2 <- function(t2){
        return( exp(t2)-(exp(1)-1) )
      }
    }
    
    
    
    
    betas <- c(3,3)
    
    regresion <- as.vector(X%*%betas) + funcion.g1(XT[,1]) + funcion.g2(XT[,2])
    Y <- regresion + eps
    
    ##########################################################
    # FUNCION EME EVALUADA EN LA MUESTRA DE COVARIABLES
    ##########################################################
    
    verdaderos <- cbind( funcion.g1(XT[,1]), funcion.g2(XT[,2]))
    
    
    ##########################################################
    # si tipo.cont=0, no contamino
    ##########################################################
    
    if(tipo.cont==0){
      Y.cont <- Y
    }
    
    ##########################################################
    # si tipo.cont=1 CONTAMINACION ESPARCIDA-----
    ##########################################################
    
    if(tipo.cont==1){ 
      
      eps.cont <- 0.1
      
      ou <- rbinom(n, size=1, prob=eps.cont) 
      error.aux <- eps
      error.aux[ ou == 1 ] <- rnorm(sum(ou),mean=0, sd= (10*desvio.epsilon) )
      
      Y.cont <- regresion + error.aux
    }
    
    
    ##########################################################
    # si tipo.cont=2 CONTAMINACION ESPARCIDA + PUNTOS ARTIFICIALES
    ##########################################################
    
    
    if(tipo.cont==2){ 
      eps.cont <- 0.15
      
      ou <- rbinom(n, size=1, prob=eps.cont) 
      error.aux <- eps
      error.aux[ ou == 1 ] <- rnorm(sum(ou),mean=15, sd= 0.1 )
      
      Y.cont <- regresion + error.aux
    }
    
    if(tipo.cont==3){ 
      aux1 <- ((1:n)[(XT[,1]<1/3) & (XT[,2]<1/3)])[1]
      aux2 <- ((1:n)[(XT[,1]>1/3) & (XT[,1]<2/3) & (XT[,2]<1/3)])[1]
      aux3 <- ((1:n)[(XT[,1]>2/3) & (XT[,2]<1/3)])[1]
      aux4 <- ((1:n)[(XT[,1]<1/3) & (XT[,2]>1/3) & (XT[,2]<2/3)])[1]
      aux5 <- ((1:n)[(XT[,1]>1/3) & (XT[,1]<2/3) & (XT[,2]>1/3) & (XT[,2]<2/3)])[1]
      aux6 <- ((1:n)[(XT[,1]>2/3) & (XT[,2]>1/3) & (XT[,2]<2/3)])[1]
      aux7 <- ((1:n)[(XT[,1]<1/3) & (XT[,2]>2/3)])[1]
      aux8 <- ((1:n)[(XT[,1]>1/3) & (XT[,1]<2/3) & (XT[,2]>2/3)])[1]
      aux9 <- ((1:n)[(XT[,1]>2/3) & (XT[,2]>2/3)])[1]
      
      aux <- c(aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9)
      print(aux)
      Y.cont <- Y
      AUX <- try( X[aux,] <- matrix( rep(20,9*2), 9, 2) )
    }
    

  if(class(AUX)!= 'try-error'){ 
    
    ######################################################
    # ESTIMACIONES
    ######################################################
    
    grilla <- seq(0,1,length=1000) 
    ngrid <- length(grilla)
    grilla.np <- matrix(grilla,ngrid,2)
    sal <- plam.cl(Y.cont,X, XT, np.point=grilla.np)
    
    regresion.hat <- sal$prediction
    #g1.hat <- sal$g.matrix[,1]
    #g2.hat <- sal$g.matrix[,2]
    nknots <- sal$nknots
    nbasis <- sal$nbasis
    alpha <- sal$alpha
    
    g1.hat.grid <- sal$np.prediction[,1]
    g2.hat.grid <- sal$np.prediction[,2]
    
    #Con la rho.función
    residuos <- regresion.hat - regresion
    coef.rob[iter,] <- sal$coef.lin
    
    ecm.g1.rob[iter] <- mean( (g1.hat.grid - funcion.g1(grilla.np[,1]))^2 ,na.rm = T) #mean( (g1.hat - funcion.g1(XT[,1]))^2 ,na.rm = T)
    ecm.g2.rob[iter] <- mean( (g1.hat.grid - funcion.g1(grilla.np[,2]))^2 ,na.rm = T) #mean( (g2.hat - funcion.g2(XT[,2]))^2 ,na.rm=T)
    ecm.global.rob[iter] <- mean( ( residuos )^2, na.rm=T)

    ecm.g1.rob.t[iter] <- tms( g1.hat.grid - funcion.g1(grilla.np[,1]), alpha=0.05)
    ecm.g2.rob.t[iter] <- tms( g2.hat.grid - funcion.g2(grilla.np[,2]), alpha=0.05)
    ecm.global.rob.t[iter] <- tms( residuos, alpha=0.05)
    
    ecm.g1.rob.tt[iter] <- tms( g1.hat.grid - funcion.g1(grilla.np[,1]), alpha=0.1)
    ecm.g2.rob.tt[iter] <- tms( g2.hat.grid - funcion.g2(grilla.np[,2]), alpha=0.1)
    ecm.global.rob.tt[iter] <- tms( residuos, alpha=0.1)
    
    #Dentro del intervalo [0.05, 0.95]
    aux <- (grilla>=0.05) & (grilla<=0.95)
    ecm.g1.rob.int[iter] <- mean( (g1.hat.grid[aux] - funcion.g1(grilla.np[aux,1]))^2 ,na.rm = T) #mean( (g1.hat - funcion.g1(XT[,1]))^2 ,na.rm = T)
    ecm.g2.rob.int[iter] <- mean( (g2.hat.grid[aux] - funcion.g2(grilla.np[aux,2]))^2 ,na.rm=T)
    
    ecm.g1.rob.int.t[iter] <- tms( g1.hat.grid[aux] - funcion.g1(grilla.np[aux,1]), alpha=0.05)
    ecm.g2.rob.int.t[iter] <- tms( g2.hat.grid[aux] - funcion.g2(grilla.np[aux,2]), alpha=0.05)
    
    ecm.g1.rob.int.tt[iter] <- tms( g1.hat.grid[aux] - funcion.g1(grilla.np[aux,1]), alpha=0.1)
    ecm.g2.rob.int.tt[iter] <- tms( g2.hat.grid[aux] - funcion.g2(grilla.np[aux,2]), alpha=0.1)
    
    
    guardo <-c(iter,coef.rob[iter,],ecm.g1.rob[iter], ecm.g2.rob[iter], ecm.global.rob[iter],
               ecm.g1.rob.t[iter], ecm.g2.rob.t[iter], ecm.global.rob.t[iter], 
               ecm.g1.rob.tt[iter], ecm.g2.rob.tt[iter], ecm.global.rob.tt[iter],
               ecm.g1.rob.int[iter], ecm.g2.rob.int[iter],
               ecm.g1.rob.int.t[iter], ecm.g2.rob.int.t[iter],
               ecm.g1.rob.int.tt[iter], ecm.g2.rob.int.tt[iter],
               degree.spline, nknots, nbasis, alpha)
    
    #Guardo:
    largo <- length(guardo)
    write(t(guardo) , file=nombreRob, ncol=largo, append=T)
    
    #Guardo g1 en la grilla
    guardo.g1 <-c(iter, g1.hat.grid)
    largo1 <- length(guardo.g1)
    write(t(guardo.g1) , file=nombre.g1, ncol=largo1, append=T)
    
    #Guardo g2 en la grilla
    guardo.g2 <-c(iter, as.vector(g2.hat.grid))
    largo2 <- length(guardo.g2)
    write(t(guardo.g2) , file=nombre.g2, ncol=largo2, append=T)
    
  }#Fin del 'try-error'
  }
}

  



