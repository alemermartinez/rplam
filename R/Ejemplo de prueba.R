library(robustbase)
library(splines)

####Ejemplos de prueba

#- Sin ceros -#

set.seed(123)

#Cantidad de muestras
n <- 100

#beta
beta <- as.matrix(c(3,2,1))

#funciones aditivas
function.g1 <- function(x1) 24*(x1-1/2)^2-(6/5)
function.g2 <- function(x2) 2*pi*sin(pi*x2)-(48/(pi^2))

z <- matrix(rnorm(n*3,0,sd=0.5),n,3)
x <- matrix(runif(n*2,0,1),n,2)

eps <- rnorm(n,0,sd=0.5)

y <- z%*%beta + function.g1(x[,1]) + function.g2(x[,2]) + eps

#Prueba de gráficos
plot(x[,1], y-z%*%beta-function.g2(x[,2]))
plot(x[,2], y-z%*%beta-function.g1(x[,1]))

#Estimación sin selección de variables
#Estimación clásica:
sal.cl <- plam.cl(y, z, x)
sal.cl$alpha
sal.cl$coef.lin
sal.cl$g.matrix
sal.cl$nknots
sal.cl$kj

#Estimación robusta:
sal.rob <- plam.rob(y, z, x)
sal.rob$alpha
sal.rob$coef.lin
sal.rob$g.matrix
sal.rob$sigma.hat
sal.rob$nknots
sal.rob$kj


#- Con ceros -#

set.seed(123)

#Cantidad de muestras
n <- 100

#beta
beta <- as.matrix(c(3,2,1,0,0,0))

#funciones aditivas
function.g1 <- function(x1) 24*(x1-1/2)^2-(6/5)
function.g2 <- function(x2) 2*pi*sin(pi*x2)-(48/(pi^2))

z <- matrix(rnorm(n*6,0,sd=0.5),n,6)
x <- matrix(runif(n*4,0,1),n,4)

eps <- rnorm(n,0,sd=0.5)

y <- z%*%beta + function.g1(x[,1]) + function.g2(x[,2]) + eps

#Prueba de gráficos
plot(x[,1], y-z%*%beta-function.g2(x[,2]))
plot(x[,2], y-z%*%beta-function.g1(x[,1]))

#Estimación sin selección de variables
#Estimación clásica:
sal.cl <- plam.cl(y, z, x)
sal.cl$alpha
sal.cl$coef.lin
sal.cl$g.matrix
sal.cl$nknots
sal.cl$kj

#Estimación robusta:
sal.rob <- plam.rob(y, z, x)
sal.rob$alpha
sal.rob$coef.lin
sal.rob$g.matrix
sal.rob$sigma.hat
sal.rob$nknots
sal.rob$kj

### Selección de variables
#Con lambdas fijos

#Estimación clásica
sal.cl.1 <- plam.cl.vs.lambdas(y, z, x, lambdas1=rep(0.2,7), lambdas2=rep(0.3,4)) #El lambda1 es q+1 porque está mirando el alpha
sal.cl.1$coef.const
sal.cl.1$coef.lin
sal.cl.1$is.zero

#Estimación robusta
sal.rob.1 <- plam.rob.vs.lambdas(y, z, x, lambdas1=rep(0.2,7), lambdas2=rep(0.3,4))
sal.rob.1$coef.const
sal.rob.1$coef.lin
sal.rob.1$is.zero

### Selección automática de lambdas

#Estimación clásica
lambdas.cl <- select.cl.lambdas(y, z, x, grid.lambda1=c(0.1,0.2,0.3,0.4), grid.lambda2=c(0.1,0.2,0.3,0.4))
lambdas.cl

#Estimación robusta
lambdas.rob <- select.rob.lambdas(y, z, x, grid.lambda1=c(0.1,0.2,0.3,0.4), grid.lambda2=c(0.1,0.2,0.3,0.4))
lambdas.rob

### Selección de variables con selección automática de lambdas

#Estimación clásica
sal.cl.2 <- plam.rob.vs(y, z, x)
sal.cl.2$coef.const
sal.cl.2$coef.lin
sal.cl.2$nknots
sal.cl.2$kj
sal.cl.2$lambda1
sal.cl.2$lambda2
sal.cl.2$is.zero

#Estimación robusta
sal.rob.2 <- plam.rob.vs(y, z, x)
sal.rob.2$sigma.hat
sal.rob.2$coef.const
sal.rob.2$coef.lin
sal.rob.2$nknots
sal.rob.2$kj
sal.rob.2$lambda1
sal.rob.2$lambda2
sal.rob.2$is.zero

