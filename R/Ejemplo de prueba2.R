#Son li Z y X?

install.packages("plm")
library(plm)

library(devtools)
install_github("alemermartinez/rplam")
library(rplam)

library(robustbase)


data(airquality)
x <- airquality
x <- x[ complete.cases(x), ]
x <- x[, c('Ozone', 'Solar.R', 'Wind', 'Temp','Month')]
y <- as.vector(x$Ozone)
X <- as.matrix(x[, c('Solar.R', 'Wind', 'Temp')])
Z <- as.matrix(x[, 'Month'])

Z <- as.factor(Z)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)

beta.hat <- fit.rob$coef.lin



n <- length(y)
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

for(i in 1:q){
  for(j in 1:3){
    detect.lindep(cbind(Z.aux[,i],X[,j]))
  }
}

for(j in 1:3){
  detect.lindep(cbind(as.numeric(Z),X[,j]))
}

x1 <- (X[,3])^2
x2 <- X[,3]
sal <- lm(as.numeric(Z)~ x1+x2)
plot(sal$residuals)
sal$coefficients
detect.lindep(cbind(as.numeric(Z), -19.34-0.00305*x1+0.5312*x2))

plot(cbind(as.numeric(Z)-( -19.34-0.00305*x1+0.5312*x2)))


###Cálculo de la varianza asintótica
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma
plot(residuos)

mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2

mean(psi.tukey(residuos)^2, trim=0.05)/(mean(psi.tukey.derivative(residuos), trim=0.05))^2


#Cuando reemplazo el promedio por el estimador de Tukey me da 0
#maxit=100
#control <- lmrob.control(trace.level = 0,         # 0
#                         nResample   =  500,      # 500 default
#                         tuning.psi = 4.685061,      # para 85% eff usar 3.443689 # para 95% eff usar 4.685061
#                         subsampling = 'simple',  #
#                         rel.tol     = 1e-5,      # 1e-7
#                         refine.tol  = 1e-5,      # 1e-7
#                         k.max       = 2e3,       # 200
#                         maxit.scale = maxit,       # 200 #2e3
#                         max.it      = maxit)       # 50 #2e3
#sal1  <- lmrob(psi.tukey(residuos)^2 ~ 1, control = control)
#sal1$coefficients
#sal2  <- lmrob(psi.tukey.derivative(residuos) ~ 1, control = control)
#sal2$coefficients

sigma1 <- mad(psi.tukey(residuos)^2)
num <- pos.est(psi.tukey(residuos)^2, sigma1, typePhi="Tukey", ini=NULL, epsilon=1e-6, iter.max=50)
sigma2 <- mad(psi.tukey.derivative(residuos))
den <- (pos.est(psi.tukey.derivative(residuos), sigma2, typePhi="Tukey", ini=NULL, epsilon=1e-6, iter.max=50))^2
coef.psi <- num/den
coef.psi

num <- pos.est(psi.tukey(residuos)^2, sigma1, typePhi="Huber", ini=NULL, epsilon=1e-6, iter.max=50)
den <- (pos.est(psi.tukey.derivative(residuos), sigma2, typePhi="Tukey", ini=NULL, epsilon=1e-6, iter.max=50))^2
num/den

#Chequeo bajo normalidad
xx <- rnorm(100000,0,1)
coef.psi <- mean((psi.tukey(xx))^2)/(mean(psi.tukey.derivative(xx)))^2

hstar <- as.matrix(colMeans(Z.aux))
q <- 4
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (t(Z.aux)[,i]-hstar)%*%t(t(Z.aux)[,i]-hstar)
}
AA <- AA/n
det(AA)

SigmaDA <- sigma^2*coef.psi*solve(AA)
SigmaDA

#Ahora calculo los IC de Wald
zalpha <- qnorm(0.975,0,1)
intervalos <- matrix(0,q,2)
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos

zalpha <- qnorm(0.95,0,1)
intervalos <- matrix(0,q,2)
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos


