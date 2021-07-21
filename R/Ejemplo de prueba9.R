#Boston Housing

datos <- read.csv("R/housing.csv", header = FALSE, sep="")
str(datos)
source("R/rplam-fn.R")

#Ejemplo PLAM de Ma Yang

x1 <- datos$V1
x2 <- datos$V3
x3 <- datos$V5
x4 <- datos$V6
x5 <- datos$V7
x6 <- datos$V8
x7 <- datos$V10
x8 <- datos$V11
x9 <- datos$V12
x10 <- datos$V13
X <- log(cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
colnames(X) <- NULL
y <- datos$V14

#Como el paper
y <- datos$V14

Z <- as.matrix(datos$V11)

X <- cbind(datos$V6, log(datos$V10), log(datos$V13))

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

library(robustbase)
set.seed(11)
nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

set.seed(111)
fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)
fit.rob$nknots

col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")


n <- length(y)
DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-1]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-2]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-3]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n),rep("Robust",n))
))
names(DF2) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1", "gx2","gx3","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )


res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:n,
  y-fit.rob$prediction
))
names(DF3) <- c("Position","Residuals")

p1 <- ggplot(DF3,aes(y=Residuals))+
  geom_boxplot(fill='#0052bb',alpha=0.5)+
  theme(
    axis.text.x=element_blank()
  )

p2 <- ggplot(DF3,aes(x=Position,y=Residuals))+
  geom_point(alpha=0.8,color="#0052bb")+
  geom_abline(intercept = 0,slope=0)+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank()
  )

p3 <- ggplot(DF3,aes(sample=Residuals))+
  geom_qq(color="#0052bb",alpha=0.8)+stat_qq_line()+
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  )

plot_grid(p1, p2, p3)

p1.info <- ggplot_build(p1)$data[[1]]

out <- p1.info$outliers[[1]]
out
length(out)
length(out)/n
outliers <- rep(0,length(out))
j <- 1
for(i in 1:n){
  if(sum(out==DF3$Residuals[i])>0){
    outliers[j] <- i
    j <- j+1
  }
}

DF3 <- as.data.frame(list(
  X,
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3]),
  rbind(fit.rob$g.matrix),
  rep("Robust",n)
))
names(DF3) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1","gx2","gx3","Fit")

DF3point <- as.data.frame(list(
  X[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3]))[outliers,]
)
)
names(DF3point) <- c("x1.x","x2.x","x3.x","re.x1","re.x2","re.x3")

ggplot(DF3,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x1.x,y=re.x1),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


out.pos <- outliers
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- as.matrix(Z[-out.pos,])

fit.del.cl <- plam.cl(y.del, Z.del, X.del)#, nknots=2)
fit.del.cl$nknots

col4 <- c(fit.del.cl$alpha,fit.del.cl$coef.lin)
tabla.b1.aux <- cbind(tabla.b1,col4)
row.names(tabla.b1.aux) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1.aux) <- c("Classical","Robust","Classical on clean data")
kable(tabla.b1.aux,caption = "Table 2")

DF4 <- as.data.frame(list(
  rbind(X.del,X),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - rowSums(fit.del.cl$g.matrix[,-1]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1])),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - rowSums(fit.del.cl$g.matrix[,-2]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2])),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - rowSums(fit.del.cl$g.matrix[,-3]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n-length(out.pos)),rep("Robust",n))
))
names(DF4) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1","gx2","gx3","Fit")

ggplot(DF4,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x1.x,y=re.x1),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


######################################
##Intervalos asumiendo independencia##
######################################


n <- length(y)
q <- dim(Z)[2]
library(plm)
for(i in 1:q){
  for(j in 1:3){
    detect.lindep(cbind(Z[,i],X[,j]))
  }
}

for(j in 1:3){
  detect.lindep(cbind(as.numeric(Z),X[,j]))
}

###Cálculo de la varianza asintótica suponiendo independencia
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma
plot(residuos)

coef.psi <- mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2

hstar <- as.matrix(colMeans(Z))
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (t(Z)[,i]-hstar)%*%t(t(Z)[,i]-hstar)
}
AA <- AA/n
det(AA)

SigmaDA <- sigma^2*coef.psi*solve(AA)
SigmaDA

#Ahora calculo los IC de Wald
zalpha <- qnorm(0.975,0,1)
intervalos <- matrix(0,q,2)
beta.hat <- fit.rob$coef.lin
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos

zalpha <- qnorm(0.90,0,1)
intervalos <- matrix(0,q,2)
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos

plot(Z,X[,1])
plot(Z,X[,2])
plot(Z,X[,3])

##-- Aproximación aditiva --##

##Si no entendí mal, A=E(Z-h*(X))(Z-h*(X))^t=E(Z-h*(X))^2 y esto se aproxima con el
#promedio de los residuos al cuadrado.
library(mgcv)
#help("gam")
fit.z <- gam(as.vector(Z)~s(X[,1])+s(X[,2])+s(X[,3]))
#plot(fit.z)
mean(fit.z$residuals^2)

#Ahora con la función que estima igual que plam pero para modelos aditivos.
degree.spline <- 1 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
nk.cl.am$nknots
fit.z <- am.cl(Z,X,degree.spline=degree.spline)
AA <- mean((Z-fit.z$prediction)^2)
AA

#Otra forma de hacer el cálculo
hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n
det(AA)

SigmaDA <- sigma^2*coef.psi*solve(AA)
SigmaDA

#Ahora calculo los IC de Wald
zalpha <- qnorm(0.975,0,1)
intervalos <- matrix(0,q,2)
beta.hat <- fit.rob$coef.lin
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos

zalpha <- qnorm(0.90,0,1)
intervalos <- matrix(0,q,2)
for(j in 1:q){
  lim_inf <- beta.hat[j]-zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  lim_sup <- beta.hat[j]+zalpha*sqrt(SigmaDA[j,j])/sqrt(n)
  intervalos[j,] <- c(lim_inf,lim_sup)
}
intervalos

###
## Usando la A y la D en el paso previo para hallar el Sigma de la DA
###

##--- Este es el que uso en el paper ---##

## Estimador robusto ##
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma
#plot(residuos)

degree.spline <- 2 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
nk.cl.am$nknots
fit.z <- am.cl(Z,X,degree.spline=degree.spline)

hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))

Atheta <- matrix(0,q,q)
for(i in 1:n){
  Atheta <- Atheta + psi.tukey.derivative(residuos[i])*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Atheta <- Atheta/n
Atheta

Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + (psi.tukey(residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

SigmaDA <- sigma^2*solve(Atheta)%*%Dtheta%*%t(solve(Atheta))
SigmaDA

#Desvío para el paper
sqrt(SigmaDA)/sqrt(n)


## Estimador clásico FULL ##
sigma <- sd(y-fit.full.cl$prediction)
residuos <- (y-fit.full.cl$prediction)/sigma
#plot(residuos)

degree.spline <- 2 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
nk.cl.am$nknots
fit.z <- am.cl(Z,X,degree.spline=degree.spline)

hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))

Atheta <- matrix(0,q,q)
for(i in 1:n){
  Atheta <- Atheta + 1*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Atheta <- Atheta/n
Atheta

Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + ((residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

SigmaDA <- sigma^2*solve(Atheta)%*%Dtheta%*%t(solve(Atheta))
SigmaDA

#Desvío para el paper
sqrt(SigmaDA)/sqrt(n)


## Estimador clásico CLEAN ##

n.del <- length(y.del)
q <- dim(Z.del)[2]

sigma.del <- sd(y.del-fit.del.cl$prediction)
residuos.del <- (y.del-fit.del.cl$prediction)/sigma.del
#plot(residuos)

degree.spline <- 2 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z.del,X.del,degree.spline=degree.spline)
nk.cl.am$nknots
fit.z.del <- am.cl(Z.del,X.del,degree.spline=degree.spline)

hstar.del <- as.matrix(fit.z.del$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))

Atheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Atheta.del <- Atheta.del + 1*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Atheta.del <- Atheta.del/n
Atheta.del

Dtheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Dtheta.del <- Dtheta.del + ((residuos.del[i]))^2*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Dtheta.del <- Dtheta.del/n
Dtheta.del

SigmaDA.del <- sigma.del^2*solve(Atheta.del)%*%Dtheta.del%*%t(solve(Atheta.del))
SigmaDA.del

#Desvío para el paper
sqrt(SigmaDA.del)/sqrt(n.del)






###Prueba estimando h* de forma robusta

degree.spline <- 1 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
nk.cl.am$nknots
fit.z <- am.rob(Z,X,degree.spline=degree.spline)
AA <- mean((Z-fit.z$prediction)^2)
AA

#Otra forma de hacer el cálculo
hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n
det(AA)

SigmaDA <- sigma^2*coef.psi*solve(AA)
SigmaDA

#Atheta y Dtheta con h* robusto
degree.spline <- 1 #Si se elige 3 da aún más grande

set.seed(1111)
fit.z <- am.rob(Z,X,degree.spline=degree.spline)
fit.z$nknots

hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$alpha+rowSums(fit.z$g.matrix))

Atheta <- matrix(0,q,q)
for(i in 1:n){
  Atheta <- Atheta + psi.tukey.derivative(residuos[i])*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Atheta <- Atheta/n
Atheta

Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + (psi.tukey(residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

SigmaDA <- sigma^2*solve(Atheta)%*%Dtheta%*%t(solve(Atheta))
SigmaDA



