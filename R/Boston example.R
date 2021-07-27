######################
#-- Boston Housing --#
######################

rm(list=ls())
# Cargamos los datos
library(MASS)
attach(Boston)
str(Boston)

# Cargamos el archivo con las funciones
source("rplam-fn.R")

# Seleccionamos las variables tal como aparecen en el ejemplo
# de Ma y Yang (2011)

y <- medv #MEDV

Z <- as.matrix(ptratio) #PTRATIO

zz=ptratio
names(zz)=1:506

outZ=unlist(names(boxplot(zz)$out))

X <- cbind(rm, log(tax), log(lstat))
# RM, TAX y LSTAT

# Usamos splines cúbicos
degree.spline <- 3

# Calculamos el estimador clásico
fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

# Si queremos ver la cantidad de nodos internos seleccionados
fit.full.cl$nknots

# Lo que es siguiente a la siguiente cantidad de sumandos
fit.full.cl$kj
# En el paper esto sería el kj-1 (y coincide con nknots + degree.spline)

# Ahora calculamos el estimador robusto
set.seed(111)
fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)

# Fijé la semilla para que de siempre lo mismo. Con otras
# semillas también da la misma cantidad de nodos internos que es

fit.rob$nknots

# Ahora construimos la tabla para ver las estimaciones
# de la constante mu y del beta
col2 <- c(fit.full.cl$coef.const,fit.full.cl$coef.lin)
col3 <- c(fit.rob$coef.const,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

# Lo que sigue es para graficar las curvas
# estimadas usando ggplot2
n <- length(y)
DF2 <- as.data.frame(list(
  rbind(X,X),
  c(rep(NA,n),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1])),
  c(rep(NA,n),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2])),
  c(rep(NA,n),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n),rep("Robust",n))
))
names(DF2) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1", "gx2","gx3","Fit")
library(ggplot2)

DF2.rojo <- as.data.frame(list(
  rbind(X),
  rbind(fit.full.cl$g.matrix),
  c(rep("Classical",n))
))
names(DF2.rojo) <- c("x1","x2","x3","gx1","gx2","gx3","Fit")

color.fit <- c("red","blue") #c("#bb0c00","blue") #c("#bb0c00","#0052bb") #c("red","blue") #c("#bb0c00","#0052bb")

# En los siguientes gráficos se pueden ver:
# En rojo el estimador clásico y sus residuos parciales
# En azul el estimador robusto y sus residuos parciales

# Curvas g1
ggplot(DF2,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
)

# The following plot is for the paper:
nombre="ggplot-boston-1-bis-big.pdf"  #"ggplot-1-bis.pdf"
pdf(nombre, bg="transparent")

ggplot(DF2,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.3, shape=16, size=6)+ #size=2 alpha=0.5
  geom_line(aes(x=x1,y=gx1),size=4,alpha=1,linetype=c(rep(2,n),rep(1,n)))+ #c(rep(2,111),rep(1,111)), size=1.8
  scale_color_manual(values = color.fit)+
  geom_line(data=DF2.rojo,aes(x=x1,y=gx1),size=4,alpha=1,linetype=c(rep(2,n))) #size=1.8
dev.off()


# Curvas g2
ggplot(DF2,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

# The following plot is for the paper:
nombre="ggplot-boston-2-bis-big.pdf"  #"ggplot-1-bis.pdf"
pdf(nombre, bg="transparent")

ggplot(DF2,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.3, shape=16, size=6)+ #size=2 alpha=0.5
  geom_line(aes(x=x2,y=gx2),size=4,alpha=1,linetype=c(rep(2,n),rep(1,n)))+ #c(rep(2,111),rep(1,111)), size=1.8
  scale_color_manual(values = color.fit)+
  geom_line(data=DF2.rojo,aes(x=x2,y=gx2),size=4,alpha=1,linetype=c(rep(2,n))) #size=1.8
dev.off()


# Curvas g3
ggplot(DF2,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

# The following plot is for the paper:
nombre="ggplot-boston-3-bis-big.pdf"  #"ggplot-1-bis.pdf"
pdf(nombre, bg="transparent")

ggplot(DF2,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.3, shape=16, size=6)+ #size=2 alpha=0.5
  geom_line(aes(x=x3,y=gx3),size=4,alpha=1,linetype=c(rep(2,n),rep(1,n)))+ #c(rep(2,111),rep(1,111)), size=1.8
  scale_color_manual(values = color.fit)+
  geom_line(data=DF2.rojo,aes(x=x3,y=gx3),size=4,alpha=1,linetype=c(rep(2,n))) #size=1.8
dev.off()


# Las mayores diferencias se pueden apreciar en
# las variables 1 y 3.

# Ahora miramos los residuos obtenidos por el
# estimador robusto
res.rob <- y-fit.rob$prediction
names(res.rob) <- 1:n

summary(res.rob)


# Haremos tres gráficos usando ggplot
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

p1

# The following plot is for the paper:
nombre="boxplot-boston-res-bis-big.pdf"
pdf(nombre, bg="transparent")

p1 <- ggplot(DF3,aes(y=Residuals))+
  stat_boxplot(geom = "errorbar", width = 0.5)+
  geom_boxplot(fill='gray70',alpha=1, outlier.size=3)+ #2.5
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
p1
dev.off()

pZ <- ggplot(DF3,aes(y=Z))+
  geom_boxplot(fill='#0052bb',alpha=0.5)+
  theme(
    axis.text.x=element_blank()
  )

pZ

# The following plot is for the paper:
nombre="boxplot-boston-z.pdf"
pdf(nombre, bg="transparent")

pz <- ggplot(DF3,aes(y=Z))+
  stat_boxplot(geom = "errorbar", width = 0.5)+
  geom_boxplot(fill='gray70',alpha=1, outlier.size=3)+ #2.5
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
pz
dev.off()

# Ahora buscamos qué datos fueron detectados
# como datos atípicos por el estimador robusto.

outliers <- as.numeric(names(boxplot(res.rob, plot=FALSE)$out))

outliers
outZ

sort(boxplot(res.rob, plot=FALSE)$out)

plot(Z,res.rob)
abline(h=res.rob[343],col="red") #mayor outlier negativo
abline(h=res.rob[409],col="red") #menor outlier positivo

abline(v=13,col="blue", lty=2) #mayor outlier de Z (SON SIEMPRE CHICOS)


DFz <- as.data.frame(list(
  Z,
  y-fit.rob$prediction
))
names(DFz) <- c("Z","Residuals")

ggplot(DFz,aes(x=Z,y=Residuals)) +
  geom_point(alpha=0.3, shape=16, size=6)+
  geom_hline(yintercept=-7.987732,col="red")+
  geom_hline(yintercept=7.578529,col="red")+
  geom_vline(xintercept=13,col="blue")+
  theme(
    axis.title.y=element_blank()
  )


plot(Z,res.rob)
abline(h=res.rob[343],col="red") #mayor outlier negativo
abline(h=res.rob[409],col="red") #menor outlier positivo

abline(v=13,col="blue", lty=2)




# Los residuos detectados como atípicos son:
outliers
# Dando un total de:
length(outliers)
# datos atípicos que representan un
length(outliers)/n*100
# porciento de los datos.


# Y ahora realizamos los gráficos.
# Sobre las curvas robustas y los correspondientes
# residuos parciales, graficamos en magenta los
# residuos correspondientes a estos datos atípicos
DF3 <- as.data.frame(list(
  X,
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3]),
  rbind(fit.rob$g.matrix),
  rep("Robust",n)
))
names(DF3) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1","gx2","gx3","Fit")

DF3point <- as.data.frame(list(
  X[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3]))[outliers,]
)
)
names(DF3point) <- c("x1.x","x2.x","x3.x","re.x1","re.x2","re.x3")

# Curva g1, residuos parciales y residuos parciales
# de los datos atípicos
ggplot(DF3,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x1.x,y=re.x1),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

# Curva g2, residuos parciales y residuos parciales
# de los datos atípicos
ggplot(DF3,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

# Curva g3, residuos parciales y residuos parciales
# de los datos atípicos
ggplot(DF3,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )





# Lo que haremos ahora será recalcular el estimador
# clásico sin esos datos atípicos. Para ello, primero
# removemos los outliers:
out.pos <- outliers
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- as.matrix(Z[-out.pos,])

# Y ahora recalculamos el estimador clásico
fit.del.cl <- plam.cl(y.del, Z.del, X.del)

# Observamos que se siguen seleccionando la misma
# cantidad de nodos internos
fit.del.cl$nknots

# Extendemos la tabla con la estimación de los
# coeficientes con una nueva columna con mu y
# beta estimados por el estimador clásico con
# los datos atípicos removidos.
col4 <- c(fit.del.cl$coef.const,fit.del.cl$coef.lin)
tabla.b1.aux <- cbind(tabla.b1,col4)
row.names(tabla.b1.aux) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1.aux) <- c("Classical","Robust","Classical on clean data")
kable(tabla.b1.aux,caption = "Table 2")
# Se puede observar que ahora los coeficientes
# se encuentran más cerca de las estimaciones
# obtenidas con el estimador robusto utilizando la
# totalidad de los datos.

# Y, por último, graficamos las nuevas curvas obtenidas
# por el estimador clásico sin los outliers y las comparamos
# con las curvas obtenidas con el estimador robusto.
# Además de calcular en rojo las curvas obtenidas con el estimador
# clásico y los residuos parciales y en azul las curvas y los
# residuos parciales obtenidos con el estimador robusto, en magenta
# se pueden observar nuevamente los residuos parciales de los
# datos atípicos.
DF4 <- as.data.frame(list(
  rbind(X.del,X),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$coef.const - rowSums(fit.del.cl$g.matrix[,-1]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1])),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$coef.const - rowSums(fit.del.cl$g.matrix[,-2]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2])),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$coef.const - rowSums(fit.del.cl$g.matrix[,-3]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n-length(out.pos)),rep("Robust",n))
))
names(DF4) <- c("x1","x2","x3","re.x1", "re.x2","re.x3", "gx1","gx2","gx3","Fit")

# Curvas g1 y residuos parciales
ggplot(DF4,aes(x=x1,y=re.x1,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x1,y=gx1),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x1.x,y=re.x1),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

# Curvas g2 y residuos parciales
ggplot(DF4,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

# Curvas g3 y residuos parciales
ggplot(DF4,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


DF3point <- as.data.frame(list(
  X[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2]))[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3]))[outliers,]
)
)
names(DF3point) <- c("x1.x","x2.x","x3.x","re.x1","re.x2","re.x3")


DF4 <- as.data.frame(list(
  rbind(X.del,X),
  c(rep(NA,n-length(outliers)),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-1])),
  c(rep(NA,n-length(outliers)),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-2])),
  c(rep(NA,n-length(outliers)),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$coef.const - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n-length(outliers)),rep("Robust",n))
))
names(DF4) <- c("x1","x2","x3","re.x1", "re.x2", "re.x3", "gx1","gx2","gx3","Fit")

DF4.rojo <- as.data.frame(list(
  rbind(X.del),
  rbind(fit.del.cl$g.matrix),
  c(rep("Classical",n-length(outliers)))
))
names(DF4.rojo) <- c("x1","x2","x3","gx1","gx2","gx3","Fit")


# The following plot is for the paper:
nombre="ggplot-boston-final-1-bis-big.pdf"
pdf(nombre, bg="transparent")

ggplot(DF4,aes(x=x1,y=re.x1,color=Fit))+
  geom_point(alpha=0.3, shape=16, size=6)+
  geom_line(aes(x=x1,y=gx1),size=4,alpha=1,linetype=c(rep(2,n-length(outliers)),rep(1,n)))+
  scale_color_manual(values = color.fit)+
  geom_point(data=DF3point,aes(x=x1.x,y=re.x1),col='black',size=6)+
  geom_line(data=DF4.rojo,aes(x=x1,y=gx1),size=4,alpha=1,linetype=c(rep(2,n-length(outliers))))
dev.off()

# The following plot is for the paper:
nombre="ggplot-boston-final-2-bis-big.pdf"
pdf(nombre, bg="transparent")

ggplot(DF4,aes(x=x2,y=re.x2,color=Fit))+
  geom_point(alpha=0.3, shape=16, size=6)+
  geom_line(aes(x=x2,y=gx2),size=4,alpha=1,linetype=c(rep(2,n-length(outliers)),rep(1,n)))+
  scale_color_manual(values = color.fit)+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='black',size=6)+
  geom_line(data=DF4.rojo,aes(x=x2,y=gx2),size=4,alpha=1,linetype=c(rep(2,n-length(outliers))))
dev.off()

# The following plot is for the paper:
nombre="ggplot-boston-final-3-bis-big.pdf"
pdf(nombre, bg="transparent")

ggplot(DF4,aes(x=x3,y=re.x3,color=Fit))+
  geom_point(alpha=0.3, shape=16, size=6)+
  geom_line(aes(x=x3,y=gx3),size=4,alpha=1,linetype=c(rep(2,n-length(outliers)),rep(1,n)))+
  scale_color_manual(values = color.fit)+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='black',size=6)+
  geom_line(data=DF4.rojo,aes(x=x3,y=gx3),size=4,alpha=1,linetype=c(rep(2,n-length(outliers))))
dev.off()


# Se pueden observar que no sólo las curvas son ahora prácticamente
# idénticas sino que también son muy parecidos los residuos parciales
# de las observaciones que ambos estimadores comparten.



############################
## Cálculo de los desvíos ##
############################

# A continuación haremos tres tablas:
# La primera estimando la varianza asintótica Sigma
# como el producto de matrices
# La segunda estimando directamente la A
# La tercera calculando la estimación de A
# para el estimador robusto utilizando pesos

####### SPLINES CUADRÁTICOS

#--- PRIMERA ---#
# Es decir, usando la B y la D en el paso previo para hallar el Sigma de la DA

# Comenzamos estimando la h* con un estimador robusto:

# Para ello, definimos primer los residuos estandarizados:
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Usaremos splines cuadráticos
degree.spline <- 2

# Como queremos ver dónde ocurre el primer mínimo
# local, graficaremos el RBIC
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)

# Si bien la cantidad de nodos que se hubiese
# elegido es:
nk.rob.am$nknots

# El primer mínimo local lo tenemos en
nknots <- 8

# Luego, modelamos la relación entre Z y las covariables
# utilizando un modelo aditivo. Estimamos usando B-splines
# cuadráticos y un estimador robusto:
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# De esta manera, h* es:
hstar <- as.matrix(fit.z.rob$prediction)
# Que sería equivalente a: as.matrix(fit.z.rob$coef.const+rowSums(fit.z.rob$g.matrix))

# Ahora calculamos la matriz B. Observemos que es de
# qxq y, en este caso, q=1
q <- 1
Btheta <- matrix(0,q,q)
for(i in 1:n){
  Btheta <- Btheta + psi.tukey.derivative(residuos[i])*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Btheta <- Btheta/n
Btheta

# Ahora calculamos la matriz D (también de qxq)
Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + (psi.tukey(residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

# Luego, matriz de covarianzas asintótica es
SigmaDA <- sigma^2*solve(Btheta)%*%Dtheta%*%t(solve(Btheta))
SigmaDA

# Y, por lo tanto, el desvío para el paper es:
desvio.rob.BD <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.BD

#############################################################
# Ahora lo que haremos es calcular el desvío obtenido
# para el estimador clásico con la totalidad de los
# datos.
#####################################################

# Comenzamos nuevamente definiendo sus residuos estandarizados.
sigma <- sd(y-fit.full.cl$prediction)

residuos <- (y-fit.full.cl$prediction)/sigma

# Observemos que estimé el sigma con el sd (por si queremos
# considerar otro estimador)

# Ahora analizamos el BIC al utilizar splines cuadráticos
degree.spline <- 2
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
plot(nk.cl.am$grid, nk.cl.am$BIC)

# Si bien la cantidad de nodos que minimiza el BIC es:
nk.cl.am$nknots

# Seleccionaremos según el gráfico:
nknots <- 6

# Luego, estimamos el modelo aditivo que modela Z en función de
# las X's utilizando un estimador clásico y splines cuadráticos
fit.z.cl <- am.cl(Z,X,degree.spline=degree.spline, nknots=nknots)

# De esta manera, h* se define como:
hstar <- as.matrix(fit.z.cl$prediction)

# Ahora calculamos la matriz de qxq B
Btheta <- matrix(0,q,q)
for(i in 1:n){
  Btheta <- Btheta + 1*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Btheta <- Btheta/n
Btheta

# Y la matriz D
Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + ((residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

# De esta manera, la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*solve(Btheta)%*%Dtheta%*%t(solve(Btheta))
SigmaDA

# Y, por lo tanto, el desvío para el paper para este caso es:
desvio.cl.full.BD <- sqrt(SigmaDA)/sqrt(n)
desvio.cl.full.BD

##################################################
# Ahora, repetimos este último procedimiento pero para
# el estimador clásico que se obtiene cuando los outliers
# fueron removidos.

# Observemos que ahora tenemos menos datos, por lo que
# tenemos otros n
n.del <- length(y.del)

# Los nuevos residuos son:
sigma.del <- sd(y.del-fit.del.cl$prediction)
residuos.del <- (y.del-fit.del.cl$prediction)/sigma.del

# Analizamos el BIC usando splines cuadráticos
degree.spline <- 2
nk.cl.del.am <- select.nknots.cl.am(Z.del,X.del,degree.spline=degree.spline)
plot(nk.cl.del.am$BIC)

# Si bien la cantidad total de nodos internos hubiese sido
nk.cl.del.am$nknots
# Vamos a considerar:
nknots <- 6

# Ahora estimamos de forma clásica el modelo aditivo que
# modela Z.del en términos de X.del
fit.z.del <- am.cl(Z.del,X.del,degree.spline=degree.spline, nknots=nknots)

# Y obtenemos la siguiente h*
hstar.del <- as.matrix(fit.z.del$prediction)

# Ahora calculamos B
Btheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Btheta.del <- Btheta.del + 1*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Btheta.del <- Btheta.del/n.del
Btheta.del
# Obsérvese que en este caso dividimos por n.del y no por n.

# Y calculamos D
Dtheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Dtheta.del <- Dtheta.del + ((residuos.del[i]))^2*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Dtheta.del <- Dtheta.del/n.del
Dtheta.del

# De esta manera, la nueva matriz de covarianzas asintótica es:
SigmaDA.del <- sigma.del^2*solve(Btheta.del)%*%Dtheta.del%*%t(solve(Btheta.del))
SigmaDA.del

# Y, por lo tanto, el desvío para el paper es:
desvio.cl.clean.BD <- sqrt(SigmaDA.del)/sqrt(n.del)
desvio.cl.clean.BD

# La primera tabla es entonces
col2 <- c(desvio.cl.full.BD, desvio.rob.BD, desvio.cl.clean.BD)
tabla.desBD <- rbind(col2) #cbind(col2)
row.names(tabla.desBD) <- c('Desvios')
colnames(tabla.desBD) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.desBD,caption = "Table 1")


#--- SEGUNDA ---#
# Es decir, usando directamente la definición de A

# Comenzamos estimando la h* con un estimador robusto:

# Calculamos los residuos
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Analizamos el RBIC como antes (de hecho, da lo mismo)
# para seleccionar la cantidad de nodos
degree.spline <- 2
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)
nknots <- 8

# Nuevamente, estimamos el modelo aditivo que relaciona Z con X
# y realizamos la estimación de forma robusta
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# Luego, el h* es, nuevamente,
hstar <- as.matrix(fit.z.rob$prediction)

# Calculamos ahora el coeficiente "v" según la notación del paper
coef.v <- mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2
coef.v

# Calculamos ahora la matriz A:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n

# Luego, la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*coef.v*solve(AA)
SigmaDA

# Y, por lo tanto, el desvío para el paper es:
desvio.rob.A <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.A

#################################################################
# Ahora repetimos el procedimiento con esl estimador
# clásico que utiliza la totalidad de los datos: con la matriz sigma^2 A
#################################################################

# Calculamos primero los residuos y luego analizamos
# el BIC tal como lo hacíamos antes:
sigma <- sd(y-fit.full.cl$prediction)


degree.spline <- 2
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
plot(nk.cl.am$grid, nk.cl.am$BIC)
nknots <- 6

# Estimamos de forma clásica la relación entre Z y X
# suponiendo un modelo aditivo
fit.z <- am.cl(Z,X,degree.spline=degree.spline, nknots=nknots)

# Luego, nuevamente h* es:
hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$coef.const+rowSums(fit.z$g.matrix))

# El coeficiente "v" en este caso se estima como:
residuos <- (y-fit.full.cl$prediction)/sigma

coef.v <- mean((residuos)^2)  #deberia dar cerca de 1
coef.v
sigma ^2

# Ahora calculamos la matriz A:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n

# La matriz de covarianzas asintótica estimada es:
SigmaDA <- sigma^2*solve(AA)
SigmaDA

# Finalmente, el desvío para este caso es:
desvio.cl.full.A <- sqrt(SigmaDA)/sqrt(n)
desvio.cl.full.A

# Luego, repetimos el análisis con el estimador
# clásico con los datos "limpios" (sin outliers)

# Calculamos los residuos
n.del <- length(y.del)
sigma.del <- sd(y.del-fit.del.cl$prediction)
residuos.del <- (y.del-fit.del.cl$prediction)/sigma.del

# Nuevamente, miramos el BIC para hallar la cantidad de
# nodos
degree.spline <- 2 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z.del,X.del,degree.spline=degree.spline)
plot(nk.cl.am$BIC)
nknots <- 6

# Realizamos el ajuste del modelo aditivo
fit.z.del <- am.cl(Z.del,X.del,degree.spline=degree.spline, nknots=nknots)

# Y entonces obtenemos el siguiente estimador de h*:
hstar.del <- as.matrix(fit.z.del$prediction)

# El coeficiente "v" es, en este caso:
coef.v.del <- mean((residuos.del)^2)
coef.v.del      #DEBERIA SER CERCANO A 1

sigma.del^2

# Calculamos la matriz A estimada:
q <- 1
AA.del <- matrix(0,q,q)
for(i in 1:n.del){
  AA.del <- AA.del + (Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
AA.del <- AA.del/n.del

# Luego, el estimador de la matriz de covarianzas asintótica es:
SigmaDA.del <- sigma.del^2*solve(AA.del)
SigmaDA.del

# Y, por lo tanto, el desvío para el paper es:
desvio.cl.clean.A <- sqrt(SigmaDA.del)/sqrt(n.del)
desvio.cl.clean.A


# La segunda tabla es entonces
col2 <- c(desvio.cl.full.A, desvio.rob.A, desvio.cl.clean.A)
tabla.des <- rbind(col2) #cbind(col2)
row.names(tabla.des) <- c('Desvios')
colnames(tabla.des) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.des,caption = "Table 1")


#--- TERCERA ---#
# Esta tabla es similar a la anterior pero la matriz
# A para el estimador robusto se estima con pesos:

# Calculamos nuevamente los residuos:
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Calculamos los pesos que utilizan estos residuos
pesos <- psi.w(residuos)

# Repetimos el procedimiento que ya hemos
# realizado en dos oportunidades. Primero seleccionamos
# la cantidad de nodos internos mirando el RBIC y luego
# estimamos:
degree.spline <- 2
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)
nknots <- 8
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# Obtenemos el siguiente estimador de h*:
hstar <- as.matrix(fit.z.rob$prediction)

# Volvemos a calcular el estimador de "v" en este caso:
coef.v <- mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2
coef.v

# Ahora calculamos el estimador de A usando pesos:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + pesos[i]*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/(sum(pesos))

# De esta manera, el estimador de la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*coef.v*solve(AA)
SigmaDA

# Y, por lo tanto, el desvío para este caso es:
desvio.rob.A.pesos <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.A.pesos

# Luego, la tercera tabla quedaría:
col2 <- c(desvio.cl.full.A, desvio.rob.A.pesos, desvio.cl.clean.A)
tabla.des <- rbind(col2) #cbind(col2)
row.names(tabla.des) <- c('Desvios')
colnames(tabla.des) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.des,caption = "Table 2,  estimo usando Sigma=sigma^2 v A^{-1}")

kable(tabla.desBD,caption = "Table 1, estimo usando Sigma=B^{-1} D B^{-1}")

# LO QUE SE DEBE REPORTAR EN EL PAPER

pirulo=rbind(tabla.b1.aux[2,], c(tabla.des [1],   tabla.desBD[2], tabla.des[3]))
tabla.paper.beta= round( pirulo, 4)
 row.names(tabla.paper.beta)=c(  "$\\widehat{\\beta}$", "\\widehat{s}_{\\widehat{\\beta}}")
colnames(tabla.paper.beta) <- c("Classical","Robust","Classical on clean data")
kable(tabla.paper.beta,caption = "Estimates and SD")

alpha <- 0.95
zalpha <- qnorm(1-(1-alpha)/2,0,1)

tabla.paper.beta[1,]-zalpha*tabla.paper.beta[2,]
tabla.paper.beta[1,]+zalpha*tabla.paper.beta[2,]

tabla.paper.beta[1,]-3*tabla.paper.beta[2,]

#             Classical                  Robust Classical on clean data
#                -0.9095                 -0.5543                 -0.7286

# EL EXTREMO del IC PARA EL ROBUSTO NO INCLUYE AL CLASICO!!

tabla.paper.beta[1,]+3*tabla.paper.beta[2,]
#            Classical                  Robust Classical on clean data
#                -0.2405                 -0.4655                 -0.2180

tabla.paper.beta[1,]+3*0.007518721


#################################
# SI USO SIGMA= sigma^2*v*A^{-1}
###################################

round(tabla.paper.beta[1,]-3*tabla.des,4)

#         Classical     Robust Classical on clean data
#        -0.9095 -0.5783                 -0.7287

# ACA CON EL DESVIO USANDO A  INCLUYE AL CLASICO! porque los devios son 50% mas grandes

round(tabla.paper.beta[1,]+3*tabla.des,4)
 #      Classical  Robust Classical on clean data
 #        -0.2405 -0.4415                 -0.2179


tabla.des/tabla.desBD
#       Classical   Robust Classical on clean data
#COCIENTES  1.181237 1.539971                1.220519


##################################################
# RESULTADOS
###################################################
kable(tabla.paper.beta,caption = "Estimates and SD")


Table: Estimates and SD

|                              | Classical|  Robust| Classical on clean data|
|:-----------------------------|---------:|-------:|-----------------------:|
|$\widehat{\beta}$             |   -0.5750| -0.5099|                 -0.4733|
|\widehat{s}_{\widehat{\beta}} |    0.1115|  0.0148|                  0.0851|


library(xtable)


print(xtable(tabla.paper.beta, digits=3, align='c|c|c|c|', caption='Estimates and SD'), type='latex')

###### CON 3 digitos
##################################
\begin{table}[ht]
\centering
\begin{tabular}{c|c|c|c|}
  \hline
 & Classical & Robust & Classical on clean data \\
  \hline
$\widehat{\beta}$   & -0.575 & -0.510 & -0.473 \\
  $\widehat{s}_{\widehat{\beta}} $ & 0.112 & 0.015 & 0.074 \\
   \hline
\end{tabular}
\caption{Estimates and SD}
\end{table}





####### SPLINES CÚBICOS

#--- PRIMERA ---#
# Es decir, usando la B y la D en el paso previo para hallar el Sigma de la DA

# Comenzamos estimando la h* con un estimador robusto:

# Para ello, definimos primer los residuos estandarizados:
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Usaremos splines cúbicos
degree.spline <- 3

# Como queremos ver dónde ocurre el primer mínimo
# local, graficaremos el RBIC
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)

# Si bien la cantidad de nodos que se hubiese
# elegido es:
nk.rob.am$nknots

# El primer mínimo local lo tenemos en
nknots <- 4

# Luego, modelamos la relación entre Z y las covariables
# utilizando un modelo aditivo. Estimamos usando B-splines
# cuadráticos y un estimador robusto:
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# De esta manera, h* es:
hstar <- as.matrix(fit.z.rob$prediction)
# Que sería equivalente a: as.matrix(fit.z.rob$coef.const+rowSums(fit.z.rob$g.matrix))

# Ahora calculamos la matriz B. Observemos que es de
# qxq y, en este caso, q=1
q <- 1
Btheta <- matrix(0,q,q)
for(i in 1:n){
  Btheta <- Btheta + psi.tukey.derivative(residuos[i])*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Btheta <- Btheta/n
Btheta

# Ahora calculamos la matriz D (también de qxq)
Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + (psi.tukey(residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

# Luego, matriz de covarianzas asintótica es
SigmaDA <- sigma^2*solve(Btheta)%*%Dtheta%*%t(solve(Btheta))
SigmaDA

# Y, por lo tanto, el desvío para el paper es:
desvio.rob.BD <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.BD

#############################################################
# Ahora lo que haremos es calcular el desvío obtenido
# para el estimador clásico con la totalidad de los
# datos.
#####################################################

# Comenzamos nuevamente definiendo sus residuos estandarizados.
sigma <- sd(y-fit.full.cl$prediction)

residuos <- (y-fit.full.cl$prediction)/sigma

# Observemos que estimé el sigma con el sd (por si queremos
# considerar otro estimador)

# Ahora analizamos el BIC al utilizar splines cúbicos
degree.spline <- 3
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
plot(nk.cl.am$grid, nk.cl.am$BIC)

# Si bien la cantidad de nodos que minimiza el BIC es:
nk.cl.am$nknots

# Seleccionaremos según el gráfico:
nknots <- 4

# Luego, estimamos el modelo aditivo que modela Z en función de
# las X's utilizando un estimador clásico y splines cuadráticos
fit.z.cl <- am.cl(Z,X,degree.spline=degree.spline, nknots=nknots)

# De esta manera, h* se define como:
hstar <- as.matrix(fit.z.cl$prediction)

# Ahora calculamos la matriz de qxq B
Btheta <- matrix(0,q,q)
for(i in 1:n){
  Btheta <- Btheta + 1*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Btheta <- Btheta/n
Btheta

# Y la matriz D
Dtheta <- matrix(0,q,q)
for(i in 1:n){
  Dtheta <- Dtheta + ((residuos[i]))^2*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
Dtheta <- Dtheta/n
Dtheta

# De esta manera, la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*solve(Btheta)%*%Dtheta%*%t(solve(Btheta))
SigmaDA

# Y, por lo tanto, el desvío para el paper para este caso es:
desvio.cl.full.BD <- sqrt(SigmaDA)/sqrt(n)
desvio.cl.full.BD

##################################################
# Ahora, repetimos este último procedimiento pero para
# el estimador clásico que se obtiene cuando los outliers
# fueron removidos.

# Observemos que ahora tenemos menos datos, por lo que
# tenemos otros n
n.del <- length(y.del)

# Los nuevos residuos son:
sigma.del <- sd(y.del-fit.del.cl$prediction)
residuos.del <- (y.del-fit.del.cl$prediction)/sigma.del

# Analizamos el BIC usando splines cúbicos
degree.spline <- 3
nk.cl.del.am <- select.nknots.cl.am(Z.del,X.del,degree.spline=degree.spline)
plot(nk.cl.del.am$grid, nk.cl.del.am$BIC)

# Si bien la cantidad total de nodos internos hubiese sido
nk.cl.del.am$nknots
# Vamos a considerar:
nknots <- 4

# Ahora estimamos de forma clásica el modelo aditivo que
# modela Z.del en términos de X.del
fit.z.del <- am.cl(Z.del,X.del,degree.spline=degree.spline, nknots=nknots)

# Y obtenemos la siguiente h*
hstar.del <- as.matrix(fit.z.del$prediction)

# Ahora calculamos B
Btheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Btheta.del <- Btheta.del + 1*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Btheta.del <- Btheta.del/n.del
Btheta.del
# Obsérvese que en este caso dividimos por n.del y no por n.

# Y calculamos D
Dtheta.del <- matrix(0,q,q)
for(i in 1:n.del){
  Dtheta.del <- Dtheta.del + ((residuos.del[i]))^2*(Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
Dtheta.del <- Dtheta.del/n.del
Dtheta.del

# De esta manera, la nueva matriz de covarianzas asintótica es:
SigmaDA.del <- sigma.del^2*solve(Btheta.del)%*%Dtheta.del%*%t(solve(Btheta.del))
SigmaDA.del

# Y, por lo tanto, el desvío para el paper es:
desvio.cl.clean.BD <- sqrt(SigmaDA.del)/sqrt(n.del)
desvio.cl.clean.BD

# La primera tabla es entonces
col2 <- c(desvio.cl.full.BD, desvio.rob.BD, desvio.cl.clean.BD)
tabla.desBD <- rbind(col2) #cbind(col2)
row.names(tabla.desBD) <- c('Desvios')
colnames(tabla.desBD) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.desBD,caption = "Table 1")


#--- SEGUNDA ---#
# Es decir, usando directamente la definición de A

# Comenzamos estimando la h* con un estimador robusto:

# Calculamos los residuos
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Analizamos el RBIC como antes (de hecho, da lo mismo)
# para seleccionar la cantidad de nodos
degree.spline <- 3
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)
nknots <- 4

# Nuevamente, estimamos el modelo aditivo que relaciona Z con X
# y realizamos la estimación de forma robusta
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# Luego, el h* es, nuevamente,
hstar <- as.matrix(fit.z.rob$prediction)

# Calculamos ahora el coeficiente "v" según la notación del paper
coef.v <- mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2
coef.v

# Calculamos ahora la matriz A:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n

# Luego, la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*coef.v*solve(AA)
SigmaDA

# Y, por lo tanto, el desvío para el paper es:
desvio.rob.A <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.A

#################################################################
# Ahora repetimos el procedimiento con esl estimador
# clásico que utiliza la totalidad de los datos: con la matriz sigma^2 A
#################################################################

# Calculamos primero los residuos y luego analizamos
# el BIC tal como lo hacíamos antes:
sigma <- sd(y-fit.full.cl$prediction)


degree.spline <- 3
nk.cl.am <- select.nknots.cl.am(Z,X,degree.spline=degree.spline)
plot(nk.cl.am$grid, nk.cl.am$BIC)
nknots <- 4

# Estimamos de forma clásica la relación entre Z y X
# suponiendo un modelo aditivo
fit.z <- am.cl(Z,X,degree.spline=degree.spline, nknots=nknots)

# Luego, nuevamente h* es:
hstar <- as.matrix(fit.z$prediction) #as.matrix(fit.z$coef.const+rowSums(fit.z$g.matrix))

# El coeficiente "v" en este caso se estima como:
residuos <- (y-fit.full.cl$prediction)/sigma

coef.v <- mean((residuos)^2)  #deberia dar cerca de 1
coef.v
sigma ^2

# Ahora calculamos la matriz A:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + (Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/n

# La matriz de covarianzas asintótica estimada es:
SigmaDA <- sigma^2*solve(AA)
SigmaDA

# Finalmente, el desvío para este caso es:
desvio.cl.full.A <- sqrt(SigmaDA)/sqrt(n)
desvio.cl.full.A

# Luego, repetimos el análisis con el estimador
# clásico con los datos "limpios" (sin outliers)

# Calculamos los residuos
n.del <- length(y.del)
sigma.del <- sd(y.del-fit.del.cl$prediction)
residuos.del <- (y.del-fit.del.cl$prediction)/sigma.del

# Nuevamente, miramos el BIC para hallar la cantidad de
# nodos
degree.spline <- 3 #Si se elige 3 da aún más grande
nk.cl.am <- select.nknots.cl.am(Z.del,X.del,degree.spline=degree.spline)
plot(nk.cl.am$grid, nk.cl.am$BIC)
nknots <- 4

# Realizamos el ajuste del modelo aditivo
fit.z.del <- am.cl(Z.del,X.del,degree.spline=degree.spline, nknots=nknots)

# Y entonces obtenemos el siguiente estimador de h*:
hstar.del <- as.matrix(fit.z.del$prediction)

# El coeficiente "v" es, en este caso:
coef.v.del <- mean((residuos.del)^2)
coef.v.del      #DEBERIA SER CERCANO A 1

sigma.del^2

# Calculamos la matriz A estimada:
q <- 1
AA.del <- matrix(0,q,q)
for(i in 1:n.del){
  AA.del <- AA.del + (Z.del[i,]-hstar.del[i,])%*%t(Z.del[i,]-hstar.del[i,])
}
AA.del <- AA.del/n.del

# Luego, el estimador de la matriz de covarianzas asintótica es:
SigmaDA.del <- sigma.del^2*solve(AA.del)
SigmaDA.del

# Y, por lo tanto, el desvío para el paper es:
desvio.cl.clean.A <- sqrt(SigmaDA.del)/sqrt(n.del)
desvio.cl.clean.A


# La segunda tabla es entonces
col2 <- c(desvio.cl.full.A, desvio.rob.A, desvio.cl.clean.A)
tabla.des <- rbind(col2) #cbind(col2)
row.names(tabla.des) <- c('Desvios')
colnames(tabla.des) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.des,caption = "Table 1")


#--- TERCERA ---#
# Esta tabla es similar a la anterior pero la matriz
# A para el estimador robusto se estima con pesos:

# Calculamos nuevamente los residuos:
sigma <- fit.rob$sigma.hat
residuos <- (y-fit.rob$prediction)/sigma

# Calculamos los pesos que utilizan estos residuos
pesos <- psi.w(residuos)

# Repetimos el procedimiento que ya hemos
# realizado en dos oportunidades. Primero seleccionamos
# la cantidad de nodos internos mirando el RBIC y luego
# estimamos:
degree.spline <- 3
set.seed(1111)
nk.rob.am <- select.nknots.rob.am(Z,X,degree.spline=degree.spline)
plot(nk.rob.am$grid, nk.rob.am$RBIC)
nknots <- 4
set.seed(1111)
fit.z.rob <- am.rob(Z,X,degree.spline=degree.spline, nknots=nknots)

# Obtenemos el siguiente estimador de h*:
hstar <- as.matrix(fit.z.rob$prediction)

# Volvemos a calcular el estimador de "v" en este caso:
coef.v <- mean(psi.tukey(residuos)^2)/(mean(psi.tukey.derivative(residuos)))^2
coef.v

# Ahora calculamos el estimador de A usando pesos:
q <- 1
AA <- matrix(0,q,q)
for(i in 1:n){
  AA <- AA + pesos[i]*(Z[i,]-hstar[i,])%*%t(Z[i,]-hstar[i,])
}
AA <- AA/(sum(pesos))

# De esta manera, el estimador de la matriz de covarianzas asintótica es:
SigmaDA <- sigma^2*coef.v*solve(AA)
SigmaDA

# Y, por lo tanto, el desvío para este caso es:
desvio.rob.A.pesos <- sqrt(SigmaDA)/sqrt(n)
desvio.rob.A.pesos

# Luego, la tercera tabla quedaría:
col2 <- c(desvio.cl.full.A, desvio.rob.A.pesos, desvio.cl.clean.A)
tabla.des <- rbind(col2) #cbind(col2)
row.names(tabla.des) <- c('Desvios')
colnames(tabla.des) <- c("Classical","Robust", "Classical on clean data")
library(knitr)
kable(tabla.des,caption = "Table 2,  estimo usando Sigma=sigma^2 v A^{-1}")

kable(tabla.desBD,caption = "Table 1, estimo usando Sigma=B^{-1} D B^{-1}")

# LO QUE SE DEBE REPORTAR EN EL PAPER

pirulo=rbind(tabla.b1.aux[2,], c(tabla.des [1],   tabla.desBD[2], tabla.des[3]))
tabla.paper.beta= round( pirulo, 4)
row.names(tabla.paper.beta)=c(  "$\\widehat{\\beta}$", "\\widehat{s}_{\\widehat{\\beta}}")
colnames(tabla.paper.beta) <- c("Classical","Robust","Classical on clean data")
kable(tabla.paper.beta,caption = "Estimates and SD")


alpha <- 0.99
zalpha <- qnorm(1-(1-alpha)/2,0,1)

round(tabla.paper.beta[1,]-zalpha*tabla.paper.beta[2,],4)
round(tabla.paper.beta[1,]+zalpha*tabla.paper.beta[2,],4)


tabla.paper.beta[1,]-3*tabla.paper.beta[2,]
tabla.paper.beta[1,]+3*tabla.paper.beta[2,]

#             Classical                  Robust Classical on clean data
#                -0.9095                 -0.5543                 -0.7286

# EL EXTREMO del IC PARA EL ROBUSTO NO INCLUYE AL CLASICO!!

tabla.paper.beta[1,]+3*tabla.paper.beta[2,]
#            Classical                  Robust Classical on clean data
#                -0.2405                 -0.4655                 -0.2180

