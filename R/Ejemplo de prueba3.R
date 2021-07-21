library(robustbase)

#install.packages("RSADBE")
library(RSADBE)
data(SPD)
str(SPD)

source("R/rplam-fn.R")

##########
#Prueba 1#
##########

y <- SPD$Y
Z <- cbind(SPD$X1,SPD$X2,SPD$X3,SPD$X4)
X <- cbind(SPD$X5,SPD$X6)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$","$\\beta_3$","$\\beta_4$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",30),rep("Robust",30))
))
names(DF2) <- c("x5","x6","re.x5", "re.x6", "gx5","gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:30,
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

#
(DF3$Residuals==out[1])[7]
(DF3$Residuals==out[2])[23]
#Posiciones 7 y 23

DF3 <- as.data.frame(list(
  X,
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1],
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2],
  rbind(fit.rob$g.matrix),
  rep("Robust",30)
))
names(DF3) <- c("x5","x6","re.x5", "re.x6", "gx5","gx6","Fit")

DF3point <- as.data.frame(list(
  X[c(7,23),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1])[c(7,23),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2])[c(7,23),]
)
)
names(DF3point) <- c("x5.x","x6.x","re.x5","re.x6")

ggplot(DF3,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


out.pos <- c(7, 23)
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- Z[-out.pos,]

fit.del.cl <- plam.cl(y.del, Z.del, X.del)
fit.del.cl$nknots

col4 <- c(fit.del.cl$alpha,fit.del.cl$coef.lin)
tabla.b1.aux <- cbind(tabla.b1,col4)
row.names(tabla.b1.aux) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$","$\\beta_3$","$\\beta_4$")
colnames(tabla.b1.aux) <- c("Classical","Robust","Classical on clean data")
kable(tabla.b1.aux,caption = "Table 2")

DF4 <- as.data.frame(list(
  rbind(X.del,X),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - fit.del.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - fit.del.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",28),rep("Robust",30))
))
names(DF4) <- c("x5","x6","re.x5", "re.x6", "gx5","gx6","Fit")

ggplot(DF4,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,28),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,28),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )



#########
#Prueba2#
#########

y <- SPD$Y
Z <- cbind(SPD$X1,SPD$X3,SPD$X4)
X <- cbind(SPD$X2,SPD$X5,SPD$X6)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$","$\\beta_3$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-1]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-2]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-3]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3])),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",30),rep("Robust",30))
))
names(DF2) <- c("x4", "x5","x6","re.x4", "re.x5", "re.x6", "gx4", "gx5","gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x4,y=re.x4,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x4,y=gx4),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:30,
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

#No tiene outliers



#########
#Prueba3#
#########

y <- SPD$Y
Z <- cbind(SPD$X1)
X <- cbind(SPD$X2,SPD$X3,SPD$X4,SPD$X5,SPD$X6)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-1]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-2]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-3]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-4]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-4])),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - rowSums(fit.full.cl$g.matrix[,-5]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-5])),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",30),rep("Robust",30))
))
names(DF2) <- c("x2","x3", "x4", "x5","x6","re.x2", "re.x3","re.x4", "re.x5", "re.x6", "gx2", "gx3", "gx4", "gx5","gx6","Fit")

library(ggplot2)

ggplot(DF2,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x4,y=re.x4,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x4,y=gx4),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )



res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:30,
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

#
(DF3$Residuals==out[1])[1]
(DF3$Residuals==out[2])[15]
(DF3$Residuals==out[3])[20]
(DF3$Residuals==out[4])[23]
(DF3$Residuals==out[5])[29]

DF3 <- as.data.frame(list(
  X,
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-4]),
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-5]),
  rbind(fit.rob$g.matrix),
  rep("Robust",30)
))
names(DF3) <- c("x2", "x3", "x4", "x5","x6","re.x2", "re.x3", "re.x4", "re.x5", "re.x6", "gx2", "gx3", "gx4", "gx5","gx6","Fit")

DF3point <- as.data.frame(list(
  X[c(1,15,20,23,29),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-1]))[c(1,15,20,23,29),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-2]))[c(1,15,20,23,29),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-3]))[c(1,15,20,23,29),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-4]))[c(1,15,20,23,29),],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-5]))[c(1,15,20,23,29),]
)
)
names(DF3point) <- c("x2.x", "x3.x", "x4.x", "x5.x","x6.x","re.x2", "re.x3", "re.x4", "re.x5","re.x6")

ggplot(DF3,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x4,y=re.x4,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x4,y=gx4),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x4.x,y=re.x4),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(1,30)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


out.pos <- c(1,15,20,23,29)
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- as.matrix(Z[-out.pos,])

fit.del.cl <- plam.cl(y.del, Z.del, X.del)
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
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - rowSums(fit.del.cl$g.matrix[,-4]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-4])),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - rowSums(fit.del.cl$g.matrix[,-5]),fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - rowSums(fit.rob$g.matrix[,-5])),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",25),rep("Robust",30))
))
names(DF4) <- c("x2", "x3", "x4", "x5","x6","re.x2", "re.x3", "re.x4", "re.x5", "re.x6", "gx2", "gx3", "gx4", "gx5","gx6","Fit")

ggplot(DF4,aes(x=x2,y=re.x2,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x2,y=gx2),size=1.3,alpha=1,linetype=c(rep(2,25),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x2.x,y=re.x2),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x3,y=re.x3,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x3,y=gx3),size=1.3,alpha=1,linetype=c(rep(2,25),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x3.x,y=re.x3),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x4,y=re.x4,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x4,y=gx4),size=1.3,alpha=1,linetype=c(rep(2,25),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x4.x,y=re.x4),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


ggplot(DF4,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,25),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,25),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

#Me quedo sólo con las obvias

y <- SPD$Y
Z <- cbind(SPD$X1)
X <- cbind(SPD$X5, SPD$X6)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",30),rep("Robust",30))
))
names(DF2) <- c("x5","x6","re.x5", "re.x6", "gx5", "gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:30,
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

#Tampoco tiene outliers

pairs(SPD)

ss <- lm(SPD$X1~SPD$X2)
summary(ss)

ss <- lm(SPD$X2~SPD$X3)
summary(ss)

ss <- lm(SPD$X4~SPD$X3) #Están bastante correlacionadas linealmente
summary(ss)

ss <- lm(SPD$X4~SPD$X6)
summary(ss)

#########
#Prueba4#
#########
y <- SPD$Y
Z <- cbind(SPD$X1,SPD$X3) #Descarto x2 por poca correlación y x4 por repetida
X <- cbind(SPD$X5,SPD$X6)

degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",30),rep("Robust",30))
))
names(DF2) <- c("x5","x6","re.x5", "re.x6", "gx5", "gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,30),rep(1,30)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

res.rob <- y-fit.rob$prediction
summary(res.rob)

library(cowplot)

DF3 <- as.data.frame(list(
  1:30,
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

