#install.packages("SemiPar")
library(SemiPar)
data(ragweed)
pairs(ragweed)

#Saco a la lluvia porque es discreta

y <- ragweed$ragweed
Z <- cbind(ragweed$wind.speed)
X <- cbind(ragweed$day.in.seas, ragweed$temperature)

plot(ragweed$day.in.seas,ragweed$ragweed)
plot(ragweed$temperature, ragweed$ragweed)
plot(ragweed$wind.speed, ragweed$ragweed)


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

n <- length(y)
DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n),rep("Robust",n))
))
names(DF2) <- c("x5","x6","re.x5", "re.x6", "gx5", "gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
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


################################
##Agrego discreta que es lluvia#
################################

y <- ragweed$ragweed
Z <- cbind(ragweed$wind.speed, ragweed$rain)
X <- cbind(ragweed$day.in.seas, ragweed$temperature)


plot(ragweed$day.in.seas,ragweed$ragweed)
plot(ragweed$temperature, ragweed$ragweed)
plot(ragweed$wind.speed, ragweed$ragweed)
plot(ragweed$rain, ragweed$ragweed)


degree.spline <- 3

nk.cl <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
nk.cl$nknots

nk.rob <- select.nknots.rob(y,Z,X,degree.spline=degree.spline)
nk.rob$nknots

fit.full.cl <- plam.cl(y,Z,X,degree.spline=degree.spline)

set.seed(123)
fit.rob <- plam.rob(y,Z,X,degree.spline=degree.spline)


col2 <- c(fit.full.cl$alpha,fit.full.cl$coef.lin)
col3 <- c(fit.rob$alpha,fit.rob$coef.lin)
tabla.b1 <- cbind(col2,col3)
row.names(tabla.b1) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$")
colnames(tabla.b1) <- c("Classical","Robust")
library(knitr)
kable(tabla.b1,caption = "Table 1")

n <- length(y)
DF2 <- as.data.frame(list(
  rbind(X,X),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.full.cl$y  - fit.full.cl$Z%*%fit.full.cl$coef.lin - fit.full.cl$alpha - fit.full.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.full.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n),rep("Robust",n))
))
names(DF2) <- c("x5","x6","re.x5", "re.x6", "gx5", "gx6","Fit")

library(ggplot2)
ggplot(DF2,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF2,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,n),rep(1,n)))+
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
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1],
  fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2],
  rbind(fit.rob$g.matrix),
  rep("Robust",n)
))
names(DF3) <- c("x5","x6","re.x5", "re.x6", "gx5","gx6","Fit")

DF3point <- as.data.frame(list(
  X[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1])[outliers,],
  (fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2])[outliers,]
)
)
names(DF3point) <- c("x5.x","x6.x","re.x5","re.x6")

ggplot(DF3,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF3,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(1,n)))+
  scale_color_manual(values = c("#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )


out.pos <- outliers
y.del <- y[-out.pos]
X.del <- X[-out.pos,]
Z.del <- Z[-out.pos,]

fit.del.cl <- plam.cl(y.del, Z.del, X.del)
fit.del.cl$nknots

col4 <- c(fit.del.cl$alpha,fit.del.cl$coef.lin)
tabla.b1.aux <- cbind(tabla.b1,col4)
row.names(tabla.b1.aux) <- c('$\\beta_0$', "$\\beta_1$", "$\\beta_2$")
colnames(tabla.b1.aux) <- c("Classical","Robust","Classical on clean data")
kable(tabla.b1.aux,caption = "Table 2")

DF4 <- as.data.frame(list(
  rbind(X.del,X),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - fit.del.cl$g.matrix[,-1],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-1]),
  c(fit.del.cl$y  - fit.del.cl$Z%*%fit.del.cl$coef.lin - fit.del.cl$alpha - fit.del.cl$g.matrix[,-2],fit.rob$y  - fit.rob$Z%*%fit.rob$coef.lin - fit.rob$alpha - fit.rob$g.matrix[,-2]),
  rbind(fit.del.cl$g.matrix,fit.rob$g.matrix),
  c(rep("Classical",n-length(out.pos)),rep("Robust",n))
))
names(DF4) <- c("x5","x6","re.x5", "re.x6", "gx5","gx6","Fit")

ggplot(DF4,aes(x=x5,y=re.x5,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x5,y=gx5),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x5.x,y=re.x5),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )

ggplot(DF4,aes(x=x6,y=re.x6,color=Fit)) +
  geom_point(alpha=0.5)+
  geom_line(aes(x=x6,y=gx6),size=1.3,alpha=1,linetype=c(rep(2,n-length(out.pos)),rep(1,n)))+
  scale_color_manual(values = c("#bb0c00","#0052bb"))+
  geom_point(data=DF3point,aes(x=x6.x,y=re.x6),col='#FF33E3',size=2)+
  theme(
    axis.title.y=element_blank()
  )



