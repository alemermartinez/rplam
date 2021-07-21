#Forestfires

datos <- read.csv("R/forestfires.csv", header = TRUE, sep=",")
str(datos)

y <- datos$area
Z <- cbind(datos$temp, datos$RH, datos$wind, datos$rain)

pairs(cbind(y,Z))

#Tiene que ser discreta
plot(datos$rain,y)
hist(datos$rain)

plot(datos$temp, y)
plot(datos$RH, y)
plot(datos$wind,y)

str(datos)

plot(datos$ISI,y, xlim=c(0,25))
plot(datos$FFMC,y)
plot(datos$DMC,y)
plot(datos$DC,y)


#No llegué a ninguna conclusión
