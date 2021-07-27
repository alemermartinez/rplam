

#############################
# Borro todo por seguridad
#############################

#rm(list=ls())

#library(devtools)
#install_github("alemermartinez/rplam")
#library(rplam)

##########################################################
#Llamamos al archivo que tiene el programa de la simulacion para grado 1
##########################################################

#source("simu_splines_clasico_new.R")
source("simu_splines_clasico_new-newm3.R")

#############################
#Tama?o de la muestra
#############################

n <- 100

############################
# EJEMPLO 1, 2, 3, 4 o 5
############################

ejemplo <- 5

#############################
#Desvio de los errores
#############################

desvio.epsilon <- 0.2

#Ejemplos 1, 4, 5, 6, 7, 8, 9, 10: 0.2
#Ejemplos 2 y 3 no los uso: 0.5

#############################
#Epsilon de corte para el calculo de las estimaciones
#############################

epsilon <- 1e-10

#############################
# Grado del polinomio
#############################

degree.spline <- 3

#############################
# FIJO EL TIPO DE CONTAMINACION: 0, 1, 2 o 3
#############################

for (tipo.cont in 0:3){
  
#############################
# FIJO CUANTAS ITERACIONES HAGO
#############################

NITER <- 500 

#############################
# FIJO DONDE TERMINO EL ANTERIOR
#############################

anterior<- 0

tiempo1<-proc.time()

try( simulacion(ejemplo, NITER,anterior,tipo.cont,n,desvio.epsilon,epsilon,degree.spline) )

tiempo2<-proc.time()

print(tiempo2-tiempo1)

}

