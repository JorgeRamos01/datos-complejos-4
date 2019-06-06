library(geoR)

#Explorando los datos
data(parana)
help(parana)
summary(parana)

plot(parana, lowess = T)

#Generamos el variograma para poder observar la senal subyacente
parana.vario <- variog(parana, max.dist=400)
plot(parana.vario)

#Ajustamos un modelo exponencial y un modelo Matern con kappa=1.5
parana.exp <- variofit(parana.vario)
parana.mat<- variofit(parana.vario, kappa=1.5)

#Graficamos los modelos contra el variograma
plot(parana.vario, main="Ajuste de modelos", xlab="Distancia", ylab="Semivarianza")
lines(parana.exp, col="blue")
lines(parana.mat, col="red")
legend(0,5000,legend=c("Exponencial", "Matern"), col=c("blue", "red"), lty=1)
