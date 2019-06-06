library(geoR)
library(RandomFields)
library(maptools)
library(rgdal)
library(proj4)

####Cargamos los datos
setwd("~/Documentos/Datos complejos/Tarea 4")
Obe<-readOGR(dsn="obesidad/casos.dbf",layer = "casos")

Poli<-readOGR(dsn="obesidad/poligono_region.dbf")

Poli2<-readOGR(dsn="obesidad/poligono_region_mun.dbf")

NL<-readOGR(dsn="obesidad/nl_ageb_urb.dbf")

polys <- Poli@polygons[[1]]@Polygons
coords.original <- coordinates(polys[[1]])
p.map <- list(Polygons(list(Polygon(coords.original)), "coords.polys"))
poly.map <- SpatialPolygons(p.map)
plot(poly.map,axes=T)

coords.borders <- poly.map@polygons[[1]]@Polygons[[1]]@coords

#############################
######Consideramos los casos primero por peso
############################

mod <- "gaussian"
c.par <- c(1,.05)
kap.par <- .5
mean.par <- 35
nug.par <- .1
set.seed(7152)
data<-as.data.frame(cbind(Obe@coords,Obe@data$Peso))

data.geo <- as.geodata(data)
data.geo<- jitterDupCoords(data.geo,1)

#Exploramos los datos
summary(data.geo)
#Los graficamos
plot(poly.map,axes=T, main="Casos: Peso")
points(data.geo$coords, col="red", lty=2)

#Obtenemos su variograma
peso.vario <- variog(data.geo, max.dist=0.1,lambda = 1)
plot(peso.vario, main="Semivariograma: Peso")
#Intentamos ajustar un modelo
peso.vario.fit <- variofit(peso.vario, cov.model = "matern", kappa=1.5)
lines.variomodel(peso.vario.fit)

#########################################
######Consideramos los casos por tallacms
#########################################
data1<-as.data.frame(cbind(Obe@coords,Obe@data$Tallacms))

data.geo1 <- as.geodata(data1)
data.geo1<- jitterDupCoords(data.geo1,1)

#Exploramos los datos
summary(data.geo1)

#Los graficamos
plot(poly.map,axes=T, main="Casos: Tallacms")
points(data.geo1$coords, col="blue", lty=2)

#Obtenemos su variograma
peso.vario1 <- variog(data.geo1, max.dist=0.1,lambda = 1)
plot(peso.vario1, main="Semivariograma: Tallacms")
#Ajustamos un modelo
peso.vario1.fit <- variofit(peso.vario1, cov.model = "matern", kappa=0.5)
lines.variomodel(peso.vario1.fit)


#########################################
######Consideramos los casos por Circintu
#########################################
data2<-as.data.frame(cbind(Obe@coords,Obe@data$Circintu))

data.geo2 <- as.geodata(data2)
data.geo2<- jitterDupCoords(data.geo2,1)

#Exploramos los datos
summary(data.geo2)

#Los graficamos
plot(poly.map,axes=T, main="Casos: Circintu")
points(data.geo2$coords, col="orange", lty=2)

#Obtenemos su variograma
peso.vario2 <- variog(data.geo2, max.dist=0.1,lambda = 1)
plot(peso.vario2, main="Semivariograma: Circintu")
#Ajustamos un modelo
peso.vario2.fit <- variofit(peso.vario2, cov.model = "matern", kappa=2)
lines.variomodel(peso.vario2.fit)

############################################# 
#realizamos Kriging para los datos de peso
##############################################

## definimos las ubicaciones para prediccion en un grid
grid <- pred_grid(coords.borders,by=.002)
par(mfrow=c(2,2))
for(i in seq(0.5,2,0.5)){
  ## parametros del modelo de prediccion (kriging ordinario)
  krige <- krige.control(type.krige="OK",trend.d="1st",trend.l="1st",
                       obj.model = NULL,cov.model="matern", kappa=i,
                       cov.pars=c.par,nugget=5,micro.scale=0,
                       dist.epsilon=1e-10)
  pred <- krige.conv(geodata=data.geo,loc=grid,krige=krige,borders=coords.borders)
  texto<-paste("Peso, Matern kappa=", i)
  image(pred, main=texto)
  contour(pred,add=TRUE)
  legend.krige(x.leg=c(-100.5,-100.45),y.leg=c(25.6,25.63),
             pred$pred,vert=TRUE,offset.leg=.3)
  points(data,col="blue")
}

par(mfrow=c(1,1))
## parametros del modelo de prediccion (kriging ordinario) con cov.model= gaussian
krige <- krige.control(type.krige="OK",trend.d="1st",trend.l="1st",
                       obj.model = NULL,cov.model="gaussian",
                       cov.pars=c.par,nugget=5,micro.scale=0,
                       dist.epsilon=1e-10)
pred <- krige.conv(geodata=data.geo,loc=grid,krige=krige,borders=coords.borders)
image(pred, main="Peso, Gaussian")
contour(pred,add=TRUE)
legend.krige(x.leg=c(-100.5,-100.45),y.leg=c(25.6,25.63),
             pred$pred,vert=TRUE,offset.leg=.3)
points(data,col="blue")

############################################# 
#realizamos Kriging para los datos de Tallacms
##############################################

## definimos las ubicaciones para prediccion en un grid
grid <- pred_grid(coords.borders,by=.002)

  ## parametros del modelo de prediccion (kriging ordinario) con cov.model= matern
krige <- krige.control(type.krige="OK",trend.d="1st",trend.l="1st",
                         obj.model = NULL,cov.model="matern", kappa=0.5,
                         cov.pars=c.par,nugget=5,micro.scale=0,
                         dist.epsilon=1e-10)
pred <- krige.conv(geodata=data.geo1,loc=grid,krige=krige,borders=coords.borders)
image(pred, main="Tallacms, Matern kappa=0.5")
contour(pred,add=TRUE)
legend.krige(x.leg=c(-100.5,-100.45),y.leg=c(25.6,25.63),
               pred$pred,vert=TRUE,offset.leg=.3)
points(data1,col="blue")


############################################# 
#realizamos Kriging para los datos de Circuncin
##############################################

## definimos las ubicaciones para prediccion en un grid
grid <- pred_grid(coords.borders,by=.002)

## parametros del modelo de prediccion (kriging ordinario) con cov.model= matern
krige <- krige.control(type.krige="OK",trend.d="1st",trend.l="1st",
                         obj.model = NULL,cov.model="matern", kappa=1.5,
                         cov.pars=c.par,nugget=5,micro.scale=0,
                         dist.epsilon=1e-10)
pred <- krige.conv(geodata=data.geo2,loc=grid,krige=krige,borders=coords.borders)
image(pred, main="Circuncin, Matern kappa=1.5")
contour(pred,add=TRUE)
legend.krige(x.leg=c(-100.5,-100.45),y.leg=c(25.6,25.63),
               pred$pred,vert=TRUE,offset.leg=.3)
points(data2,col="blue")



