require(wesanderson)
colores <- wes_palette("Cavalcanti", 3)

plot(c(0,5), c(0,2), type="n", axes = FALSE, xlab = "Intensidad", ylab = "Densidad Exponencial")
axis(1)
axis(2)
curve(dexp(x, rate = 2, ), col=colores[3], lwd=3, from=0, to=5, add=TRUE)
curve(dexp(x, rate = 1), col=colores[1], lwd=3, add=TRUE)
curve(dexp(x, rate = 1/3, ), col=colores[2], lwd=3, add=TRUE)
legend("topright", lwd = c(3,3,3), col = colores, legend = c("Media 1", "Media 3", "Media 0.1"))


plot(c(0,5), c(10^(-6),10), type="n", axes = FALSE, log="y", xlab = "Intensidad", ylab = "Densidad Exponencial",
     main = "Densidades Exponenciales - Escala Semilogarítmica")
axis(1)
axis(2)
curve(dexp(x, rate = 2), col=colores[3], lwd=3, from=0.001, to=5, add=TRUE)
curve(dexp(x, rate = 1), col=colores[1], lwd=3, add=TRUE)
curve(dexp(x, rate = 1/3), col=colores[2], lwd=3, add=TRUE)
legend("bottomleft", lwd = c(3,3,3), col = colores, legend = c("Media 1", "Media 3", "Media 0.1"))

dgammaSAR <- function(x, Looks, mean, log=FALSE) {
  dgamma(x, shape=Looks, rate=Looks/mean, log=log)
}

plot(c(0,5), c(0,1), type="n", axes = FALSE, xlab = "Intensidad", ylab = "Densidades Gama de Media Unitaria",
     main = "Densidades del Speckle - Escala Lineal")
axis(1)
axis(2)
curve(dgammaSAR(x, 1, 1), col=colores[3], lwd=3, from=0, to=5, add=TRUE)
curve(dgammaSAR(x, 3, 1), col=colores[2], lwd=3, add=TRUE)
curve(dgammaSAR(x, 5, 1), col=colores[1], lwd=3, add=TRUE)
legend("topright", lwd = c(3,3,3), col = colores[c(3,2,1)], legend = c("L = 1", "L = 3", "L = 5"))

plot(c(0,5), c(10^(-6),1), type="n", axes = FALSE, log = "y", xlab = "Intensidad", ylab = "Densidades Gama de Media Unitaria",
       main = "Densidades del Speckle - Escala Logarítmica")
axis(1)
axis(2)
curve(dgammaSAR(x, 1, 1), col=colores[3], lwd=3, from=0, to=5, add=TRUE)
curve(dgammaSAR(x, 3, 1), col=colores[2], lwd=3, add=TRUE)
curve(dgammaSAR(x, 5, 1), col=colores[1], lwd=3, add=TRUE)
legend("topright", lwd = c(3,3,3), col = colores[c(3,2,1)], legend = c("L = 1", "L = 3", "L = 5"))

dGI0 <- function(z, p_alpha, p_gamma, p_Looks, log=FALSE) {
  
  if(log==TRUE) {
    return(
      (p_Looks*log(p_Looks) + lgamma(p_Looks-p_alpha) + (p_Looks-1)*log(z) ) - 
        (p_alpha*log(p_gamma) + lgamma(-p_alpha) + lgamma(p_Looks) + 
           (p_Looks-p_alpha)*log(p_gamma + z*p_Looks) ) 
    )   
  }
  else { return( 
    ( p_Looks^p_Looks * gamma(p_Looks-p_alpha) * z^(p_Looks-1) ) / 
      (p_gamma^p_alpha * gamma(-p_alpha) * gamma(p_Looks) * (p_gamma + z*p_Looks)^(p_Looks-p_alpha)) 
  )
  }
}

plot(c(0,5), c(0,3), type="n", axes = FALSE, xlab = "Intensidad", ylab = "Densidades GI0 de Media Unitaria",
     main = "Densidades del Retorno - Escala Lineal")
axis(1)
axis(2)
curve(dGI0(x, -1.5, 0.5, 1), col=colores[3], lwd=3, n = 500, from=0, to=5, add=TRUE)
curve(dGI0(x, -1.5, 0.5, 3), col=colores[2], lwd=3, n = 500, add=TRUE)
curve(dGI0(x, -1.5, 0.5, 8), col=colores[1], lwd=3, n = 500, add=TRUE)
legend("topright", lwd = c(3,3,3), col = colores[c(3,2,1)], legend = c("L = 1", "L = 3", "L = 5"))

plot(c(0,5), c(10^(-2),10), type="n", axes = FALSE, log = "y", xlab = "Intensidad", ylab = "Densidades GI0 de Media Unitaria",
     main = "Densidades del Retorno - Escala Semilogarítmica")
axis(1)
axis(2)
curve(dGI0(x, -1.5, 0.5, 1), col=colores[3], lwd=3, n = 500, from=0, to=5, add=TRUE)
curve(dGI0(x, -1.5, 0.5, 3), col=colores[2], lwd=3, n = 500, add=TRUE)
curve(dGI0(x, -1.5, 0.5, 8), col=colores[1], lwd=3, n = 500, add=TRUE)
legend("topright", lwd = c(3,3,3), col = colores[c(3,2,1)], legend = c("L = 1", "L = 3", "L = 5"))

source("./imagematrix.R")
if(!require(tiff)){install.packages("tiff"); require(tiff)}

barras <- readTIFF("strips.tif") # Leemos la imagen de barras
dim(barras)
range(barras) # Verificamos el rango de sus valores
plot(imagematrix(barras)) # La visualizamos

barras <- barras / max(barras) * (.9 - .1) + .1
barras_1Look <- barras * matrix(data=rexp(n = 256*256, rate = 1), ncol=256, nrow = 256)
plot(imagematrix(normalize(barras_1Look)))
plot(imagematrix(equalize(barras_1Look)))

plot(barras_1Look[, 214], type="n", xlab = "Línea", ylab = "Valores observados", axes = FALSE,
     main = "Valores en la Columna 214: verdad y speckle 1-Look")
axis(1)
axis(2)
lines(barras_1Look[, 214], col=colores[1])
lines(barras[, 214],col=colores[2])
