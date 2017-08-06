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

