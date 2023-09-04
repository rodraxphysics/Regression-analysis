install.packages("carData")
install.packages("ggpubr")
library(ggpubr)
library(corrplot)
library(readr)
library(car)
library(ggplot2)
library(ggplot)
library(readr)

#Importacion del dataset
datos <- read_csv("C:/Users/Rodrigo Sandoval/Documents/R/AG1-ESTADISTICA/materials.csv")

#Transformacion variables booleanas a numericas
datos <- transform(datos, is_stable = as.numeric(datos$is_stable))
datos <- transform(datos, is_magnetic = as.numeric(datos$is_magnetic))

#Eliminacion de variables cualitativas
datos <- datos[, -c(1,2)]

#Adicion del band gap como variable booleana
datos$band_gap_bin <- ifelse(datos$band_gap != 0, 1, 0)

############################### REGRESION LINEAL #################################
linear_model = lm(energy_per_atom~uncorrected_energy_per_atom, data=datos)
summary(linear_model) # R-squared:  0.9794   p-value: < 2.2e-16

#Intervalos de confianza
predicciones <- predict(linear_model, newdata = datos, interval = "confidence")
predicciones

# Grafico Datos vs Regresion Lineal
plot(energy_per_atom ~ uncorrected_energy_per_atom, data = datos,
     xlab = "Uncorrected energy per atom",
     ylab = "Energy per atom", main = "Datos vs Regresion Lineal",
     pch = 20, cex = 2,
     col = "grey")
abline(linear_model, col = "red")
legend("topleft", legend = c("Datos", "Regresión lineal"), 
       col = c("grey", "red"), pch = c(20, NA), lty = c(NA, 1), cex =1.1)

# Graficos de Diagnostico
par(mfrow=c(2,2))
plot(linear_model, which=1)
plot(linear_model, which=2)
plot(linear_model, which=3)
plot(linear_model, which=4)
par(mfrow=c(1,1))

########################### REGRESION MULTILINEAL #################################

#Analisis de correlaciones 
cor_matrix <- cor(datos, method = "pearson")
cor_df <- as.data.frame(as.table(cor_matrix))
cor_df <- subset(cor_df, Var1 != Var2)
cor_df <- cor_df[order(-abs(cor_df$Freq)),]
cor_df

# Grafico visual de correlaciones
datos_corr=datos
colnames(datos_corr) <- c("y_1", "x_a", "x_b", "x_c", "x_d", "x_e", "x_f", "x_g", "x_h", "x_i", "x_j", "x_k", "x_l", "x_m", "x_n", "x_o", "x_p", "x_q", "x_r", "x_s", "y_2")
M = cor(datos_corr)
corrplot.mixed(M, order = "AOE")

#Regresion Multilineal

# Modelo de 10 variables independientes (con el mayor coeficiente de determinacion)
modelo <- lm(band_gap ~ nelements + density +density_atomic + crystal_symmetry + +sides_abc +formation_energy_per_atom +efermi +total_magnetization +num_magnetic_sites +atomic_number, data = datos)
summary(modelo) #R-squared:  0.56  p-value: < 2.2e-16

# Modelo de 10 variables independientes con abreviacion
modelo <- lm(y_1 ~ x_b + x_d + x_e + x_f  + x_h + x_l + x_o + x_q + x_r + x_s, data = datos_corr)
summary(modelo)

# Calculo de coeficientes VIF de las 13 variables
vif_values <- vif(modelo)
barplot(vif_values , main = "VIF Values", horiz = TRUE, col = "steelblue")
abline(v = 5, lwd = 3, lty = 2)

# Modelo de 2 variables independientes (menos propenso a overfitting)
modelo <- lm(band_gap ~ efermi + formation_energy_per_atom , data = datos)
summary(modelo) #R-squared:  0.5159  p-value: < 2.2e-16

#Plots de Diagnostico
par(mfrow=c(2,2))
plot(modelo, which=1)
plot(modelo, which=2)
plot(modelo, which=3)
plot(modelo, which=4)
par(mfrow=c(1,1))

#Grafica Datos vs Regresion Multilineal
rango_efermi <- range(datos$efermi)
nuevos_valores_efermi <- seq(from = rango_efermi[1], to = rango_efermi[2], length.out = 20)

rango_formation_energy_per_atom <- range(datos$formation_energy_per_atom)
nuevos_valores_formation_energy_per_atom <- seq(from = rango_formation_energy_per_atom[1], to = rango_formation_energy_per_atom[2], length.out = 20)

predicciones <- outer(X = nuevos_valores_efermi, Y = nuevos_valores_formation_energy_per_atom,
                      FUN = function(efermi, formation_energy_per_atom) {
                        predict(object = modelo,
                                newdata = data.frame(efermi, formation_energy_per_atom))
                      })


superficie <- persp(x = nuevos_valores_efermi, y = nuevos_valores_formation_energy_per_atom,
                    z = predicciones,
                    theta = 20, phi = 15,
                    col="#4682B4", shade = 0.1,
                    zlim = range(-15,25),
                    xlab = "efermi", ylab = "formation_energy_per_atom", zlab = "band_gap",
                    ticktype = "detailed",
                    main = "Datos vs Regresion Multilineal"
)

observaciones <- trans3d(datos$efermi, datos$formation_energy_per_atom, datos$band_gap, superficie)
error <- trans3d(datos$efermi, datos$formation_energy_per_atom, fitted(modelo), superficie)
points(observaciones, col = "red", pch = 20)
segments(observaciones$x, observaciones$y, error$x, error$y)

########################### REGRESION LOGISTICA #################################
log_model <- glm(band_gap_bin ~ efermi+formation_energy_per_atom +density + +total_magnetization +num_magnetic_sites, data=datos, family = binomial)
summary(log_model) 

# Pseudo coeficiente de determinacion
model_deviance <- log_model$deviance
null_deviance <- log_model$null.deviance
pseudo_R2_mcfadden <- 1 - (model_deviance / null_deviance)
print(paste("Pseudo R^2 de McFadden: ", pseudo_R2_mcfadden)) #0.47

# Regresino logistica con una variable independiente
log_model <- glm(band_gap_bin ~ efermi, data=datos, family = binomial)
summary(log_model) 

#Plot Datos vs Regresion Logistica
ggplot(datos, aes(x=efermi, y=band_gap_bin)) +
  geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE) +
  ggtitle("Datos vs Regresion Logistica") +
  xlab("efermi") +
  ylab("band_gap_bin")

# Pseudo coeficiente de determinacion
model_deviance <- log_model$deviance
null_deviance <- log_model$null.deviance
pseudo_R2_mcfadden <- 1 - (model_deviance / null_deviance)
print(paste("Pseudo R^2 de McFadden: ", pseudo_R2_mcfadden)) #0.357


