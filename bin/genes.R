#Definiendo los paquetes a necesitar
instala_paquetes<-function(){
  pkgs_CRAN<-c("R.utils", "readr")
  pkgs_Bioc<-c("limma", "edgeR")
  chooseCRANmirror(ind = 55)
  #Ahora iteramos por cada paquete
  
  if(!("BiocManager" %in% installed.packages()[,"Package"])){
    writeLines("El paquete BiocManager se va a instalar")
    install.packages("BiocManager",dependencies = T)
  }

  BiocManager::install(version = "3.14")
  pkgs_faltantes<-pkgs_CRAN[!(pkgs_CRAN %in% installed.packages()[,"Package"])] 
  if(length(pkgs_faltantes)>0){
    writeLines("Instalando paquetes necesarios")
    install.packages(pkgs_faltantes, dependencies = T)
  }
  
  pkgs_faltantes<-pkgs_Bioc[!(pkgs_Bioc %in% installed.packages()[,"Package"])] 
  if(length(pkgs_faltantes)>0){
    writeLines("Instalando paquetes necesarios")
    install.packages(pkgs_faltantes, dependencies = T)
  }
  #cargando todos los paqutes
  invisible(lapply(c(pkgs_CRAN,pkgs_Bioc), library, character.only = TRUE))
}
#Sys.setlocale(category = "LC_COLLATE", locale = "C")
Sys.setlocale(category = "LC_CTYPE", locale = "C")
#Sys.setlocale(category = "LC_MONETARY", locale = "C")
#Sys.setlocale(category = "LC_NUMERIC", locale = "C")
#Sys.setlocale(category = "LC_TIME", locale = "C")
instala_paquetes()
writeLines("
Hola, bienvenido al script para el analisis de RNA-Seq desarrollado por:
                          Alina y Uriel.
")

defW<-getOption("warn")
options(warn = -1)

cat("El directorio de trabajo actual es: ", getwd()," Es el correcto[s/n]: ")
#readLines(file("stdin"),1) nos permite obtener inputs del usuario en forma
#interactiva, algo similar a python.
raw <- readLines(file("stdin"),1)
if(raw=="n"){
  cat("Ingresa el directorio correcto: ")
  dr<-raw <- readLines(file("stdin"),1)
  setwd(dr)
}
writeLines("\n")
cat("Los archivos en formato v?lido que estan disponibles son: \n")
rfiles<-dir()[grep("\\.csv$|\\.txt$",dir())]
rfilesAv<-paste(paste(1:length(rfiles),".",sep = ""), rfiles, sep = " ")
cat(rfilesAv, sep = "\n")
writeLines("Selecciona un archivo escribiendo el n?mero que le antecede: ")
selFile <- readLines(file("stdin"),1)
selFile<-as.numeric(selFile)
cat(c("Archivo",rfiles[selFile]),sep=": ")
writeLines("\n")

#data<-read.table(rfiles[selFile], header = T)
data<- readr::read_delim(file = rfiles[selFile], delim = ",", col_names = T,
                         show_col_types = F)
str(data)
writeLines("\n")

rm_na<-function(x){
  #Looking for rows (genes) with NAs
  genes_NA<-unique(which(is.na(x))%%dim(x)[1])
  
  #Removing genes with NAs (if there is at least one): 
  if(length(genes_NA)>1){
    #Retrieving and showing names of genes with NA,
    #it is assumed that column 1 stores genes' names:
    writeLines("Los siguientes genes contienen NAs y han sido removidos: ")
    cat(dplyr::pull(x[genes_NA,],1), sep = ", ")
    writeLines("\n")
    x<-na.omit(x)
  }
}


data<-rm_na(data)
#Una vez limpios de NAs, guardamos la informaci�n de los ID de genes en un
#vector y dejamos solamente las columna que corresponden a muestras.
genesID<-data[,1]
data<-data[,-1]



#head(data)
#str(data)
#summary(data)
#print(getwd())



#cat("La matriz de conteos encontrada en el directorio es: ", dir()[1])


#print(wel)
#d<-dir()[grep("\\.R$",dir())]
#print(d)
#getwd()  
#Con esto obtenemos inputs del usuario
#argus <- commandArgs(trailingOnly = T)
#print(argus)

# RNAseq


### Paquetes necesarios
#--------------------
library(limma)      #
#library(Glimma)     #
library(edgeR)      #
library(R.utils)    #
#--------------------

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
## Añadiendo información de grupos de interés 
group <- as.factor(c("WT", "WT", "ult1", "ult1", "ult1ult2", "ult1ult2"))
samplenames<-c("WT", "WT", "ult1", "ult1", "ult1ult2", "ult1ult2")
## Pasando de números a id de genes como nombre de filas
rownames(x) <- x[,1]
## Eliminando columna redundante
x <- x[,-1]
## Extrayendo de nuevo los numbres de los genes
geneid <- rownames(x)
## Creando la estructura de datos que nos permitirá interactuar con los 
## paquetes cargados previamente.
x <- DGEList(counts = x, group = group, genes=geneid)
## Información de la estructura
class(x)
head(x)

#-------------------------------------------------------------------------------
## obteniendo "Conteos por millón" que servirán para los gráficos exploratorios
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

#-------------------------------------------------------------------------------
# Revisando por genes con conteos cero en todas las muestras 
table(rowSums(x$counts==0)==9)
## Filtrado de genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#-------------------------------------------------------------------------------
## Gráfico comparativos de las densidades de los conteos en el log de los datos
## crudos y en el log de los datos filtrados arriba
lcpm.cutoff <- log2(10/M + 2/L) # log de la media del tamaño de la librería, 
# ajustado para evitar log(0).

library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(x$counts), text.col=col, bty="n")

lcpm <- cpm(x, log=TRUE) #actualizando los log(cpm) usando los datos filtrados
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(x$counts), text.col=col, bty="n")

#-------------------------------------------------------------------------------
# Gráficos de caja para ver la dispersión de las muestras antes y después de
# normalizar
boxplot(lcpm,las=2, col=col)
title(main="A. Unnormalised data",ylab="Log-cpm")
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B.Normalised data",ylab="Log-cpm")

# Nota: Debido a que observamos unos factores de normalización muy cercanos a 1,
# se probó para ver si en realidad el valor esperado es 1

#Primero normalidad
shapiro.test(x$samples$norm.factors) #Nos inclinamos hacia la nula de N(mu,sigma)
#Ahora si la de interés. 
t.test(x$samples$norm.factors, alternative = "two.sided",mu=1)

#-------------------------------------------------------------------------------
# Reducción de dimensionalidad con MDS (método basado en distancias)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")

#-------------------------------------------------------------------------------
# Creando nuestra matriz diseño para el modelo lineal
design <- model.matrix(~0+group)
# Buscando y reemplazando patrones: ¿suena a grep?
colnames(design) <- gsub("group", "", colnames(design))
design

#-------------------------------------------------------------------------------
# Definiendo la matriz de contrastes de interés
contr.matrix <- makeContrasts(
  WTvsult1 = WT-ult1, 
  WTvsult1ult2 = WT - ult1ult2, 
  levels = colnames(design))
contr.matrix

#-------------------------------------------------------------------------------
# Calcular pesos para ajuste de modelo lineal usando limma, al parecer
# es una especie de MC Ponderados
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

#-------------------------------------------------------------------------------
vfit <- lmFit(v, design) # Ajustando modelos lineales a cada gen
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) #Ajustando los contrastes
efit <- eBayes(vfit) # Usando EB para inferencia (en este caso expresión dif)
plotSA(efit, main="Final model: Mean-variance trend") # Mean-variance trend

#-------------------------------------------------------------------------------
# resumen de los resultados de los contrastes
et<-decideTests(efit)
summary(et)

#-------------------------------------------------------------------------------
# usando otro método (también de EB) para expresión dif dado un logfc mínimo
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
str(tfit)
# Extrayendo los índices de los genes que se encuentran dif bajo ambas 
# condiciones 
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$genes[de.common])

#-------------------------------------------------------------------------------
# Diagrama de Venn para mostrar la cantidad de genes en cada condición
par(mfrow=c(1,1))
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

#-------------------------------------------------------------------------------
# Lista (ordenada de menor a mayor según p-value) de genes para los constrastes
# de interés
WT.vs.ult1 <- topTreat(tfit, coef=1, n=Inf)
WT.vs.ult1ult2 <- topTreat(tfit, coef=2, n=Inf)
head(WT.vs.ult1)
head(WT.vs.ult1ult2)

#-------------------------------------------------------------------------------
# Estas son unas gráficas que aún nos faltan por estudiar
par(mfrow=c(1,1))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13)) # Investigar qué hace esta función.
#-------------------------------------------------------------------------------


options(warn = defW)