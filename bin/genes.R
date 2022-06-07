
instala_paquetes<-function(){
  pkgs_CRAN<-c("R.utils", "readr", "dplyr","RColorBrewer")
  pkgs_Bioc<-c("limma", "edgeR")
  chooseCRANmirror(ind = 55)
  #Ahora iteramos por cada paquete
  
  if(!require("BiocManager", quietly = T, character.only = T)){
    writeLines("El paquete BiocManager se va a instalar")
    install.packages("BiocManager")
  }
  if(BiocManager::version()!='3.14') {
    BiocManager::install(version = "3.14", update = F)
  }
  for(pkg in 1:length(pkgs_CRAN)){
    if(!require(pkgs_CRAN[pkg], quietly = T, character.only = T)){
      writeLines("Instalando paquetes necesarios")
      install.packages(pkgs_CRAN[pkg], dependencies = T, quiet = T)
    }
  }
  
  for(pkg in 1:length(pkgs_Bioc)){
    if(!require(pkgs_Bioc[pkg], quietly = T, character.only = T)){
      writeLines("Instalando paquetes necesarios")
      BiocManager::install(pkgs_Bioc[pkg], update = F)
    }
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
writeLines("Selecciona el archivo que contiene los datos de cada muestra escribiendo el numero que le antecede: ")
selFile <- as.numeric(readLines(file("stdin"),1))

cat(c("Archivo",rfiles[selFile]),sep=": ")
writeLines("\n")


countData<- readr::read_delim(file = rfiles[selFile], delim = ",", col_names = T,
                         show_col_types = F)
str(countData)
writeLines("\n")

rm_na<-function(x){
  #Looking for rows (genes) with NAs
  genes_NA<-unique(which(is.na(x))%%dim(x)[1])
  
  #Removing genes with NAs (if there is at least one): 
  if(length(genes_NA)>0){
    #Retrieving and showing names of genes with NA,
    #it is assumed that column 1 stores genes' names:
    writeLines("Los siguientes genes contienen NAs y han sido removidos: ")
    cat(dplyr::pull(x[genes_NA,],1), sep = ", ")
    writeLines("\n")
    x<-na.omit(x)
  }
  return(x)
}



countData<-rm_na(countData)
writeLines("viendo datos después de filtar NAs")
str(countData)
#Una vez limpios de NAs, guardamos la informacion de los ID de genes en un
#vector y dejamos solamente las columna que corresponden a muestras, ie, eliminamos
#la primer columna.
genesID<-countData[,1]
countData<-countData[,-1]
writeLines("Viendo datos después de quitar columna de genes")
str(countData)





library(dplyr)


#-------------------------------------------------------------------------------
#Función para obtener los datos del diseño experimental
getED<- function(validFiles,sampleFileNumber) {  # nolint # nolint
  avaliableFiles<- validFiles[-sampleFileNumber] # nolint
  writeLines("Selecciona el archivo que contiene la informacion del diseño experimental para cada 
  \nmuestra escribiendo el numero que le antecede: ")
  rfilesAv<- paste(paste(1:length(avaliableFiles),".",sep = ""), avaliableFiles, sep = " ")
  cat(rfilesAv, sep = "\n")
  selectedFile <- as.numeric(readLines(file("stdin"),1))
  EDFile<- readr::read_csv(avaliableFiles[selectedFile],show_col_types = F)
  #Validando que no haya muestras con NAs
  sampsWNAs<-unique(which(is.na(EDFile))%%dim(EDFile)[1])
  if(length(sampsWNAs)>0){
    writeLines("Se detectó que las muestras listadas carecen de información para alguna\n
    variable o esta es invalida.\n
    Por favor revisa tu archivo e intentalo de nuevo.")
    cat(dplyr::pull(EDFile[sampsWNAs,],1), sep = ", ")
    quit(save = "no")
  }
  EDFile <- EDFile %>% mutate(across(-1,as.factor))
  return(EDFile)
}
EDMetaData<-getED(rfiles,selFile)


#-------------------------------------------------------------------------------
### Paquetes necesarios
#--------------------
library(limma)      #
library(edgeR)      #
library(R.utils)    #
#--------------------
countData <- DGEList(counts = countData, group = pull(EDMetaData,2), genes = genesID)
## Información de la estructura
#-------------------------------------------------------------------------------
#
writeLines("El analisis esta por comenzar y los resultados seran almacenados
en la carpeta \"Resultados_RNASeq\", en el directorio de trabajo que
usted proporciono al principio.")

folderResultados<- "Resultados_RNASeq"
if(!file.exists(folderResultados)) {
 dir.create(folderResultados)
}



#-------------------------------------------------------------------------------
## obteniendo "Conteos por millón" que servirán para los gráficos exploratorios
cpm <- cpm(countData)
lcpm <- cpm(countData, log=TRUE)
L <- mean(countData$samples$lib.size) * 1e-6
M <- median(countData$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)


library(RColorBrewer)

fileNameGenerator<-function(folder,fileName){
  paste(c(getwd(),folder,fileName), sep = "", collapse = "/")
}

plotGenerator<-function(f,folderOut,fileNameOut,width,height,dpi,...) {
file <- fileNameGenerator(folderOut,fileNameOut) 
png(file = file, width = width, height = height, res = dpi)
f(...)
dev.off()
}

densityPlotCounts<-function(counts) {
  data<-counts
  nsamples <- ncol(data)
  col <- brewer.pal(nsamples, "Paired")
  L <- mean(data$samples$lib.size) * 1e-6
  M <- median(data$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  par(mfrow=c(1,2))
  lcpm <- cpm(data, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="A. Datos Crudos", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", colnames(data$counts), text.col=col, bty="n")
  keep.exprs <- filterByExpr(data, group=data$sample$group)
  data <- data[keep.exprs,, keep.lib.sizes=FALSE]
  lcpm <- cpm(data, log=TRUE) #actualizando los log(cpm) usando los datos filtrados

  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="B. Datos Filtrados", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", colnames(data$counts), text.col=col, bty="n")
}

plotGenerator(densityPlotCounts, folderResultados, "01_densidad_Crudos_vs_Filtrados.png", 1000*2.1, 1800, 300, counts = countData)

#-------------------------------------------------------------------------------
# Revisando por genes con conteos cero en todas las muestras 
table(rowSums(countData$counts==0)==9)
## Filtrado de genes
keep.exprs <- filterByExpr(countData, group=collapse(EDMetaData,2))
countData <- countData[keep.exprs,, keep.lib.sizes=FALSE]
#dim(countData)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Gráficos de caja para ver la dispersión de las muestras antes y después de
# normalizar
boxPlotNormalization<-function(counts){
  data<-counts
  nsamples <- ncol(data)
  col <- brewer.pal(nsamples, "Paired")
  lcpm <- cpm(data, log=TRUE)
  par(mfrow=c(1,2))
  boxplot(lcpm,las=2, col=col)
  title(main="A. Datos crudos",ylab="Log-cpm")
  data <- calcNormFactors(data, method = "TMM")
  data$samples$norm.factors
  lcpm <- cpm(data, log=TRUE)
  boxplot(lcpm, las=2, col=col, main="")
  title(main="B. Datos normalizados",ylab="Log-cpm")
}
plotGenerator(boxPlotNormalization, folderResultados, "02_boxPlot_Crudos_vs_Normalizados.png", 1000*2.1, 1800, 300, counts = countData)

# Nota: Debido a que observamos unos factores de normalización muy cercanos a 1,
# se probó para ver si en realidad el valor esperado es 1

#Primero normalidad
#shapiro.test(x$samples$norm.factors) #Nos inclinamos hacia la nula de N(mu,sigma)
#Ahora si la de interés. 
#t.test(x$samples$norm.factors, alternative = "two.sided",mu=1)

#-------------------------------------------------------------------------------
# Reducción de dimensionalidad con MDS (método basado en distancias)
writeLines("Viendo lo que manda el pull del DE")
collapse(EDMetaData,2)

dimensionReductionPlot<-function(counts) {
  x<-counts
  lcpm <- cpm(x, log=TRUE)
  par(mfrow=c(1,1))
  col.group <- pull(EDMetaData,2)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Paired")
  col.group <- as.character(col.group)
  plotMDS(lcpm, labels=pull(EDMetaData,1), col=col.group)
  title(main="Sample groups")
}

plotGenerator(dimensionReductionPlot,folderResultados, "03_MDS_Dimension_Reduction.png", 1000*2.1, 1800, 300, counts = countData)
#-------------------------------------------------------------------------------

fullSampleInfo <- data.frame(EDMetaData[,-1], countData$samples[,-1], check.names = F)
names(fullSampleInfo) <- c("group", names(fullSampleInfo)[-1])
countData$samples <- fullSampleInfo
# Creando nuestra matriz diseño para el modelo lineal
formulaED <- paste("~0+",paste(names(countData$samples[1:dim(EDMetaData)[2]-1]),sep="",collapse = "+"),sep="") 
formulaED <- as.formula(formulaED)
names(EDMetaData)<-c(names(EDMetaData)[1],names(countData$samples[1:dim(EDMetaData)[2]-1]))
EDMetaData
design <- model.matrix(formulaED, EDMetaData)
# Buscando y reemplazando patrones: ¿suena a grep?
#colnames(design) <- gsub("group", "", colnames(design))s
writeLines("---------------------------------Hasta aqui nos interesa por ahora--------------------------------------------")
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