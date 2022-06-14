
defW<-getOption("warn")
options(warn = -1)
instala_paquetes<-function(){
  pkgs_CRAN<-c("R.utils", "readr", "dplyr","RColorBrewer","stringr")
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
invisible(Sys.setlocale(category = "LC_CTYPE", locale = "C"))
#Sys.setlocale(category = "LC_MONETARY", locale = "C")
#Sys.setlocale(category = "LC_NUMERIC", locale = "C")
#Sys.setlocale(category = "LC_TIME", locale = "C")
invisible(instala_paquetes())
writeLines("
Hola, bienvenido al script para el analisis de RNA-Seq desarrollado por:
                          Alina y Uriel.
")



cat("El directorio de trabajo actual es: ", getwd()," Es el correcto[s/n]: ")

raw <- readLines(file("stdin"),1)
if(raw=="n"){
  cat("Ingresa el directorio correcto: ")
  dr<-raw <- readLines(file("stdin"),1)
  setwd(dr)
}

writeLines("\n")
cat("Los archivos en formato valido que estan disponibles son: \n")
rfiles<-dir()[grep("\\.csv$|\\.txt$",dir())]
rfilesAv<-paste(paste(1:length(rfiles),".",sep = ""), rfiles, sep = " ")
cat(rfilesAv, sep = "\n")
writeLines("Selecciona el archivo que contiene los datos de cada muestra escribiendo el numero que le antecede: ")
selFile <- as.numeric(readLines(file("stdin"),1))

cat(c("Archivo",rfiles[selFile]),sep=": ")
writeLines("\n")


countData<- readr::read_delim(file = rfiles[selFile], delim = ",", col_names = T,
                         show_col_types = F)

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

getED<- function(validFiles,sampleFileNumber) {  # nolint # nolint
  avaliableFiles<- validFiles[-sampleFileNumber] # nolint
  writeLines("Selecciona el archivo que contiene la informacion del disenio experimental para cada 
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
  reflevel <- unique(pull(EDFile,2))[1]
  EDFile <- EDFile %>% mutate(across(-1, factor))
  #print(paste(levels(EDFile$group)),"Aqui deben aparecer los niveles")
  EDFile <- mutate_at(EDFile, 2, relevel, ref = reflevel)
  return(EDFile)
}

EDMetaData<-getED(rfiles,selFile)


countData<-rm_na(countData)


genesID<-countData[,1]
countData<-countData[,-1]





#-------------------------------------------------------------------------------
#Función para obtener los datos del diseño experimental


countData <- DGEList(counts = countData, group = pull(EDMetaData,2), genes = genesID)

writeLines("El analisis esta por comenzar y los resultados seran almacenados
en la carpeta \"Resultados_RNASeq\", en el directorio de trabajo que
usted proporciono al principio.")

folderResultados<- "Resultados_RNASeq"
if(!file.exists(folderResultados)) {
 dir.create(folderResultados)
}



cpm <- cpm(countData)
lcpm <- cpm(countData, log=TRUE)
L <- mean(countData$samples$lib.size) * 1e-6
M <- median(countData$samples$lib.size) * 1e-6



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
  lcpm <- cpm(data, log=TRUE) 

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

keep.exprs <- filterByExpr(countData, group=collapse(EDMetaData,2))
countData <- countData[keep.exprs,, keep.lib.sizes=FALSE]


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

formulaED <- paste("~0+",paste(names(countData$samples[1:dim(EDMetaData)[2]-1]),sep="",collapse = "+"),sep="") 
formulaED <- as.formula(formulaED)
names(EDMetaData)<-c(names(EDMetaData)[1],names(countData$samples[1:dim(EDMetaData)[2]-1]))
design <- model.matrix(formulaED, EDMetaData)



nivelesDisp <- levels(EDMetaData$group)
colnames(design)[1:length(nivelesDisp)] <- nivelesDisp
constrastes <- c()
for(i in 1:(length(nivelesDisp)-1)){
  n<- length(nivelesDisp)
  constrastes<-c(constrastes,
                 paste(nivelesDisp[i],"-", nivelesDisp[-(1:i)]))
}  
contr.matrix <- makeContrasts(
  contrasts = constrastes, 
  levels = colnames(design))


v <- voom(countData, design, save.plot = T)
mean_varPlot_voom <- function(data){
  plot(y~x, data, xlab = "log2(conteo + 0.5)", ylab = "Error estandar", pch = 20,
  main = "voom: Tendencia Media-Varianza")
}
plotGenerator(mean_varPlot_voom,folderResultados, "04_Relacion_Media_Varianza_Pre.png", 1000*2.1, 1800, 300, data = v$voom.xy)


vfit <- lmFit(v, design) 
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 
efit <- eBayes(vfit) 
plotGenerator(plotSA,folderResultados, "05_Relacion_Media_Varianza_MFinal.png", 1000*2.1, 
1800, 300, fit = efit, xlab = "Promedio log-expresion", ylab = "Error estandar",
main="Modelo final: Tendencia Media-Varianza")

et<-decideTests(efit)

extractDown<-function(contraste,geneInfo){
idDown <- which(contraste < 0)
genesDown <- data.frame(pull(geneInfo,1)[idDown])
colnames(genesDown) <- colnames(contraste) 
return(genesDown)
}

extractUp<-function(contraste,geneInfo){
idDown <- which(contraste > 0)
genesDown <- data.frame(pull(geneInfo,1)[idDown])
colnames(genesDown) <- colnames(contraste) 
return(genesDown)
}

saveGenesDE <- function(DE_Genes,up_down_no){
  nombreContraste <- colnames(DE_Genes)
  file <- paste(nombreContraste,"_",up_down_no,"_DEGenes_Nombres.csv",sep = "") 
  fileDir <- fileNameGenerator(folderResultados,fileName = file)
  write.csv(x = DE_Genes, file = fileDir, row.names = F)
}

write.down.DEFiles <- function(constrastes,geneInfo, up_down_no = "Down") {
  nContrastes <- ncol(constrastes)
  for(i in 1:nContrastes){
    down_et <- extractDown(constrastes[,i], geneInfo)
    saveGenesDE(down_et, up_down_no)
  }
}

write.up.DEFiles <- function(constrastes,geneInfo, up_down_no = "Up") {
  nContrastes <- ncol(constrastes)
  for(i in 1:nContrastes){
    down_et <- extractUp(constrastes[,i], geneInfo)
    saveGenesDE(down_et, up_down_no)
  }
}


write.down.DEFiles(et, genesID)
write.up.DEFiles(et, genesID)

write.down.details <- function(contrastes, fit, geneInfo){
  nContrastes <- ncol(contrastes)
  for(i in 1:nContrastes){
    down_et <- extractDown(contrastes[,i], geneInfo)
    details <- topTreat(fit, coef = i, n = Inf, genelist = geneInfo)
    allGenes<- details$gene_id
    commonsID <- match(down_et[,1], allGenes)
    nombreContraste <- colnames(down_et)
    file <- paste(nombreContraste,"_Details_Down","_DEGenes.csv", sep = "")
    fileDir <- fileNameGenerator(folderResultados, fileName = file)
    write.csv(x = details[commonsID,], file = fileDir, row.names = F)
  }
}

write.up.details <- function(contrastes, fit, geneInfo){
  nContrastes <- ncol(contrastes)
  for(i in 1:nContrastes){
    up_et <- extractUp(contrastes[,i], geneInfo)
    details <- topTreat(fit, coef = i, n = Inf, genelist = geneInfo)
    allGenes<- details$gene_id
    commonsID <- match(up_et[,1], allGenes)
    nombreContraste <- colnames(up_et)
    file <- paste(nombreContraste,"_Details_Up","_DEGenes.csv", sep = "")
    fileDir <- fileNameGenerator(folderResultados, fileName = file)
    write.csv(x = details[commonsID,], file = fileDir, row.names = F)
  }
}

write.down.details(et, efit, genesID)
write.up.details(et, efit, genesID)

##-------------------------------------------------------------------------------
## Estas son unas gráficas que aún nos faltan por estudiar
#par(mfrow=c(1,1))
#plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
#       xlim=c(-8,13)) # Investigar qué hace esta función.
##-------------------------------------------------------------------------------
#
#
#options(warn = defW)