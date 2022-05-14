#Sys.setlocale(category = "LC_COLLATE", locale = "C")
Sys.setlocale(category = "LC_CTYPE", locale = "C")
#Sys.setlocale(category = "LC_MONETARY", locale = "C")
#Sys.setlocale(category = "LC_NUMERIC", locale = "C")
#Sys.setlocale(category = "LC_TIME", locale = "C")
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
data




options(warn = defW)
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