#############################
# Util functions for data processing
# Author: Tao Lin
#############################
# ------preprocessing
read.origin.pheno_geno <- function(path){
  return (read.cross("csvr", , path, genotypes=c("-1", "1", "0", "9"), alleles=c("-1", "1", "0", "9")))
}

read.gene_expression <- function(path){
  return (read.csv(path, sep="\t", header=TRUE, na.strings=c("","NA")))
}

# ------phewas
load.dataset <- function(tmp.in.path, file.name){
  path.in.file <- paste(tmp.in.path, file.name, sep='/')
  return(readRDS(path.in.file))
}

create.dir <- function(tmp.path, force=TRUE){
  if (dir.exists(tmp.path) && force){
    print(paste("Remove the path: ", tmp.path))
    unlink(tmp.path, recursive=TRUE)
  }
  dir.create(tmp.path, showWarnings=FALSE)
}

remove.file <- function(tmp.path){
  file.remove(tmp.path)
}

extract.filename <- function(filename){
  filename.base <- basename(filename)
  filename.base.list <- strsplit(filename.base, "_")[[1]]
  filename.base.list.length <- length(filename.base.list)
  filename.basename.without.type <- paste(filename.base.list[1: filename.base.list.length - 1], collapse='_')
  filename.basename.with.type <- paste(filename.base.list[1: filename.base.list.length], collapse='_')
  filename.basename.with.type <- strsplit(filename.basename.with.type, "\\.")[[1]][1]
  return(list("basename.without.type"=filename.basename.without.type, "basename.with.type"=filename.basename.with.type))
}

extract.filename.set <- function(tmp.in.path){
  # remove geno and pheno from the filename list, and obtain an unique file name list
  tmp.file.paths <- list.files(tmp.in.path)
  filename.list <- data.frame()
  for(file.name in tmp.file.paths){
    filename.list <- rbind(filename.list, data.frame(extract.filename(file.name)))
  }
  return(filename.list)
}

# ------qtlmapping
pullout.phenotypes.name <- function(in.data){
  name <- names(in.data$pheno)
  len_name <- length(name)
  name <- name[1: len_name - 1]
  return (name)
}
