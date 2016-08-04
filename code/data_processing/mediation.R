#############################
# mediation
# Author: Tao Lin
#############################
# import packages
library(qtl)
library(stringr)
library(data.table)
library(intermediate)
library(parallel)
library(foreach)
library(doParallel)

############################# Functions #############################
# prepare the mediation
prepare.single.mediation.get.prob <- function(v){
  # -1 - BB, 1 - DD, 0 - BD, 9 - NA
  if (!is.numeric(v)){
    return(array(data = c(0.0, 0.0)))
  }
  if (is.na(v)) {
    return(array(data = c(0.0, 0.0)))
  }
  if ( v == -1) { # BB
    return(array(data = c(1, 0)))
  }
  if (v == 1){ # DD
    return(array(data = c(0, 1)))
  }
  if (v == 0){ # BD
    return(array(data = c(0.5, 0.5)))
  }
  if (v == 9){ # NA
    return(array(data = c(0.0, 0.0)))
  }
  return(array(data = c(0.0, 0.0)))
}

prepare.single.mediation.build.prob_matrix <- function(filename, geno){
  # build a empty 3d matrix
  dim.first <- dim(geno)[1]
  dim.second <- dim(geno)[2]
  dim.third <- 2
  name.dim.first <- rownames(geno)
  name.dim.second <- colnames(geno)
  name.dim.third <- c("B", "D")
  prob_matrix <- array(NA, c(dim.first, dim.second, dim.third), dimnames = list(name.dim.first, name.dim.second, name.dim.third))

  # init the cluster
  # num.cores <- detectCores() - 2
  # cl <- makeCluster(num.cores, outfile="")
  # registerDoParallel(cl)

  # fill the prob to the prob_matrix
  # foreach(ind.row=1: dim.first, .export="prepare.single.mediation.get.prob") %dopar% {
  message(paste("Build probability matrix: ", filename))
  
  for(ind.row in seq(1, dim.first)) {
    for(ind.col in seq(1, dim.second)) {
      prob <- prepare.single.mediation.get.prob(geno[ind.row, ind.col])
      prob_matrix[ind.row, ind.col, 1:2] <- prob
      # print(paste(filename, "ind.row", ind.row, "ind.col", ind.col))
    }
  }
  # stopCluster(cl)
  return(prob_matrix)
}

prepare.get.info.file <- function(in.path){
  splitted_path  <- strsplit(in.path, split = '/')[[1]]
  len_path       <- length(splitted_path)
  # protein
  name.protein   <- splitted_path[len_path - 1]
  name.protein.list <- strsplit(name.protein, "_")[[1]]
  name.protein.length <- length(name.protein.list)
  name.protein <- paste(name.protein.list[2: name.protein.length], collapse='_')
  # gene
  name.gene      <- splitted_path[len_path]
  name.gene      <- extract.filename(name.gene)$basename.with.type
  return(data.frame("symbol"=name.gene, "protein_name"=name.protein))
}

prepare.get.info.peak <- function(in.data){
  info.peak <- in.data[which(max(in.data$lod) == in.data$lod), ]
  return(info.peak)
}

prepare.get.bestqtl <- function(in.path){
  # load qtl info dataset
  qtl.info   <- readRDS(in.path)
  # get the info.
  info.file <- prepare.get.info.file(in.path)
  info.peak <- prepare.get.info.peak(qtl.info)
  return(data.table("protein_name"=info.file$protein_name, "symbol"=info.file$symbol, "snp"=rownames(info.peak), "peak_chr"=info.peak$chr, "peak_pos"=info.peak$pos, "peak_lod"=info.peak$lod))
}

prepare.single.mediation.get.bestqtl <- function(filename, tmp.in.path.qtl){
  # get whole list of qtl files
  qtl.filelist <- list.files(tmp.in.path.qtl)
  qtl.filelist.length <- length(qtl.filelist)
  # do the statistics
  err_count <- 0

  info <- data.table()
  for (ind in seq(1, qtl.filelist.length)){
    # get file name and file path.
    qtlfile <- qtl.filelist[ind]
    qtlfile.path <- paste(tmp.in.path.qtl, qtlfile, sep="/")
    message(paste0("Processing ", ind, "'th (", ind / qtl.filelist.length * 100 , "%), " , qtlfile, " in ", filename))

    bestqtl <- tryCatch({
      prepare.get.bestqtl(qtlfile.path)
    }, error=function(err){
      message(paste("ERROR: ", err))
      err_count <- err_count + 1
      return(data.frame())
    })
    
    if (length(info) == 0) info <- bestqtl
    else info <- rbind(info, bestqtl)
  }
  
  message(paste("......Existing", err_count, "errors......"))
  
  return(info)
}

prepare.all.mediation <- function(filename.set, tmp.in.path.geno_pheno, tmp.in.qtlmapping, tmp.out.path){
  # define index
  tmp.file.from <- 1
  tmp.file.to <- length(filename.set)

  # config parallel setting
  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores, outfile="")
  registerDoParallel(cl)
  
  foreach(tmp.index=tmp.file.from: tmp.file.to, .export = c("load.dataset", "create.dir", "data.table",
                                                            "prepare.single.mediation.build.prob_matrix",
                                                            "prepare.single.mediation.get.prob",
                                                            "prepare.single.mediation.get.bestqtl",
                                                            "prepare.get.bestqtl", "prepare.get.info.file",
                                                            "extract.filename", "prepare.get.info.peak")) %dopar% {
  # for (tmp.index in seq(tmp.file.from, tmp.file.to)){
    # define file and path
    filename <- filename.set[tmp.index]
    path.in.geno <- paste0(filename, "_geno.rds")
    path.in.qtlmapping <- paste(tmp.in.qtlmapping, filename, sep="/")
    path.out <- paste(tmp.out.path, filename, sep="/")
    create.dir(path.out, force=FALSE)
    message(paste0("Create dir ", path.out, " ..."))

    path.out.geno <- paste(path.out, "genomatrix.rds", sep="/")
    path.out.bestqtl <- paste(path.out, "bestqtl.rds", sep="/")

    # load dataset
    message(paste("Preparing organism ", filename, "for the mediation..."))
    geno <- load.dataset(tmp.in.path.geno_pheno, path.in.geno)

    # get qtl information
    if (!file.exists(path.out.bestqtl)) {
      bestqtl <- prepare.single.mediation.get.bestqtl(filename, path.in.qtlmapping)
      saveRDS(bestqtl, path.out.bestqtl)
    }

    # using geno to build a 3d matrix
    if (!file.exists(path.out.geno)) {
    prob_matrix <- prepare.single.mediation.build.prob_matrix(filename, geno)
    saveRDS(prob_matrix, path.out.geno)
    }
  }
  stopCluster(cl)
}

# do the mediation
extract.ordered.pheno.info <- function(info.pheno, info.bestqtl.clean){
  symbol.names <- info.bestqtl.clean$symbol
  ordered.pheno.info <- data.frame()
  
  for (name in symbol.names){
    tmp <- info.pheno[name] 
    if (length(ordered.pheno.info) == 0) ordered.pheno.info <- tmp
    else ordered.pheno.info <- cbind(ordered.pheno.info, tmp)
  }
  return(ordered.pheno.info)
}
build.vaild.info <- function(info.bestqtl, info.gene_position, info.pheno){
  # rebuild gene position
  colnames(info.gene_position) <- c("cleansymbol", "chr", "pos")
  info.gene_position <- data.table(info.gene_position)
  info.gene_position <- info.gene_position[!duplicated(info.gene_position$cleansymbol), ]
  # rebuild gene, i.e., remove useless information from gene symbol, and add it as a new column.
  info.bestqtl <- data.table(info.bestqtl)
  info.bestqtl.clean_list <- lapply(info.bestqtl[, symbol], function(x){
    strsplit(as.character(x), "_")[[1]][1]
  })
  info.bestqtl <- info.bestqtl[, cleansymbol:=unlist(info.bestqtl.clean_list)]
  # join based on the same symbol, and built annotation.
  annot <- info.gene_position[info.bestqtl, on=c("cleansymbol"), nomatch=0]
  # rebuild bestqtl
  gene_list <- data.table("symbol"=unique(annot$symbol))
  info.bestqtl.clean <- info.bestqtl[gene_list, on="symbol", allow.cartesian=TRUE, nomatch=0]
  # extract the protein
  info.valid.pheno.ordered <- extract.ordered.pheno.info(info.pheno, info.bestqtl.clean)
  # build covar
  info.faked_covar <- matrix(1, nrow=dim(info.valid.pheno.ordered)[1], ncol=1)
  return(list("pheno"=info.valid.pheno.ordered, "faked_covar"=info.faked_covar, "annot"=annot, "bestqtl"=info.bestqtl.clean))
}
do.single.mediation <- function(tmp.out.path, filename, tmp.in.path.geno_pheno, tmp.in.path.qtl_mediation, info.gene_position){
  # define path and load dataset
  path.qtl_mediation <- paste(tmp.in.path.qtl_mediation, filename, sep="/")
  info.pheno <- load.dataset(tmp.in.path.geno_pheno, paste0(filename, "_pheno.rds"))
  info.prob_matrix <- load.dataset(path.qtl_mediation, "genomatrix.rds")
  info.bestqtl <- load.dataset(path.qtl_mediation, "bestqtl.rds")
  path.output <- file.path(tmp.out.path, filename)
  create.dir(path.output, force=FALSE)

  # Get valid information.
  info.valid <- build.vaild.info(info.bestqtl, info.gene_position, info.pheno)
  info.valid.pheno.matrix <- matrix(as.numeric(unlist(info.valid$pheno)), nrow=nrow(info.valid$pheno))
  info.valid.faked_covar <- info.valid$faked_covar
  # run mediation
  err_count <- 0
  tmp.index.from <- 1
  tmp.index.to <- length(info.valid$pheno)
  # config parallel setting
  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores, outfile="")
  registerDoParallel(cl)
  
  foreach(tmp.index=tmp.index.from: tmp.index.to, .export = c("mediation.scan", "data.table")) %dopar% {
    # define path
    info.valid.pheno.ind <- info.valid$pheno[tmp.index]
    path.output.result <- file.path(path.output, paste0(colnames(info.valid.pheno.ind), ".rds"))
    # processing
    if (!file.exists(path.output.result)){
      info.bestqtl.gene <- info.valid$bestqtl[symbol == colnames(info.valid.pheno.ind)]
      
      tryCatch({
        med <- mediation.scan(target = as.numeric(as.matrix(info.valid.pheno.ind)), mediator = info.valid.pheno.matrix, annotation = info.valid$annot,
                       qtl.geno = info.prob_matrix[, info.bestqtl.gene$snp, ], covar = info.valid.faked_covar, verbose = FALSE)
        saveRDS(med, path.output.result)
      }, error=function(err){
        message(paste("ERROR: ", err))
        err_count <- err_count + 1
      })
    }
    message(paste("Do the single mediation for", filename, "current index:", tmp.index, ",total:", tmp.index.to))
  }
  stopCluster(cl)
}
do.all.mediation <- function(filename.set, tmp.in.path.geno_pheno, tmp.in.path.qtl_mediation, tmp.in.path.global, tmp.out.path){
  # define path index
  tmp.file.from <- 1
  tmp.file.to <- length(filename.set)

  # load the dataset
  gene_position_aligner <- load.dataset(tmp.in.path.global, "gene_position_aligner.rds")

  for(tmp.index in seq(tmp.file.from, tmp.file.to)){
    filename <- filename.set[tmp.index]
    message(paste0("Do mediation on ", filename, "..."))
    do.single.mediation(tmp.out.path, filename, tmp.in.path.geno_pheno, tmp.in.path.qtl_mediation, gene_position_aligner)
  }
}

############################# Main Entry #############################
# ------Define paths
path.root <- "/data/Dropbox/sv/gene_expression/new"
path.root <- "/Users/hali/Dropbox/sv/gene_expression/new"
path.root <- "/home/gene_expression"

path.code <- paste(path.root, "code", sep="/")
path.code.data_processing <- paste(path.code, "data_processing", sep="/")

path.data <- paste(path.root, "data", sep="/")
path.data.in <-paste(path.data, "input", sep="/")
path.data.out <- paste(path.data, "output", sep="/")

path.data.in.global <- paste(path.data.out, "global", sep="/")
path.data.in.geno_pheno <- paste(path.data.out, "geno_pheno_splitted_rds", sep="/")
path.data.in.qtlmapping <- paste(path.data.out, "qtlmapping", sep="/")
path.data.out.global <- paste(path.data.out, "global", sep="/")
path.data.out.qtl_mediation <- paste(path.data.out, "qtlmapping_mediation", sep="/")
path.data.out.mediation <- paste(path.data.out, "mediation", sep="/")

# ------Load functions/parameters from other files
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------Parse the filename set from the filename list
filename.set <- unique(extract.filename.set(path.data.in.geno_pheno)$basename.without.type)

filename.set.new <- list()
for (ind in seq(1, length(filename.set))){
  name <- filename.set[ind]
  if (grepl("Metabolite", name)) next
  filename.set.new[[length(filename.set.new) + 1]] <- name
}
filename.set <- unlist(filename.set.new)

# ------Debug
# filename = filename.set[4]
# tmp.in.path.geno_pheno = path.data.in.geno_pheno
# tmp.in.path.qtl_mediation = path.data.out.qtl_mediation
# tmp.in.path.global = path.data.in.global
# tmp.out.path = path.data.out.mediation
# info.gene_position <- load.dataset(tmp.in.path.global, "gene_position_aligner.rds")
# tmp.index = 1

# ------Run on Server
# prepare.all.mediation(filename.set, path.data.in.geno_pheno, path.data.in.qtlmapping, path.data.out.qtl_mediation)
do.all.mediation(filename.set, path.data.in.geno_pheno, path.data.out.qtl_mediation, path.data.in.global, path.data.out.mediation)
