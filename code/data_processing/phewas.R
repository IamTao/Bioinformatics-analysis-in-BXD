#############################
# Phewas
# Author: Tao Lin
#############################
# import packages
library(qtl)
library(emma)
library(plyr)
library(biomaRt)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)

############################# Functions #############################
# Calculate p value
calculate.single.pvalue.assistant <- function(pheno, geno, ind, pval_all, pheno.length){
  message(paste0("Processing ", ind, "-th pheno, existing ", pheno.length, " phenotypes (", ind / pheno.length * 100, "%)"))

  pheno.at_index <- pheno[, ind]
  pheno.at_index.not_na <- !is.na(pheno.at_index)
  pheno.at_index <- pheno.at_index[pheno.at_index.not_na]
  geno.at_index  <- geno[pheno.at_index.not_na, ]

  if(length(pheno.at_index) < 15){
    pheno.name <- colnames(pheno)[ind]
    message(paste0("phenotype ID ", pheno.name, " has ", length(pheno.at_index), " strains, not enough for analysis!"))
  } else{
    geno_imp.at_index <- as.matrix(geno.at_index)
    average.id <- colMeans(geno_imp.at_index, na.rm = T)

    for (j in 1: ncol(geno_imp.at_index)){
      geno_imp.at_index[is.na(geno_imp.at_index[, j]), j] <- average.id[j]
    }

    stdev <- apply(geno_imp.at_index, 2, sd)
    stdev[which(stdev == 0)] <- 1 # used to fix the bug when stdev is 0, because the remaining strains have the same genotype
    geno_stand.at_index <- sweep(sweep(geno_imp.at_index, 2, average.id, "-"), 2, stdev, "/")
    K_mat.i <- (geno_stand.at_index %*% t(geno_stand.at_index)) / ncol(geno_stand.at_index)

    pheno.at_index <- as.numeric(as.character(pheno.at_index))
    mygwas <- mlmm(Y = pheno.at_index, X = geno_imp.at_index, K = K_mat.i, nbchunks = 2, maxsteps = 3)

    pval <- mygwas$pval_step[[1]]$out
    colnames(pval)[2] <- ind
    if (length(pval_all) == 0){
      pval_all <- pval
    }
    else{
      pval_all <- join(pval_all, pval, by = "SNP", type = "left", match = "all")
    }
  }
  return(pval_all)
}

calculate.single.pvalue <- function(tmp.in.path, filename, tmp.out.path, bxd.mapping){
  path.pheno <- paste0(filename, "_pheno.rds")
  path.geno <- paste0(filename, "_geno.rds")
  pheno <- load.dataset(tmp.in.path, path.pheno)
  pheno.length <- dim(pheno)[2]
  geno <- load.dataset(tmp.in.path, path.geno)

  pval_all <- data.frame()
  err_count <- 0
  for (ind in 1: pheno.length){
    pval_all <- tryCatch({
      calculate.single.pvalue.assistant(pheno, geno, ind, pval_all, pheno.length)
    }, error=function(err){
      message(paste("ERROR: ", err))
      err_count <- err_count + 1
      return(pval_all)
    })
  }
  message(paste("......Existing", err_count, "errors......"))

  pval_all <- data.table(pval_all)
  path.out.workspace <- paste(tmp.out.path, paste0(filename, ".rds"), sep="/")
  saveRDS(pval_all, file=path.out.workspace)
}

calculate.all.pvalue.parallel <- function(tmp.in.path.code, tmp.in.path.data, tmp.in.path.global, filename.set, tmp.out.path, force=FALSE){
  # Load global dataset
  bxd.mapping <- load.dataset(tmp.in.path.global, file.name="map_BXD.rds")

  tmp.file.from <- 1
  tmp.file.to <- length(filename.set)

  # start parallel
  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores, outfile="")
  registerDoParallel(cl)
  foreach(tmp.index=tmp.file.from: tmp.file.to, .export=c("calculate.single.pvalue", "load.dataset",
                                                          "calculate.single.pvalue.assistant", "join", "data.table")) %dopar% {
    source(paste(tmp.in.path.code, "emma.r", sep="/"))
    source(paste(tmp.in.path.code, "mlmm.r", sep="/"))
    
    filename <- filename.set[tmp.index]                         
    path.out.workspace <- paste(tmp.out.path, paste0(filename, ".rds"), sep="/")
    if (force || !file.exists(path.out.workspace)) {
      message(paste("Processing", filename, "..."))
      calculate.single.pvalue(tmp.in.path.data, filename, tmp.out.path, bxd.mapping)
    }
    else{
      message(paste("Existing processed", filename, "continue..."))
    }
  }
  stopCluster(cl)
}

calculate.all.pvalue.sequential <- function(tmp.in.path.data, tmp.in.path.global, filename.set, tmp.out.path, force=FALSE){
  # Load global dataset
  bxd.mapping <- load.dataset(tmp.in.path.global, file.name="map_BXD.rds")

  tmp.file.from <- 1
  tmp.file.to <- length(filename.set)

  for(tmp.index in seq(tmp.file.from: tmp.file.to)) {
    filename <- filename.set[tmp.index]                         
    path.out.workspace <- paste(tmp.out.path, paste0(filename, ".rds"), sep="/")
    if (force || !file.exists(path.out.workspace)) {
      message(paste("Processing", filename, "..."))
      calculate.single.pvalue(tmp.in.path.data, filename, tmp.out.path, bxd.mapping)
    }
    else{
      message(paste("Existing processed", filename, "continue..."))
    }
  }
}

# Calculate phewas
calculate.single.phewas.assistant <- function(gene.info, bxd.mapping, phenotypes.id_aligner, pvalue){
  # get/define the snp info of the gene (nearest for now)
  region   <- 2
  N_eff <- 2000
  significantline_bonf<- -log10(0.05/dim(pvalue)[2])
  significantline_eff <- -log10(0.05/N_eff)
  rotation <- 40

  # get gene information
  gene.name <- gene.info$GenaName
  gene.chromosome_name <- gene.info$chromosome_name
  gene.start_postition <- gene.info$start_position
  gene.end_postition <- gene.info$end_position

  # snps in a near region
  snp.info <- bxd.mapping[Chr == gene.chromosome_name][Pos >= (gene.start_postition - region)][Pos <= (gene.start_postition + region)]
  # snps within the gene
  # snp.info <- bxd.mapping[Chr == gene.chromosome_name][Pos >= min(gene.start_postition, gene.end_postition)][Pos <= max(gene.start_postition, gene.end_postition)]

  # calculate y
  filterdata <- pvalue[snp.info, allow.cartesian=TRUE, on="SNP", nomatch=0]

  set(filterdata,j=c('SNP', 'Pos', 'Chr'), value=NULL)
  pval  <- apply(filterdata, 2, min)
  pfinal <- data.table(pval)
  pfinal[, PhenoID.1:=seq(1, length(pval))]
  pfinal[, logp:=-log10(pval)]
  y <-  phenotypes.id_aligner[pfinal, on="PhenoID.1", nomatch=0]
  return(y)
}

calculate.single.phewas <- function(tmp.in.path, filename, tmp.out.path, bxd.mapping, phenotypes.id_aligner, gene.position_aligner, force){
  # define path or create path
  tmp.out.parentpath <- paste(tmp.out.path, filename, sep="/")
  filename.withextension <- paste0(filename, ".rds")
  create.dir(tmp.out.parentpath, force)

  # load dataset
  pvalue <- load.dataset(tmp.in.path, file.name=filename.withextension)

  # define index
  tmp.index.from <- 1
  tmp.index.to <- nrow(gene.position_aligner)

  # start parallel
  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores, outfile="")
  registerDoParallel(cl)

  foreach(tmp.index=tmp.index.from: tmp.index.to, .export=c("calculate.single.phewas.assistant", "data.table", "set")) %dopar% {
    # define path
    gene.info <- gene.position_aligner[tmp.index, ]
    filename.withextension <- paste0(gene.info$GeneName, ".rds")
    tmp.out.path <- paste(tmp.out.parentpath, filename.withextension, sep="/")
    # processing
    if (force || !file.exists(tmp.out.path)){
      message(paste("Processing", filename, "for", tmp.index, "-th gene"))
      calculated.y_info <- calculate.single.phewas.assistant(gene.info, bxd.mapping, phenotypes.id_aligner, pvalue)
      saveRDS(calculated.y_info, file=tmp.out.path)
    }
    else{
      message(paste("Exsiting", filename, "for", tmp.index, "-th gene. Ignore..."))
    }

  }
  stopCluster(cl)
}

calculate.all.phewas <- function(tmp.in.path.data, tmp.in.path.global, filename.set, tmp.out.path, force=FALSE){
  # Load global dataset
  phenotypes.id_aligner <- data.table(load.dataset(tmp.in.path.global, file.name="phenotypes_id_aligner.rds"))
  bxd.mapping <- data.table(load.dataset(tmp.in.path.global, file.name="map_BXD.rds"))
  gene.position_aligner <- load.dataset(tmp.in.path.global, file.name="gene_position_aligner.rds")

  tmp.file.from <- 1
  tmp.file.to <-  length(filename.set)

  for(tmp.index in seq(tmp.file.from: tmp.file.to)) {
    calculate.single.phewas(tmp.in.path.data, filename.set[tmp.index], tmp.out.path, bxd.mapping, phenotypes.id_aligner, gene.position_aligner, force)
  }
}
############################# Main Entry #############################
# ------Define paths
path.root <- "/data/Dropbox/sv/gene_expression/new"
path.root <- "/Users/hali/Dropbox/sv/gene_expression/new"
# path.root <- "/home/gene_expression"
path.code <- paste(path.root, "code", sep="/")
path.code.data_processing <- paste(path.code, "data_processing", sep="/")

path.data <- paste(path.root, "data", sep="/")
path.data.in <-paste(path.data, "input", sep="/")
path.data.out <- paste(path.data, "output", sep="/")

path.data.in.global <- paste(path.data.out, "global", sep="/")
path.data.in.geno_pheno <- paste(path.data.out, "geno_pheno_splitted_rds", sep="/")
path.data.out.global <- paste(path.data.out, "global", sep="/")
path.data.out.pvalue <- paste(path.data.out, "pvalue", sep="/")
path.data.out.phewas <- paste(path.data.out, "phewas", sep="/")

# ------Load functions/parameters from other files
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------Load functions from other files
source(paste(path.code.data_processing, "emma.r", sep="/"))
source(paste(path.code.data_processing, "mlmm.r", sep="/"))
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------Parse the filename set from the filename list
filename.set <- unique(extract.filename.set(path.data.in.geno_pheno)$basename.without.type)

# ------Debug
#path.pheno <- paste0(filename.set[[4]], "_pheno.rds")
#path.geno <- paste0(filename.set[[4]], "_geno.rds")
#pheno <- load.dataset(path.data.in.geno_pheno, path.pheno)
#geno <- load.dataset(path.data.in.geno_pheno, path.geno)

# calculate.single.pvalue(path.data.in.geno_pheno, filename.set[[4]], path.data.out.pvalue)
# calculate.single.phewas(path.data.out.pvalue, filename.set[[4]], path.data.out.phewas, bxd.mapping, phenotypes.id_aligner, gene.position_aligner)

# ------Run on Server
# calculate.all.pvalue.sequential(path.data.in.geno_pheno, path.data.in.global, filename.set, path.data.out.pvalue, force=FALSE)
calculate.all.pvalue.parallel(path.code.data_processing, path.data.in.geno_pheno, path.data.in.global, filename.set, path.data.out.pvalue, force=FALSE)
# calculate.all.phewas(path.data.out.pvalue, path.data.out.global, filename.set, path.data.out.phewas, force=FALSE)
