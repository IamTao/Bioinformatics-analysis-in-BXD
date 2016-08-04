#############################
# Preprocessing the dataset for the subseqent analysis(e.g., phewas, qtl, mediation)
# Author: Tao Lin
#############################
# import packages
library(data.table)
library(qtl)
library(parallel)
library(foreach)
library(doParallel)

############################# Functions #############################
# Deal with genotype_expression
generate.rds.output.splitted <- function(tmp.out.path, filename){
  filename <- extract.filename(filename)
  basename <- filename$basename.with.type
  path.out.tmpfile <- paste(tmp.out.path, basename, sep="/")
  return(paste0(path.out.tmpfile, ".rds"))
}

deal.single.genotype_expression.splitted <- function(tmp.in.path, filename, tmp.out.path){
  path.in.file <- paste(tmp.in.path, filename, sep="/")
  message(paste0("Processing file from path ", path.in.file))
  path.out.file <- generate.rds.output.splitted(tmp.out.path, filename)
  message(paste0("Generate file to path ", path.out.file))

  in.data <- read.gene_expression(path.in.file)
  in.data.name <- in.data[, 1]
  rownames(in.data) <- in.data.name
  in.data[, 1] <- NULL
  in.data.tranpose <- as.data.frame(t(in.data))
  saveRDS(in.data.tranpose, file=path.out.file)
  return(in.data.tranpose)
}

filter.invalid.data.original <- function(in.data){
  pheno <- in.data$pheno
  pheno.colname <- colnames(pheno)
  
  # remove 'affy_' gene
  matching.condition <- lapply(seq(1, length(pheno.colname)), function(ind){
    name <- pheno.colname[ind]
    ismatch <- grepl("affy_", tolower(name))
  })
  new.pheno <- pheno[!unlist(matching.condition)]
  # While use read.cross to read the data, the 'Krtap9-3_A_51_P100624' will become 'Krtap9.3_A_51_P100624'. The following code fix this issue.
  new.pheno.colname <- lapply(colnames(new.pheno), function(name){
    if (grepl("\\.", name)) name <- gsub("\\.", "-", name)
    return (name)
  })
  new.pheno.colname <- unlist(new.pheno.colname)
  colnames(new.pheno) <- new.pheno.colname
  in.data$pheno <- new.pheno
  return(in.data)
}

generate.rds.output.original <- function(tmp.out.path, filename){
  basename <- extract.filename(filename)$basename.with.type
  path.out.tmpfile <- paste(tmp.out.path, basename, sep="/")
  return(paste0(path.out.tmpfile, ".rds"))
}

deal.single.genotype_expression.original <- function(tmp.in.path, filename, tmp.out.path){
  path.in.file <- paste(tmp.in.path, filename, sep="/")
  message(paste0("Processing file from path ", path.in.file))
  path.out.file <- generate.rds.output.original(tmp.out.path, filename)
  message(paste0("Generate file to path ", path.out.file))

  in.data <- read.origin.pheno_geno(path.in.file)
  in.data <- filter.invalid.data.original(in.data)
  saveRDS(in.data, file=path.out.file)
}

deal.all.genotype_expression <- function(tmp.in.geno_pheno, tmp.out.geno_pheno, tmp.in.geno_pheno_origin, tmp.out.geno_pheno_origin){
  # define path
  tmp.file.paths <- list.files(tmp.in.geno_pheno)
  tmp.file.from <- 1
  tmp.file.to <- length(tmp.file.paths)
  
  # begin parallel
  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  foreach(tmp.index=tmp.file.from: tmp.file.to, .export=c("deal.single.genotype_expression.splitted", "read.gene_expression", "generate.rds.output.splitted", "extract.filename")) %dopar% {
    message(paste0("Processing file from path ", tmp.file.paths[tmp.index]))
    deal.single.genotype_expression.splitted(tmp.in.geno_pheno, tmp.file.paths[tmp.index], tmp.out.geno_pheno)
  }

  tmp.file.paths <- list.files(tmp.in.geno_pheno_origin)
  tmp.file.from <- 1
  tmp.file.to <- length(tmp.file.paths)
  foreach(tmp.index=tmp.file.from: tmp.file.to, .export=c("deal.single.genotype_expression.original", "data.table", "read.origin.pheno_geno",
                                                          "generate.rds.output.original", "read.cross", "extract.filename", "filter.invalid.data.original")) %dopar% {
    message(paste0("Processing file from path ", tmp.file.paths[tmp.index]))
    deal.single.genotype_expression.original(tmp.in.geno_pheno_origin, tmp.file.paths[tmp.index], tmp.out.geno_pheno_origin)
  }
  stopCluster(cl)
}

# Deal with normal csv file
deal.mapping <- function(tmp.in.path, file.name, tmp.out.path){
  path.in.file <- paste(tmp.in.path, file.name, sep="/")
  message(paste0("Processing file from path ", path.in.file))
  path.out.file <- generate.rds.output.original(tmp.out.path, file.name)

  in.data <- read.table(path.in.file, sep = "\t", header = T)
  saveRDS(in.data, file=path.out.file)
}

# Deal with phenotypes_id_aligner
deal.phenotypes_id_aligner <- function(tmp.in.path, file.name, tmp.out.path){
  path.in.file <- paste(tmp.in.path, file.name, sep="/")
  message(paste0("Processing file from path ", path.in.file))
  path.out.file <- generate.rds.output.original(tmp.out.path, file.name)

  in.data <- read.csv(path.in.file, sep = ",", header = T)
  in.data.todatatable <- data.table(in.data)
  saveRDS(in.data, file=path.out.file)
}

# Deal with gene_aligner
deal.gene_aligner <- function(tmp.in.path, file.name, tmp.out.path){
  path.in.file <- paste(tmp.in.path, file.name, sep="/")
  message(paste0("Processing file from path ", path.in.file))
  path.out.file <- generate.rds.output.original(tmp.out.path, file.name)

  in.data <- read.csv(path.in.file, sep = ",", header = T)
  in.data <- in.data[-seq(1, 27, by=1), ]
  in.data.todatatable <- data.table(in.data)
  saveRDS(in.data, file=path.out.file)
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
path.data.in.global <- paste(path.data.in, "global", sep="/")
path.data.in.geno_pheno_origin <- paste(path.data.in, "geno_pheno", sep="/")
path.data.in.geno_pheno <- paste(path.data.out, "geno_pheno_splitted_csv", sep="/")

path.data.out.geno_pheno <- paste(path.data.out, "geno_pheno_splitted_rds", sep="/")
path.data.out.geno_pheno_origin <- paste(path.data.out, "geno_pheno_origin_rds", sep="/")
path.data.out.global <- paste(path.data.out, "global", sep="/")
path.data.out.pvalue <- paste(path.data.out, "pvalue", sep="/")
path.data.out.phewas <- paste(path.data.out, "phewas", sep="/")
path.data.out.qtlmapping.rds <- paste(path.data.out, "qtlmapping", sep="/")
path.data.out.qtlmapping.summary <- paste(path.data.out, "qtlmapping_summary", sep="/")
path.data.out.qtl_mediation <- paste(path.data.out, "qtlmapping_mediation", sep="/")
path.data.out.mediation <- paste(path.data.out, "mediation", sep="/")

# ------Load functions/parameters from other files
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------ create some necessary paths beforehand
create.dir(path.data.out.geno_pheno)
create.dir(path.data.out.geno_pheno_origin)
create.dir(path.data.out.global)
create.dir(path.data.out.pvalue)
create.dir(path.data.out.phewas)
create.dir(path.data.out.qtlmapping.rds)
create.dir(path.data.out.qtlmapping.summary)
create.dir(path.data.out.qtl_mediation)
create.dir(path.data.out.mediation)

# ------Debug
tmp.in.geno_pheno <-  path.data.in.geno_pheno
tmp.out.geno_pheno  <- path.data.out.geno_pheno
tmp.in.geno_pheno_origin <- path.data.in.geno_pheno_origin
tmp.out.geno_pheno_origin <- path.data.out.geno_pheno_origin
# files <- list.files(path.data.in.geno_pheno)
# t <- deal.single.genotype_expression.splitted(path.data.in.geno_pheno, files[10], path.data.out.geno_pheno) # focus on eye.
# files <- list.files(path.data.in.geno_pheno_origin)
# t <- deal.single.genotype_expression.original(path.data.in.geno_pheno_origin, files[10], path.data.out.geno_pheno_origin) # focus on eye.

# ------ Run on Server
deal.all.genotype_expression(path.data.in.geno_pheno, path.data.out.geno_pheno, path.data.in.geno_pheno_origin, path.data.out.geno_pheno_origin)
deal.mapping(path.data.in.global, "map_BXD.txt", path.data.out.global)
deal.phenotypes_id_aligner(path.data.in.global, "phenotypes_id_aligner.csv", path.data.out.global)
deal.gene_aligner(path.data.in.global, "gene_position_aligner.csv", path.data.out.global)