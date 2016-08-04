#############################
# qtlmapping
# Author: Tao Lin
#############################
# import packages
library(qtl)
library(stringr)
library(parallel)
library(foreach)
library(doParallel)

############################# Functions #############################
# qtl mapping
calculate.single.qtlmapping <- function(filename, phenotypes.name, ind, in.data, tmp.path.out){
  phenotype.name <- phenotypes.name[ind]
  phenotype.name.withextension <- paste0(phenotype.name, ".rds")
  path.out.rds <- paste(tmp.path.out, phenotype.name.withextension, sep="/")
  if (!file.exists(path.out.rds)) {
    message(paste0("Processing: ", filename, ": ", ind, "'th pheno; its name: ", phenotype.name))
    # Skips traits that have fewer than 8 strains' worth of data (minimum for QTL mapping)
    #if(length(which(!is.na(in.data$pheno[ , ind]))) < 8) {
    #  chr <- c("Insufficient data for QTL calculations")
    #  y <- cbind(phenotypes.name[ind], chr)
    #  return(y)
    #}
  
    # Calculates the QTL with NORMAL ASSUMPTION
    out.i <- scanone(in.data, pheno.col=ind, method="hk", model="normal")
    # Calculates the significance threshold for the QTL
    # operm <- scanone(in.data, pheno.col=ind, method="hk", n.perm=1000)
    # added to correct for the error caused by infinite lod score after scanone
    # out.i$lod[is.infinite(out.i$lod)] = NA
    # Pulls out all markers with p values < alpha
    # x <- summary(out.i, perms=operm, alpha=0.67, pvalues=TRUE)
  
    # If no results with alpha < the number selected above, initiate blank results; REQUIRED
  #   if(identical(x[[2]], numeric(0))) {
  #     chr <- c("No suggestive or significant results")
  #     pos <- c("")
  #     names(chr) <- c(" ")
  #     x <- cbind(chr, pos)
  #   }
  
    # Puts the phenotype name in the structure with the significant markers (for all)
    # y <- cbind(phenotypes.name[ind], x[])
  
    # Output the qtl rds under the normal assumption
    saveRDS(out.i, path.out.rds)
  }
  else{
    message(paste0("Already Processed: ", filename, ": ", ind, "'th pheno; its name: ", phenotype.name))
  }
}

calculate.all.qtlmapping <- function(filename.set, tmp.in.path, tmp.out.path.rds, tmp.out.path.summary){
  # define path index
  tmp.file.from <- 1
  tmp.file.to <- length(filename.set)

  num.cores <- detectCores() - 2
  cl <- makeCluster(num.cores, outfile="")
  registerDoParallel(cl)
  
  # foreach(tmp.index = tmp.file.from: tmp.file.to, .export = c("calculate.single.qtlmapping", "scanone", "load.dataset", "pullout.phenotypes.name", "create.dir")) %do%{
  foreach(tmp.index = tmp.file.from: tmp.file.to) %do%{
      
    # get path and get file
    filename <- filename.set[tmp.index]
    message(paste0("Processing file ", filename, "......"))
    filename.withextension.rds <- paste0(filename, ".rds")
    filename.withextension.csv <- paste0(filename, ".csv")
    in.data <- load.dataset(tmp.in.path, filename.withextension.rds)

    # get the name of phenotypes, and define the start/end index.
    phenotypes.name <- pullout.phenotypes.name(in.data)
    index.phenotypes.name.from <- 1
    index.phenotypes.name.to <- length(phenotypes.name)

    # define path
    path.out.rds <- file.path(tmp.out.path.rds, filename)
    create.dir(path.out.rds, force=FALSE)
    path.out.summary <- file.path(tmp.out.path.summary, filename.withextension.rds)
    
    # output
    output.df <- data.frame()
    foreach(ind = index.phenotypes.name.from: index.phenotypes.name.to, .export = c("calculate.single.qtlmapping", "scanone", "load.dataset", "pullout.phenotypes.name", "create.dir")) %dopar% {
      calculate.single.qtlmapping(filename, phenotypes.name, ind, in.data, path.out.rds)
# 
#       names(out.i)[3] <- phenotypes.name[ind]
#       if(length(output.df) == 0) {
#         output.df <- as.data.frame(out.i)
#         output.df.rownames <- data.frame(rownames(output.df))
#         rownames(output.df) <- NULL
#         output.df <- cbind(output.df.rownames, output.df)
#       }
#       else {
#         tmp <- as.data.frame(out.i[3])
#         output.df <- cbind(output.df, tmp)
#       }
    }
    message(paste("Output the total result to the path: ", path.out.summary))
    saveRDS(output.df, path.out.summary)
  }
  stopCluster(cl)
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
path.data.in.geno_pheno <- paste(path.data.out, "geno_pheno_origin_rds", sep="/")
path.data.out.global <- paste(path.data.out, "global", sep="/")
path.data.out.phewas <- paste(path.data.out, "phewas", sep="/")
path.data.out.qtlmapping.rds <- paste(path.data.out, "qtlmapping", sep="/")
path.data.out.qtlmapping.summary <- paste(path.data.out, "qtlmapping_summary", sep="/")

# ------Load functions/parameters from other files
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------Load functions from other files
source(paste(path.code.data_processing, "utils.R", sep="/"))

# ------Parse the filename set from the filename list
filename.set <- unique(extract.filename.set(path.data.in.geno_pheno)$basename.with.type)

# ------Debug
# in.data.t <- load.dataset("/Users/hali/Dropbox/sv/gene_expression/new/data/output/geno_pheno_origin", "Kidney_male_BXD.rds")

# ------Run on Server
calculate.all.qtlmapping(filename.set, path.data.in.geno_pheno, path.data.out.qtlmapping.rds, path.data.out.qtlmapping.summary)
