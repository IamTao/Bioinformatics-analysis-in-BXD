###############################################
# 2016-07-18
# Author: Tao Lin
###############################################
# Load package.
library(shiny)
library(qtl)
library(qtlcharts)
library(intermediate)
library(biomaRt)
library(plyr)
library(ggplot2)
library(data.table)
library(shinyjs)
library(shinysky)
library(stringr)

# define path
path.in.data <- "/Users/hali/Documents/gene_expression/new/data/output"
# path.in.data <- "/datadisk/gene_expression/data/output"
path.in.global <- file.path(path.in.data, "global")

# Load dataset
gene_aligner <- readRDS(file.path(path.in.global, "gene_position_aligner.rds"))
gene_aligner <- data.table(gene_aligner)

# Server.
shinyServer(function(input, output, session){
  # Gene symbol input, default from URL address if given
  output$selectTextInput <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$gene
    if (is.null(default)) default <- "Ace"
    textInput("gene_id", "Gene Symbol", value=default)
    # select2Input("gene_id","Gene Symbol", choices=c(my_autocomplete_list), type = c("input", "select"))
  })
  # select the data source.
  output$selectDataSource <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$data_source
    if (is.null(default)) default <- "Midbrain"
    selectInput("data_source", "Data Source", c("BAT" = "BAT", "Eye" = "Eye", "Heart" = "Heart", "Hippocampus" = "Hippocampus", "Kidney" = "Kidney", "Liver" = "Liver",
                                                "Liver protein" = "LiverProt", "Metabolite" = "Metabolite", "Midbrain" = "Midbrain", "Muscle" = "Muscle", "NAc" = "NAc",
                                                "PFC" = "PFC", "scWAT" = "scWAT", "Spleen" = "Spleen", "Striatum" = "Striatum"), default)
  })
  # select the type of diet.
  output$selectDiet <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$diet
    if (is.null(default)) default <- "CD"
    radioButtons("diet", "Diet", c("CD" = "CD", "HFD" = "HFD"), inline=TRUE, selected = default)
  })
  # select the type of Analysis
  output$selectAnalysis <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$analysis
    if (is.null(default)) default <- "PheWAS"
    radioButtons("analysis", "Analysis",
                 c("PheWAS", "QTLMapping", "Mediation"), selected = default)
  })
  # select the type of interactive.
  output$selectInteractive <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$diet
    if (is.null(default)) default <- 0
    radioButtons("interactive", "Interactive plots", c("off" = 0, "on" = 1), inline=TRUE, selected = default)
  })
  
  # reset the value of interactive plot.
  observeEvent(input$analysis == 'QTLMapping', {
    reset("interactive")
  })
  
  # plot
  output$lodPlot <- renderPlot({create.plot(input$gene_id, input$data_source, input$diet, input$analysis)})
  # iplot
  output$interactivePlot <- iplot_render({create.iplot(input$gene_id, input$data_source, input$diet, input$analysis)})

  # download
  output$downloadData <- downloadHandler(
    filename = function() {paste0(input$gene_id, '_', input$data_source, '_', input$diet, '_', input$analysis, '.Rds')},
    content = function(file) {
      plot_type     <- decide.plot_type(input$data_source, input$diet, input$analysis)
      file_location <- get.file_location(input$gene_id, plot_type)
    
      if (!file.exists(file_location))
        return(NULL)
      tmp <- readRDS(file_location)
      saveRDS(tmp, file)
    }
  )

  output$downloadPlot <- downloadHandler(
    # specify the file name
    filename = function() {paste0(input$gene_id, '_', input$data_source, '_', input$diet, '_', input$analysis, '.pdf')},
    content = function(file){
      # open the device
      pdf(file, width = 10, height = 7)
      # make the plot
      create.plot(input$gene_id, input$data_source, input$diet, input$analysis)
      # close the device
      dev.off()
    }
  )
  
})

# update auto-complete list
update.autolist <- function(){
  all_files  <- list.files(data_path,  recursive = TRUE)
  useless <- c("gene_position_aligner", "pheno_phewas", "phenotypes_id_aligner", "snp_location_aligner", "best_protein_qtl", "processed_protein_info")
  autolist <- laply(all_files, function(x){
    splitted <- strsplit(x, split = '/')
    tmp <- splitted[[1]][length(splitted[[1]])]
    final <- substr(tmp, start = 1, stop = nchar(tmp) - 4)
    if (final %in% useless) ""
    else final
  })
  return(unique(autolist))
}
my_autocomplete_list <- list() # update.autolist()

# No data plot
nodata <- function(info) {
  plot(0, type = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", main=info)
  text(1, 0, "No data", cex=5)
}

# Decide the plot type based on selected analysis approach.
decide.plot_type <-function(analysis, data_source, diet){
  if (tolower(analysis) == "phewas")
    return (file.path(tolower(analysis), paste(data_source, diet, sep="_")))
  else if (tolower(analysis) == "qtlmapping")
    return(file.path(tolower(analysis), paste(data_source, diet, sep="_")))
  else if (tolower(analysis) == "mediation")
    return(file.path(tolower(analysis), paste(data_source, diet, sep="_")))
}

# get the location of a file based on its gene name and plot type.
get.file_location <- function(gene_name, plot_type){
  return(paste(path.in.data, plot_type, paste0(gene_name, ".Rds"), sep="/"))
}

# Given the name of gene, get the information of gene.
get.gene_info <- function(gene_name){
  # The format of input: a character string
  # The format of output: either NULL/Dataframe(name, chr, start_pos)
  if (is.null(gene_name)){
    return(NULL)
    }
  else{
    gene_info    <- gene_aligner[GeneName == gene_name]
    return(list("name"=as.character(gene_info$GeneName), "chr"=as.numeric(gene_info$chromosome), "pos"=gene_info$start_position))
    }
}

# Create the plot.
create.plot <- function(gene_name, data_source, diet, analysis){
  if (is.null(gene_name)) return(NULL)
  
  plot_type     <- decide.plot_type(analysis, data_source, diet)
  file_location <- get.file_location(gene_name, plot_type)
  gene_info     <- get.gene_info(gene_name)
  
  if (!file.exists(file_location)) {
    nodata(paste0("Path: '", file_location, " 'not found"))
    return(NULL)
  } 
  data_table    <- readRDS(file_location)

  if (grepl("phewas", plot_type)){
    rotation          <- 40
    sigline           <- 3.5
    given_gene_y      <- data_table     
    given_gene_y$Category_ID  <- as.numeric(given_gene_y$Category_ID); given_gene_y$lod  <- as.numeric(given_gene_y$lod)
    y_limits          <- ceiling(max(given_gene_y$lod, na.rm = TRUE))
    x_limits          <- 4466
    cat               <- levels(given_gene_y$Category)
    pos_labels        <- vector(length=length(cat))
    for (i in 1 : length(cat)){
      a <- cat[i]
      b <- subset(given_gene_y, given_gene_y$Category == a)
      pos_labels[i] <- ((max(b$Category_ID, na.rm = TRUE) - min(b$Category_ID, na.rm = TRUE))/2) + min(b$Category_ID, na.rm = TRUE)
      rm(a,b)
    }     
    plotmain=paste("PheWAS Manhattan Plot for ", gene_info$name, " (on Chr", gene_info$chr, ": ", gene_info$start_pos, " Mb)", sep="")     
    # make normal plot
    par(mar=c(12,6,5,3))
    plot(given_gene_y$Category_ID[(given_gene_y$lod >= 0) & (given_gene_y$lod < 2)], 
         given_gene_y$lod[(given_gene_y$lod >= 0) & (given_gene_y$lod < 2)], 
         cex = 1, pch = 20, ylim = c(0, y_limits), xlim = c(0,x_limits), 
         col = c("#3953A4","#93C83D","#F68C1E","#6FCCDD")[given_gene_y$color[(given_gene_y$lod >= 0) & (given_gene_y$lod < 2)]], 
         ylab = "LOD",main = plotmain, cex.main = 2, cex.lab = 1.8, cex.axis = 1.8, xlab = "", xaxt = "n")
    par(new=TRUE)
    plot(given_gene_y$Category_ID[(given_gene_y$lod >= 2) & (given_gene_y$lod < 3)], 
         given_gene_y$lod[(given_gene_y$lod >= 2) & (given_gene_y$lod < 3)], 
         cex = 1.5, pch = 20, ylim = c(0, y_limits), xlim = c(0, x_limits),
         col = c("#3953A4","#93C83D","#F68C1E","#6FCCDD")[given_gene_y$color[(given_gene_y$lod>=2) & (given_gene_y$lod < 3)]], 
         xlab = "", ylab ="", xaxt = "n", yaxt = "n")
    par(new=TRUE)
    plot(given_gene_y$Category_ID[(given_gene_y$lod >= 3) & (given_gene_y$lod < 4)], 
         given_gene_y$lod[(given_gene_y$lod >= 3) & (given_gene_y$lod < 4)], 
         cex = 2, pch = 20, ylim = c(0,y_limits), xlim = c(0, x_limits),
         col=c("#3953A4","#93C83D","#F68C1E","#6FCCDD")[given_gene_y$color[(given_gene_y$lod >= 3) & (given_gene_y$lod < 4)]], 
         xlab = "", ylab = "", xaxt = "n",yaxt = "n")
    par(new=TRUE)
    plot(given_gene_y$Category_ID[given_gene_y$lod >= 4], 
         given_gene_y$lod[given_gene_y$lod >= 4], 
         cex = 3, pch = 20, ylim = c(0,y_limits), xlim = c(0,x_limits),
         col = c("#3953A4","#93C83D","#F68C1E","#6FCCDD")[given_gene_y$color[given_gene_y$lod >= 4]], 
         xlab = "", ylab = "", xaxt = "n", yaxt = "n")      
    text(pos_labels, par("usr")[3] - 0.25, srt = rotation, adj = 1,labels = levels(given_gene_y$Category), xpd = TRUE, cex = 1.5)
    abline(h = sigline, col = "red", lty = 2) # we need to calculate the threshold for significance     
    for (i in 1:max(given_gene_y$Category_ID, na.rm = TRUE)){  # add the labels for significant phenotypes
      if (given_gene_y$lod[i] > sigline){
        labels<-toString(given_gene_y$Shown_pheno[i])
        text(given_gene_y$Category_ID[i], given_gene_y$lod[i], labels=toString(given_gene_y$Shown_pheno[i]), cex=1.2, pos=4)
      }
    }

  }
  else if (grepl("qtlmapping", plot_type)){
    plot(data_table, main=paste0("QTL Plot for ", gene_info$name, ", ", plot_type))
  }
  else if (grepl("mediation", plot_type)){
    med <- data_table[data_table$chr %in% c(as.character(1:19),"X"),] # to avoid problems with Y, M and NA
    plot(med, main=paste0("Mediation Plot for ", "gene_info$name", ", ", "plot_type"))
  }
  else{
    nodata("Plot type not specified")
    return(NULL)
  }

 }

# Create the iplot (interactive plot)
create.iplot <- function(gene_name, data_source, diet, analysis){
  if (is.null(gene_name)) return(NULL)
  
  plot_type     <- decide.plot_type(data_source, diet, analysis)
  file_location <- get.file_location(gene_name, plot_type)
  gene_info     <- get.gene_info(gene_name)

  if (!file.exists(file_location)) {
    return(NULL)
  } 

  data_table    <- readRDS(file_location)

  if (grepl("phewas", plot_type)){
    plotmain=paste("PheWAS Manhattan Plot for ", gene_info$name, " (on Chr", gene_info$chr, ": ", gene_info$start_pos, " Mb)", sep="")     

    # make interactive plot by iplot
    given_gene_y              <- data_table     
    given_gene_y$Category_ID  <- as.numeric(given_gene_y$Category_ID); given_gene_y$lod  <- as.numeric(given_gene_y$lod)
    p0 <- iplot(x = given_gene_y$Category_ID, y = given_gene_y$lod, group = given_gene_y$color, indID = given_gene_y$Shown_pheno,
                chartOpts = list(xlab="Category", xlab = "", xaxt = "n", ylab = "LOD", title = plotmain, height = 700, width = 1000, pointcolor = c("#3953A4","#93C83D","#F68C1E","#6FCCDD")))
    return(p0)
  }
  else if (grepl("mediation", plot_type)){
    med <- data_table[data_table$chr %in% c(as.character(1:19),"X"),]
    kplot(med)
  }
  else{
    return(NULL)    
  }
}

