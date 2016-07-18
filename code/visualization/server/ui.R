###############################################
# 2016-07-18
# Author: Tao Lin
###############################################
# Load package
library(shiny)
library(shinyjs)
library(qtlcharts)
library(shinysky)

# Design UI layout
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Analysis and Visualization in BXD"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    useShinyjs(),
    div(
      # name of gene
      htmlOutput("selectTextInput"),
      # source of data
      htmlOutput("selectDataSource"),
      # type of diet
      htmlOutput("selectDiet"),
      # type of Analysis
      htmlOutput("selectAnalysis"),
      # type of interactive
      htmlOutput("selectInteractive")
    ),
    # interactive
    br(),
    downloadButton(outputId = "downloadData", label = "Download CSV"),
    downloadButton(outputId = "downloadPlot", label = "Download plot")
  ),
  
  # Show a plot.
  mainPanel(
    conditionalPanel(condition="input.interactive == 0", plotOutput("lodPlot", height=530)),
    conditionalPanel(condition="input.interactive == 1", iplot_output("interactivePlot", height=530))
  )

))
