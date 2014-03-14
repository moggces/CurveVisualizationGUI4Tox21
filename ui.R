
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

tabPanelAbout <- source("./source/about.R")$value

############################################

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Tox21 Concentration-Response Data Visualization"),
  
  sidebarPanel(
    h4('Mode'),
    radioButtons("mode", "Select a pathway display mode:",
                 choices = list("parallel"="parallel", "overlay"="overlay")),
    tags$br(),
    
    h4('Compound loader'),
    tags$textarea(id="cmpds", rows=3, cols=1, ""),
    helpText("please input either CAS, NCGC ID, or Tox21 ID"),

    
    tags$hr(),
    
    sliderInput("widthpx", 
                "width pixel/colmum", min = 150, max = 600, value = 300, step=50),
    sliderInput("heightpx", 
                "height pixel/colmum", min = 150, max = 600, value = 300, step=50),
    tags$br(),
    
    
    h4('Pathway'),
    wellPanel (
      uiOutput("pathways")
    ),
    
    h4('Pathway readout options'),
    wellPanel (
      uiOutput("options")
      #checkboxInput("isOneAssay", "multiplex cytotoxicity as one assay", TRUE)
    ),
    
    h4('Curve plotting options'),
    wellPanel (
      uiOutput("plot_options"),
      checkboxInput("showOutlier", "cross outliers", TRUE)
    ),

    br(),
    downloadButton('downloadPlot', 'Save Plot')
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel( 'Input compounds', dataTableOutput('contents')),
      tabPanel( "Plot", plotOutput("plot", height="auto", width="500%")),
      tabPanel("Data", dataTableOutput('qhts_data')),
      tabPanelAbout()
    )
  )
))
