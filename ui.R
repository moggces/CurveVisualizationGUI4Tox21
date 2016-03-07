
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

#tabPanelAbout <- source("./source/about.R")$value


############################################

shinyUI(
  
  
  pageWithSidebar(
  
  # Application title
  headerPanel("Tox21 Curve Browser"),
  
  sidebarPanel(
    
    tags$head(tags$meta(`http-equiv`="pragma", content="no-cache"), 
              tags$meta(`http-equiv`="Cache-control", content="no-cache, no-store")),
  
    h4('Assays'),
    wellPanel (
      uiOutput("pathways")
    ),
    

    
    h4('Mode'),
    radioButtons("mode", "Select a assay display mode:",
                 choices = list("assay parallel"="parallel", "assay overlaid"="overlay", 
                                "assay parallel + cmpd overlaid"="mixed")),
    tags$hr(),
    
    h4('Compound loader'),
    tags$textarea(id="cmpds", rows=3, cols=1, ""),
    helpText("please input either CAS, NCGC ID, or Tox21 ID"),
    
    tags$hr(),

    h4('Assay readout options'),
    wellPanel (
      uiOutput("options")
      #checkboxInput("isOneAssay", "multiplex cytotoxicity as one assay", TRUE)
    ),
    
    
    h4('Curve fitting options'),
    wellPanel (
      uiOutput("plot_options"),
      checkboxInput("showOutlier", "cross outliers", FALSE)
    ),
    
    tags$hr(),
    
    h4('Others'),
    checkboxInput("report_format", "use template for reporting", FALSE), 
    
    tags$hr(),
    
    h4('Resize plot image'),
    sliderInput("widthpx", 
                "width pixel/colmum", min = 150, max = 1000, value = 300, step=50),
    sliderInput("heightpx", 
                "height pixel/colmum", min = 150, max = 1000, value = 300, step=50),
    
    tags$hr(),
    
    h4('Resize axes'),
    
    checkboxInput("xaxisLogical", "manually adjust x axis", FALSE),
    
    
    conditionalPanel(
      condition = "input.xaxisLogical !== false",
      sliderInput("xaxis", 
                  "x-axis", min = -9, max = -3, value = c(-8,-4), step=0.5)
    ),
    
    checkboxInput("yaxisLogical", "manually adjust y axis", FALSE),
    
    conditionalPanel(
      condition = "input.yaxisLogical !== false",
      sliderInput("yaxis", 
                  "y-axis", min = -500, max = 500, value = c(0,100), step=50)
    ),
    
    
    tags$hr(),
    
    downloadButton('downloadRData', 'Save Melted Rdata'),
    
    tags$br(),
    
    downloadButton('downloadPlot', 'Save Plot')
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel( 'Input compounds', dataTableOutput('contents')),
      tabPanel( "Plot", plotOutput("plot", height="auto", width="500%")),
      tabPanel("Data", dataTableOutput('qhts_data')),
      tabPanel( 'Assays', dataTableOutput('assay_info')),
      tabPanel('About', includeHTML("README.html"))
    )
  )
))
