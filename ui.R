library(shiny)

shinyUI(
  bootstrapPage(
    tags$head(
      tags$meta(`http-equiv`="pragma", content="no-cache"),
      tags$meta(`http-equiv`="Cache-control", content="no-cache, no-store")
    ),
    tags$head(tags$script(src="extra.js")),
    tags$head(tags$link(rel="stylesheet", type="text/css", href="headerfooter.css")),
    tags$head(tags$link(rel="stylesheet", type="text/css", href="custom.css")),
    includeHTML(path="./www/header.html"),
    tabsetPanel(
      tabPanel('About', tags$div(includeMarkdown("README.Rmd"), class='container')),
      tabPanel('Application',
        pageWithSidebar(
          headerPanel('', windowTitle="Tox21 Curve Browser"),
          sidebarPanel(

            h3('Required inputs'),
            wellPanel(
              h4('Compound loader'),
              tags$textarea(id="cmpds", rows=5, class='col-xs-12', ""),
              helpText("Input either CAS, NCGC ID, or Tox21 ID"),
              tags$hr(),

              h4('Assays'),
              uiOutput("pathways")
            ),

            h3('Other settings'),
            wellPanel (
              radioButtons("mode", "Select an assay display mode:",
                choices = list(
                  "assay parallel"="parallel",
                  "assay overlaid"="overlay",
                  "assay parallel + cmpd overlaid"="mixed")),
              tags$hr(),

              uiOutput("options"),
              tags$hr(),

              uiOutput("plot_options"),
              checkboxInput("showOutlier", "Cross outliers", FALSE),
              tags$hr(),

              checkboxInput("report_format", "Use template for reporting", FALSE)
            ),

            h3('Formatting & outputs'),
            wellPanel (
              h4('Resize plot image'),
              sliderInput("widthpx",
                          "width pixel/colmum", min = 150, max = 1000, value = 300, step=50),
              sliderInput("heightpx",
                          "height pixel/colmum", min = 150, max = 1000, value = 300, step=50),

              tags$hr(),

              h4('Resize axes'),
              checkboxInput("xaxisLogical", "Manually adjust x-axis", FALSE),

              conditionalPanel(
                condition = "input.xaxisLogical !== false",
                sliderInput("xaxis",
                            "x-axis", min = -9, max = -3, value = c(-8,-4), step=0.5)
              ),

              checkboxInput("yaxisLogical", "Manually adjust y-axis", FALSE),
              conditionalPanel(
                condition = "input.yaxisLogical !== false",
                sliderInput("yaxis",
                            "y-axis", min = -500, max = 500, value = c(0,100), step=50)
              ),

              tags$hr(),
              downloadButton('downloadRData', 'Save Melted Rdata'),
              tags$span(' '),
              downloadButton('downloadPlot', 'Save Plot')
            )
          ),
          mainPanel(
            tabsetPanel(
              tabPanel('Input compounds', dataTableOutput('contents')),
              tabPanel("Plot", plotOutput("plot", height="auto", width="500%")),
              tabPanel("Data", dataTableOutput('qhts_data')),
              tabPanel('Assays', dataTableOutput('assay_info'))
            )
          )
        )
      )
    ),
    includeHTML(path="./www/footer.html")
  )
)

# aasssssssssssss
