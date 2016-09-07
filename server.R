
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(grid)
library(Cairo)

#library(googleVis)

source( "./source/get.R",  local=TRUE)
source("./source/load.R",  local=TRUE)
source("./source/theme_complete_bw.R", local=TRUE)

#mapping_file <- './data/tox21_id_12807_cas_structure_4v4_080813_v4.txt'
#mapping_df <- read.delim(mapping_file, stringsAsFactors=FALSE)

# load chemical information
profile_file <- './data/tox21_mapping_v5a5.txt' #colunm name has to be GSID
mapping_df <- load_profile(profile_file) # global, dataframe output

# load assay related parameters
logit_para_file <- './data/tox21_assay_collection.txt'
assay_names <- load_profile(logit_para_file) # global, dataframe output



################

shinyServer(function(input, output) {

  # load chemicals
  data_chemical <- reactive({
    result <- NULL
    ids <- input$cmpds
    result <- unlist(strsplit(ids, "\n", perl=TRUE))
    return(result)
  })

  # load main pathway data
  pathway_data <- reactive({

    chemicals <- data_chemical()
    pathways <- input$paths

    result <- get_qhts_data_wrap(chemicals, pathways=pathways)

    return(result)
  })

  # load ch2/ch1 .. pathway data
  pathway_add_data <- reactive({

    chemicals <- data_chemical()
    result <- pathway_data()

    options <- input$opts
    options <- grep("^(?!run[1-3]$)", options,  perl=TRUE, value = TRUE)
    pathways <- input$paths

    result <- rbind.fill(result, get_qhts_data_wrap(chemicals, pathways=pathways, options=options))

    return(result)

  })

  # filter by unwanted channels
  data_filter <- reactive ({

    result <- pathway_add_data()
    options <- input$opts

    pattern <- NULL
    options1 <- grep("run[1-3]", options,  perl=TRUE, value = TRUE)
    options2 <- grep("blue|green|red", options,  perl=TRUE, value = TRUE)
    options1 <- sub("run", '.', options1)
    if (sum(options1 != '') > 0 ) pattern <- c(pattern, paste("\\", options1, sep="", collapse="|"))
    if (sum(options2 != '') > 0 ) pattern <- c(pattern, paste(options2, sep="", collapse="|"))

    pattern <- paste(pattern, sep="", collapse="|")

    if (is.null(pattern) | pattern == '')
    {
      result <- data.frame()
    } else
    {
      result <- result[grep(pattern, result$readout, perl=TRUE),]
    }

    # specifically for the situation of autoflu assays that there are no Cmpd_Library
    result <- get_relevant_cmpd_library(result)

    return(result)
  })

  # put the data into ggplot styple
  data_melter <- reactive ({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')

    qhts <- data_filter()

    result <- get_melt_data(qhts, resp_type=unique(c('raw', plot_options)))

    # get the purity rating here
    result <- join(result, subset(mapping_df, select=c(Tox21.ID, Purity_Rating_T0, Purity_Rating_T4)), by="Tox21.ID")

    if (nrow(result) > 0)
    {
      result[, "display_name"] <- paste(result$Chemical.Name,"\n",  result$CAS, "\n",
                                        "Purity:", paste(result$Purity_Rating_T0, result$Purity_Rating_T4, sep="|"), "\n",
                                        result$Tox21AgencyID, sep="")
      result <- result[order(result$display_name),]
    }

    return(result)
  })


  getVarHeight <- reactive({
    qhts <- data_filter()
    #qhts[, "display_name"] <- paste(qhts$CAS, "|\n", qhts$Tox21AgencyID, sep="")
    #qhts[, "display_name"] <- paste(qhts$Chemical.Name,"|\n",  qhts$CAS, "|\n", qhts$Tox21AgencyID, sep="")
    #nrow <- length(unique(qhts$display_name))
    nrow <- length(unique(qhts$Chemical.ID))
    mode <- input$mode
    heightpx <- input$heightpx
    if (mode == 'overlay')
    {
      return(nrow * heightpx) # 350
    } else if (mode == 'parallel'| mode == 'mixed')
    {
      return(nrow * heightpx) # 300
    }
  })

  getVarWidth <- reactive({
    qhts <- data_filter()
    mode <- input$mode
    ncol <- length(unique(qhts$pathway))
    widthpx <- input$widthpx
    if (mode == 'overlay')
    {
      #return("auto")
      return(1000)
    } else if (mode == 'parallel'| mode == 'mixed')
    {
      return(ncol * widthpx)
    }
  })

  # illustrate the pathway names to display in interface
  output$pathways <- renderUI({
    selectizeInput("paths", "Select one or more assays to show:",
                choices  = list(
                                stress_response = c(
                                'activation_AP1' = 'tox21-ap1-agonist-p1_ratio',
                                "activation_ATAD5"="tox21-elg1-luc-agonist-p1_luc",
                                'activation_EndoRS' = 'tox21-esre-bla-p1_ratio',
                                "activation_H2AX"="tox21-h2ax-cho-p2_ratio",
                                "activation_HIF1A"="tox21-hre-bla-agonist-p1_ratio",
                                'activation_HSP'='tox21-hse-bla-p1_ratio',
                                'inhibition_MMP'='tox21-mitotox-p1_ratio',
                                'activation_NFkb' = 'tox21-nfkb-bla-agonist-p1_ratio',
                                'activation_Nrf2' = 'tox21-are-bla-p1_ratio',
                                "activation_p53"="tox21-p53-bla-p1_ratio",
                                "activation_DNA_damage/srf"="tox21-dt40-srf-p1_luc", "activation_DNA_damage/dsb"="tox21-dt40-dsb-p1_luc",
                                "inhibition_aromatase"="tox21-aromatase-p1_luc"),
                                nuclear_receptor = c(
                                'agonism_AhR' = 'tox21-ahr-p1_luc',
                                'agonism_AR/partial'='tox21-ar-bla-agonist-p1_ratio',
                                'agonism_AR/full'='tox21-ar-mda-kb2-luc-agonist-p1_luc',
                                'antagonism_AR/partial'='tox21-ar-bla-antagonist-p1_ratio',
                                'antagonism_AR/full/run1'='tox21-ar-mda-kb2-luc-antagonist-p1_luc',
                                'antagonism_AR/full/run2'='tox21-ar-mda-kb2-luc-antagonist-p2_luc',
                                'agonism_CAR'='tox21-car-agonist-p1_luc',
                                'antagonism_CAR'='tox21-car-antagonist-p1_luc',
                                'agonism_ER/partial'='tox21-er-bla-agonist-p2_ratio',
                                'agonism_ER/full'='tox21-er-luc-bg1-4e2-agonist-p2_luc',
                                'antagonism_ER/partial'='tox21-er-bla-antagonist-p1_ratio',
                                'antagonism_ER/full/run1'='tox21-er-luc-bg1-4e2-antagonist-p1_luc',
                                'antagonism_ER/full/run2'='tox21-er-luc-bg1-4e2-antagonist-p2_luc',
                                'agonism_GR' = 'tox21-gr-hela-bla-agonist-p1_ratio',
                                'antagonism_GR' = 'tox21-gr-hela-bla-antagonist-p1_ratio',
                                'agonism_FXR'='tox21-fxr-bla-agonist-p2_ratio',
                                "antagonism_FXR"="tox21-fxr-bla-antagonist-p1_ratio",
                                'agonism_PPARg' = 'tox21-pparg-bla-agonist-p1_ratio',
                                'antagonism_PPARg' = 'tox21-pparg-bla-antagonist-p1_ratio',
                                'agonism_PPARd' = 'tox21-ppard-bla-agonist-p1_ratio',
                                'antagonism_PPARd' = 'tox21-ppard-bla-antagonist-p1_ratio',
                                'antagonism_RORr' = 'tox21-ror-cho-antagonist-p1_luc',
                                'agonism_RSP' = 'tox21-rar-agonist-p1_luc',
                                'agonism_RXR' = 'tox21-rxr-bla-agonist-p1_ratio',
                                'agonism_TR' = 'tox21-gh3-tre-agonist-p1_luc',
                                "antagonism_TR"="tox21-gh3-tre-antagonist-p1_luc",
                                'agonism_VDR' = 'tox21-vdr-bla-agonist-p1_ratio',
                                'antagonism_VDR' = 'tox21-vdr-bla-antagonist-p1_ratio'),
                                counter_screen = c(
                                'luciferase_toxicity' = 'tox21-luc-biochem-p1_luc',
                                'autofluor_hek293/cell'='tox21-spec-hek293-p1_cell-_main',
                                'autofluor_hek293/medium'='tox21-spec-hek293-p1_medi-_main',
                                'autofluor_hepg2/cell'='tox21-spec-hepg2-p1_cell-_main',
                                'autofluor_hepg2/medium'='tox21-spec-hepg2-p1_medi-_main')
                                #'autofluor_hek293/cell'='spec-hek293_cell_main'
                                #'autofluor_hek293/medium'='spec-hek293_medi_main',
                                #'autofluor_hepg2/cell'='spec-hepg2_cell_main',
                                #'autofluor_hepg2/medium'='spec-hepg2_medi_main'
                                ),
                multiple = TRUE)
  })

  output$options <- renderUI({
    selectInput("opts", "Select readout options:",
                choices  = list( "run1"="run1","run2"="run2", "run3"="run3",
                                 "cytotoxicity"="via",
                                 "BLA-signal"="ch2",
                                 "BLA-background"="ch1",
                                 "AutoF-blue"="blue-n",
                                 "AutoF-green"="green",
                                 "AutoF-red"="red"
                                 ),
                selected=c("run1", "run2", "run3"),
                multiple = TRUE)
  })

  output$plot_options <- renderUI({
    selectInput("plt_opts", "Select line plotting options:",
                choices  = list( "raw"="raw", "curvep"="curvep",
                                 "Hill 4-point"="hill"),
                selected=c("raw"),
                multiple = TRUE)
  })


  output$contents <- renderDataTable({

    if ( length(data_chemical()) > 0 )
    {
      get_mapping_data(data_chemical(), mapping_df)
    }
  })

  output$assay_info <- renderDataTable({

    col_n <- c('assay','common_name','technology','cell_type','species','abbreviation', 'PubChem AID')
    result <- assay_names[, colnames(assay_names) %in% col_n]
    return(result)

  })

  output$qhts_data <- renderDataTable({
    data_filter()
  })

  output$plot <- renderPlot({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')

    report_format <- input$report_format

    # generate plotting parameters (for future)
    rm_raw_color <- FALSE #input$rmRawColor
    rm_raw_line <- FALSE #input$rmRawLine
    hl_pod <- FALSE #input$hlpod
    hd_error_bar <- FALSE #input$hdErrorBar
    xaxis_range <- NULL
    yaxis_range <- NULL
    if(input$xaxisLogical) xaxis_range <- input$xaxis
    if(input$yaxisLogical) yaxis_range <- input$yaxis

    paras <- list(report_format=report_format, rm_raw_color=rm_raw_color,
                  rm_raw_line=rm_raw_line,hl_pod=hl_pod, hd_error_bar=hd_error_bar,
                  xaxis_range=xaxis_range, yaxis_range=yaxis_range)


    # melt the data
    result <- data_melter()
    p <- get_plot(result, mode=mode, plot_options=plot_options, fontsize=20, pointsize=3, paras=paras)

    if (mode == 'overlay')
    {
      p <- p  + theme_bw(base_size = 20) + facet_wrap(~ display_name  , ncol=2)
    } else if (mode == 'parallel')
    {
      p <- p + theme_bw(base_size = 20) + facet_grid(display_name ~ pathway)
    } else if (mode == 'mixed')
    {
      p <- p  + theme_bw(base_size = 20) + facet_wrap(~ pathway  , ncol=2)
    }
    if (report_format) p <- p + theme_complete_bw(base_size = 20) + guides(colour=FALSE)
    print(p)
    #print(select_plot())
  }, height=getVarHeight, width=getVarWidth)

  output$downloadRData <- downloadHandler(
    filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".Rdata", sep="") },
    content = function(file) {
      result <- data_melter()
      result <- get_published_data_only(result, assay_names)
      save(result, file=file)
    })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".pdf", sep="") },
    content = function(file) {
      result <- data_melter()
      #result <- get_published_data_only(result, assay_names)
      mode <- input$mode
      plot_options <- input$plt_opts
      show_outlier <- input$showOutlier
      if (show_outlier) plot_options <- c(plot_options, 'mask')
      report_format <- input$report_format

      # generate plotting parameters (for future)
      rm_raw_color <- FALSE #input$rmRawColor
      rm_raw_line <- FALSE #input$rmRawLine
      hl_pod <- FALSE #input$hlpod
      hd_error_bar <- FALSE #input$hdErrorBar
      xaxis_range <- NULL
      yaxis_range <- NULL
      if(input$xaxisLogical) xaxis_range <- input$xaxis
      if(input$yaxisLogical) yaxis_range <- input$yaxis

      paras <- list(report_format=report_format, rm_raw_color=rm_raw_color,
                    rm_raw_line=rm_raw_line,hl_pod=hl_pod, hd_error_bar=hd_error_bar,
                    xaxis_range=xaxis_range, yaxis_range=yaxis_range)

      # pdf file is the output file
      pdf(file=file)

      if (mode == 'overlay')
      {
        n_page <- 6
        result <- get_blank_data(result, n_page)
        nn <- unique(result$display_name)
        pages <- split(nn, ceiling(seq_along(nn)/n_page))
        lapply(names(pages), function (x) {
          sub <- result[result$display_name %in% pages[[x]],]
          p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1, paras=paras)
          p <- p  + theme_bw(base_size = 8) + facet_wrap(~ display_name  , ncol=2, nrow=3)
          if (report_format) p <- p + theme_complete_bw(base_size = 8) + guides(colour=FALSE)
          print(p)
        })

      } else if (mode == 'parallel')
      {
        n_page <- 6
        result <- get_blank_data(result, n_page)
        nn <- unique(result$display_name)
        pages <- split(nn, ceiling(seq_along(nn)/n_page))
        lapply(names(pages), function (x) {
          sub <- result[result$display_name %in% pages[[x]],]
          p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1, paras=paras)
          p <- p + theme_bw(base_size = 8) + facet_grid(display_name ~ pathway)
          if (report_format) p <- p + theme_complete_bw(base_size = 8) + guides(colour=FALSE)
          print(p)
        })

      } else if (mode == 'mixed')
      {

        n_page <- 6
        result <- get_blank_data(result, n_page, 'pathway')
        nn <- unique(result$pathway)
        pages <- split(nn, ceiling(seq_along(nn)/n_page))
        lapply(names(pages), function (x) {
          sub <- result[result$pathway %in% pages[[x]],]
          p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1, paras=paras)
          p <- p  + theme_bw(base_size = 8) + facet_wrap(~ pathway  , ncol=2)
          if (report_format) p <- p + theme_complete_bw(base_size = 8) + guides(colour=FALSE)
          print(p)
        })

      }

      dev.off()
    })
})

