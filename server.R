
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
#library(googleVis)

source( "./source/get.R",  local=TRUE)
source("./source/load.R",  local=TRUE)

#mapping_file <- './data/tox21_id_12807_cas_structure_4v4_080813_v4.txt'
#mapping_df <- read.delim(mapping_file, stringsAsFactors=FALSE)

# load chemical information (will include purity later)
profile_file <- './data/tox21_compound_id_v5a.txt' #colunm name has to be GSID
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
    options2 <- grep("cell|medi", options,  perl=TRUE, value = TRUE)
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
    if (nrow(result) > 0)
    {
      result[, "display_name"] <- paste(result$Chemical.Name,"|\n",  result$CAS, "|\n", result$Tox21AgencyID, sep="")
      #result[, "display_name"] <- paste(result$CAS, "|\n", result$Tox21AgencyID, sep="")
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
    selectInput("paths", "Select pathways to show:", 
                choices  = list("activation_ATAD5"="elg1-luc-agonist_luc", "activation_p53"="p53-bla_ratio",
                                "activation_DNA_damage/srf"="dt40-srf_luc", "activation_DNA_damage/dsb"="dt40-dsb_luc",
                                'agonism_AhR' = 'ahr_luc',
                                'agonism_AR/partial'='ar-bla-agonist_ratio',
                                'agonism_AR/full'='ar-mda-kb2-luc-agonist_luc',
                                'antagonism_AR/partial'='ar-bla-antagonist_ratio',
                                'antagonism_AR/full'='ar-mda-kb2-luc-antagonist_luc',
                                'agonism_ER/partial'='er-bla-agonist_ratio',
                                'agonism_ER/full'='er-luc-bg1-4e2-agonist_luc',
                                'antagonism_ER/partial'='er-bla-antagonist_ratio',
                                'antagonism_ER/full'='er-bla-antagonist_ratio',
                                'agonism_GR' = 'gr-hela-bla-agonist_ratio',
                                'antagonism_GR' = 'gr-hela-bla-antagonist_ratio',
                                'agonism_FXR'='fxr-bla-agonist_ratio',
                                "antagonism_FXR"="fxr-bla-antagonist_ratio",
                                'agonism_PPARg' = 'pparg-bla-agonist_ratio',
                                'antagonism_PPARg' = 'pparg-bla-antagonist_ratio',
                                'agonism_PPARd' = 'ppard-bla-agonist_ratio',
                                'antagonism_PPARd' = 'ppard-bla-antagonist_ratio',
                                'agonism_TR' = 'gh3-tre-agonist_luc',
                                "antagonism_TR"="gh3-tre-antagonist_luc",
                                "inhibition_aromatase"="aromatase_luc",
                                'inhibition_MMP'='mitotox_ratio',
                                'activation_Nrf2' = 'are-bla_ratio',
                                'activation_HSP'='hse-bla_ratio',
                                'activation_EndoRS' = 'esre-bla_ratio',
                                'activation_NFkb' = 'nfkb-bla-agonist_ratio',
                                'agonism_RSP' = 'rar-agonist_luc',
                                'agonism_RXR' = 'rxr-bla-agonist_ratio',
                                'agonism_VDR' = 'vdr-bla-agonist_ratio',
                                'antagonism_VDR' = 'vdr-bla-antagonist_ratio',
                                'antagonism_RORr' = 'ror-cho-antagonist_luc',
                                'activation_AP1' = 'ap1-agonist_ratio',
                                'autofluor_hek293/cell'='spec-hek293_cell',
                                'autofluor_hek293/medium'='spec-hek293_medi',
                                'autofluor_hepg2/cell'='spec-hepg2_cell',
                                'autofluor_hepg2/medium'='spec-hepg2_medi'
                                ), 
                multiple = TRUE)
  })
  
  output$options <- renderUI({
    selectInput("opts", "Select readout options:", 
                choices  = list( "run1"="run1","run2"="run2", "run3"="run3",
                                 "cytotoxicity"="via", 
                                 "BLA-signal"="ch2", 
                                 "BLA-background"="ch1",
                                 "AutoF-blue"="blue",
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
    
    col_n <- c('common_name','technology','cell_type','species','abbreviation')
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
    
    # generate plotting parameters (for future)
    rm_raw_color <- FALSE #input$rmRawColor
    rm_raw_line <- FALSE #input$rmRawLine
    hl_pod <- FALSE #input$hlpod
    hd_error_bar <- FALSE #input$hdErrorBar
    paras <- list(rm_raw_color=rm_raw_color, rm_raw_line=rm_raw_line,hl_pod=hl_pod, hd_error_bar=hd_error_bar )
    
    
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
    print(p)
    #print(select_plot())
  }, height=getVarHeight, width=getVarWidth)
  
  output$downloadRData <- downloadHandler(
    filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".Rdata", sep="") },
    content = function(file) {
      result <- data_melter()
      save(result, file=file)
    })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".pdf", sep="") },
    content = function(file) {
      result <- data_melter()
      mode <- input$mode
      plot_options <- input$plt_opts
      show_outlier <- input$showOutlier
      if (show_outlier) plot_options <- c(plot_options, 'mask')
      
      # generate plotting parameters (for future)
      rm_raw_color <- FALSE #input$rmRawColor
      rm_raw_line <- FALSE #input$rmRawLine
      hl_pod <- FALSE #input$hlpod
      hd_error_bar <- FALSE #input$hdErrorBar
      paras <- list(rm_raw_color=rm_raw_color, rm_raw_line=rm_raw_line,hl_pod=hl_pod, hd_error_bar=hd_error_bar )
      
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
          print(p)
        })
        
      }
      
      dev.off()
    }
  )
  
#   output$plot <- renderPlot({
#     mode <- input$mode
#     plot_options <- input$plt_opts
#     show_outlier <- input$showOutlier
#     if (show_outlier) plot_options <- c(plot_options, 'mask')
#     
#     result <- data_melter()
#     p <- get_plot(result, mode=mode, plot_options=plot_options, fontsize=20, pointsize=3)
#     
#     if (mode == 'overlay')
#     {
#       p <- p  + facet_wrap(~ display_name  , ncol=2)
#     } else if (mode == 'parallel')
#     {
#       p <- p + facet_grid(display_name ~ pathway)
#     }
#     print(p)
#     #print(select_plot())
#   }, height=getVarHeight, width=getVarWidth)
# 
# 
# 
#   
#   output$downloadPlot <- downloadHandler(
#          filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".pdf", sep="") },
#          content = function(file) {
#            result <- data_melter()
#            mode <- input$mode
#            plot_options <- input$plt_opts
#            show_outlier <- input$showOutlier
#            if (show_outlier) plot_options <- c(plot_options, 'mask')
#       
#            pdf(file)
#            
#            if (mode == 'overlay')
#            {
#              n_page <- 6
#              result <- get_blank_data(result, n_page)
#              nn <- unique(result$display_name)
#              pages <- split(nn, ceiling(seq_along(nn)/n_page))
#              lapply(names(pages), function (x) {
#                sub <- result[result$display_name %in% pages[[x]],]
#                p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
#                p <- p  + facet_wrap(~ display_name  , ncol=2, nrow=3)
#                print(p)
#              })
#              
#            } else if (mode == 'parallel')
#            {
#              n_page <- 6
#              result <- get_blank_data(result, n_page)
#              nn <- unique(result$display_name)
#              pages <- split(nn, ceiling(seq_along(nn)/n_page))
#              lapply(names(pages), function (x) {
#                sub <- result[result$display_name %in% pages[[x]],]
#                p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
#                p <- p + facet_grid(display_name ~ pathway)
#                print(p)
#              })
#              
#            }
#            
#            dev.off()
#          }
#   )
  
})
