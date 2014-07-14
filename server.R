
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(plyr)
library(reshape2)
library(ggplot2)
#library(googleVis)

source( "./source/get.R",  local=TRUE)
source("./source/load.R",  local=TRUE)
mapping_file <- './data/tox21_id_12807_cas_structure_4v4_080813_v4.txt'
mapping_df <- read.delim(mapping_file, stringsAsFactors=FALSE)



################

shinyServer(function(input, output) {
  
  data_chemical <- reactive({
    result <- NULL
    ids <- input$cmpds
    result <- unlist(strsplit(ids, "\n", perl=TRUE))
    return(result)
  })
  
  pathway_data <- reactive({
    

    chemicals <- data_chemical()
    pathways <- input$paths
    
    result <- get_qhts_data_wrap(chemicals, pathways=pathways)
    

    return(result)
  })
  
  pathway_add_data <- reactive({
    
    chemicals <- data_chemical()
    result <- pathway_data()
    
    options <- input$opts
    options <- grep("^(?!run[1-3]$)", options,  perl=TRUE, value = TRUE)
    pathways <- input$paths 
    
    result <- rbind.fill(result, get_qhts_data_wrap(chemicals, pathways=pathways, options=options))
  
    return(result)
    
  })
  
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
      
    result <- get_relevant_cmpd_library(result)
    
    return(result)
  })
  
  data_melter <- reactive ({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')
    
    qhts <- data_filter()
    
    result <- get_melt_data(qhts, resp_type=unique(c('raw', plot_options)))
    if (nrow(result) > 0)
    {
      result[, "display_name"] <- paste(result$CAS, "|\n", result$Chemical.Name, "|\n", result$Tox21AgencyID, sep="")
      #result[, "display_name"] <- paste(result$CAS, "|\n", result$Tox21AgencyID, sep="")
      result <- result[order(result$display_name),]
    }
    
    return(result)
  })
  
   
  getVarHeight <- reactive({
    qhts <- data_filter()
    qhts[, "display_name"] <- paste(qhts$CAS, "|\n", qhts$Tox21AgencyID, sep="")
    nrow <- length(unique(qhts$display_name))
    mode <- input$mode
    heightpx <- input$heightpx
    if (mode == 'overlay')
    {
      return(nrow * heightpx) # 350
    } else if (mode == 'parallel')
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
    } else if (mode == 'parallel')
    {
      return(ncol * widthpx)
    }
  })

  output$pathways <- renderUI({
    selectInput("paths", "Select pathways to show:", 
                choices  = list("ATAD5"="atad5_act_main", "p53"="p53_act_main",
                                "DNA Damage (SRF)"="srf_inh_main", "DNA Damage (DSB)"="dsb_inh_main",
                                'AhR agonism' = 'ahr_act_main',
                                'AR agonism (p-length)'='arpartial_act_main',
                                'AR agonism (f-length)'='arfull_act_main',
                                'AR antagonism (p-length)'='arpartial_inh_main',
                                'AR antagonism (f-length)'='arfull_inh_main',
                                'ER agonism (p-length)'='erpartial_act_main',
                                'ER agonism (f-length)'='erfull_act_main',
                                'ER antagonism (p-length)'='erpartial_inh_main',
                                'ER antagonism (f-length)'='erfull_inh_main',
                                'GR agonism' = 'gr_act_main',
                                'GR antagonism' = 'gr_inh_main',
                                'FXR agonism'='fxr_act_main',
                                "FXR antagonism"="fxr_inh_main",
                                'PPARG agonism' = 'pparg_act_main',
                                'PPARG antagonism' = 'pparg_inh_main',
                                'PPARD agonism' = 'ppard_act_main',
                                'PPARD antagonism' = 'ppard_inh_main',
                                'TR agonism' = 'tr_act_main',
                                "TR antagonism"="tr_inh_main",
                                "aromatase antagonism"="aromatase_inh_main",
                                'mitochondrial toxicity'='mito_inh_main',
                                'Nrf2/ARE' = 'are_act_main',
                                'HSE'='hse_act_main',
                                'autofluorescence (Hek293)'='autofluo_hek293_main'
                                ), 
                multiple = TRUE)
  })
  
  output$options <- renderUI({
    selectInput("opts", "Select readout options:", 
                choices  = list( "run1"="run1","run2"="run2", "run3"="run3",
                                 "cytotoxicity"="via", 
                                 "ch2 (in bla)"="ch2", 
                                 "ch1 (in bla)"="ch1",
                                 "medium blue (in autofluo)"="medi_blue",
                                 "cell blue (in autofluo)"="cell_blue"),
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
  

  
  output$plot <- renderPlot({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')
    
    result <- data_melter()
    p <- get_plot(result, mode=mode, plot_options=plot_options, fontsize=20, pointsize=3)
    
    if (mode == 'overlay')
    {
      p <- p  + facet_wrap(~ display_name  , ncol=2)
    } else if (mode == 'parallel')
    {
      p <- p + facet_grid(display_name ~ pathway)
    }
    print(p)
    #print(select_plot())
  }, height=getVarHeight, width=getVarWidth)


  output$qhts_data <- renderDataTable({
    data_filter()
  })
  
  output$downloadPlot <- downloadHandler(
         filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".pdf", sep="") },
         content = function(file) {
           result <- data_melter()
           mode <- input$mode
           plot_options <- input$plt_opts
           show_outlier <- input$showOutlier
           if (show_outlier) plot_options <- c(plot_options, 'mask')
      
           pdf(file)
           
           if (mode == 'overlay')
           {
             n_page <- 6
             result <- get_blank_data(result, n_page)
             nn <- unique(result$display_name)
             pages <- split(nn, ceiling(seq_along(nn)/n_page))
             lapply(names(pages), function (x) {
               sub <- result[result$display_name %in% pages[[x]],]
               p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
               p <- p  + facet_wrap(~ display_name  , ncol=2, nrow=3)
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
               p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
               p <- p + facet_grid(display_name ~ pathway)
               print(p)
             })
             
           }
           
           dev.off()
         }
  )
  
})
