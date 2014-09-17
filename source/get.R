
####### dependency: 
# id: 
get_id_type <- function (id)
{
  result <- 'unique'
  for (i in id)
  {
    if (grepl("^Tox21_[0-9]{6}$", i, perl=TRUE))
    {
      result <- 'tox21'
    } else if (grepl("^NCGC_[0-8]{8,]\\-[0-9]{2}$",i, perl=TRUE)) 
    {
      result <- 'ncgc'
    } else if ( (length(strsplit(i, "-")[[1]]) == 3 & ! grepl("[a-z|A-Z]", i, perl=TRUE) ) | grepl('NOCAS', i) )
    {
      result <- 'cas'
    } 
  }
  return(result)
}

####### dependency: get_id_type
# id: 
# cebs: 
get_qhts_data <- function (id, cebs)
{
  id_type <- get_id_type(id)
  result <- data.frame()
  
  if (id_type == 'cas')
  {
    result <- rbind(result, cebs[cebs$CAS  %in% id, ])
  } else if (id_type == 'tox21')
  {
    result <- rbind(result, cebs[ cebs$Tox21.ID %in% id, ])
  } else if (id_type == 'ncgc')
  {
    result <- rbind(result, cebs[cebs$Sample.ID %in% id, ])
  } else if (id_type == 'unique' )
  {
    result <- rbind(result, cebs[cebs$Chemical.ID %in% id, ])
  }
  
  return(result)
  
}

####### dependency: get_qhts_data
# id: 
# pathways: 
get_qhts_data_wrap <- function(chemicals, pathways, options=NULL)
{
  result <- data.frame()
  
  
  for (name in pathways)
  {
    if (is.null(options))
    {
      rda <- paste("./data/", name, ".RData", sep="")
      if (file.exists(rda)) 
      {
        load(rda) # The data frame is cebs
        result <- rbind.fill(result, get_qhts_data(chemicals, cebs))
        rm(cebs)
      }
    } else
    {
      base <- sub("_main", "", name)
      for (name2 in options)
      {
        rda <- paste("./data/", base, "_", name2, ".RData", sep="")
        if (file.exists(rda)) 
        {
          load(rda) # The data frame is cebs
          result <- rbind.fill(result, get_qhts_data(chemicals, cebs))
          rm(cebs)
        }
      }
    }
  }
  return(result)
  
}



###http://stackoverflow.com/questions/5177846/equivalent-of-curve-for-ggplot
####### dependency: get_qhts_data
get_melt_data <- function (qhts, resp_type=c('raw', 'curvep', 'hill', 'mask'))
{
  col_names <- colnames(qhts)
  basic_cols <- c('CAS', 'uniqueID', 'Tox21AgencyID', 'Chemical.ID', 'Chemical.Name', 'Tox21.ID', 'Sample.ID', 'StructureID', 'readout', 'pathway')
  basic_cols <- intersect(col_names, basic_cols)
  
  x_cols <- grep("conc[0-9]+", col_names, value = TRUE)
  #x_cols <- c(paste('conc', "", seq(0,14), sep=""))
  
  result <- qhts[, c(basic_cols, x_cols)]
  result <- melt(result, id.var=basic_cols, value.name='x', variable.name='concs')
  
  for (type in resp_type)
  {
    if (type == 'raw')
    {
      y_cols <- grep("resp[0-9]+", col_names,  value = TRUE)
      #y_cols <- c(paste('resp', "", seq(0,14), sep=""))
      
      temp <- qhts[, c(basic_cols, y_cols)]
      temp <- melt(temp, id.var=basic_cols, value.name='raw', variable.name='raw_resps')
      result$raw <- temp$raw
      
      
      
    } else if (type == 'curvep')
    {
      #y_cols <- c(paste('curvep_r', "", seq(0,14), sep=""))
      y_cols <- grep("curvep_r[0-9]+", col_names, value = TRUE)
      
      temp <- qhts[, c(basic_cols, y_cols)]
      temp <- melt(temp, id.var=basic_cols, value.name='curvep', variable.name='curvep_resps')
      result$curvep <- temp$curvep
      
      
      
    } else if (type == 'hill')
    {
      #y_cols <- c(paste('resp', "", seq(0,14), sep=""))
      y_cols <- grep("resp[0-9]+", col_names, value = TRUE)
      
      
      hill_resps <- lapply(1:nrow(qhts), function(x) {
        sapply(qhts[x, x_cols], function (y)
          if (! is.na(qhts[x,]$Hill.Coef) )
          {
            qhts[x,]$Zero.Activity + (qhts[x,]$Inf.Activity - qhts[x,]$Zero.Activity)/(1+10^((qhts[x,]$LogAC50-y)*qhts[x,]$Hill.Coef))
          } else
          {
            if (! is.na(y)) 
            {
              mean(unlist(qhts[x, y_cols]), na.rm=TRUE)
            } else NA
          }
        )
      }
      )
      hill_resps <- do.call("rbind", hill_resps)
      colnames(hill_resps) <- sub("conc", "hill", colnames(hill_resps))
      temp <- cbind(qhts[, basic_cols], hill_resps)
      temp <- melt(temp, id.var=basic_cols, value.name='hill', variable.name='hill_resps')
      result$hill <- temp$hill
      
      
    } else if ( type == 'mask')
    {
      #y_cols <- c(paste('resp', "", seq(0,14), sep=""))
      y_cols <- grep("resp[0-9]+", col_names, value = TRUE)
      
      mask_resps <- qhts[, y_cols]
      qhts[, 'mask'] <- qhts[, get_mask_column(qhts)]
      
      mask_resps <- lapply(1:nrow(qhts), function (x) {
        if (qhts[x, 'mask'] != '')
        {
          m <- ! as.logical(as.numeric((unlist(strsplit(qhts[x, 'mask'], " ")))))
          mask_resps[x, ][which(m)] <- NA
        } else
        {
          mask_resps[x, ] <- NA
        }
        return(mask_resps[x, ])
      }
      )
      mask_resps <- as.data.frame(do.call("rbind", mask_resps))
#       mask_resps <- qhts[, y_cols]
#       for ( x in 1:nrow(qhts) )
#       {
#         if( ! is.null (qhts[x,]$Mask.Flags) )
#         {
#           if ( qhts[x,]$Mask.Flags != ""  )
#           {
#             m <- ! as.logical(as.numeric((unlist(strsplit(qhts[x,]$Mask.Flags, " ")))))
#             #mask_resps[x, ][! is.na(mask_resps[x, ])][which(m)] <- NA
#             mask_resps[x, ][which(m)] <- NA
#           } 
#         } else if (! is.null (qhts[x,]$curvep_mask) )
#         {
#           if (qhts[x,]$curvep_mask != "")
#           {
#             m <- ! as.logical(as.numeric((unlist(strsplit(qhts[x,]$curvep_mask, " ")))))
#             #mask_resps[x, ][! is.na(mask_resps[x, ])][which(m)] <- NA
#             mask_resps[x, ][which(m)] <- NA
#           }
#         } else
#         {
#           mask_resps[x, ] <- NA
#         }
#       }
      colnames(mask_resps) <- sub("resp", "mask", colnames(mask_resps))
      temp <- cbind(qhts[, basic_cols], mask_resps)
      temp <- melt(temp, id.var=basic_cols, value.name='mask', variable.name='mask_resps')
      result$mask <- temp$mask
    }
  }
  
  return(result)
  
}



## ggplot(temp, aes(x=x, y=raw, color=readout)) + geom_line() + facet (uniqueID ~ pathway)


####### dependency: data_chemical(reactive, vector); get_id_type
get_mapping_data <- function (input, mapping)
{
  id_type <- get_id_type(input)
  id_type <- switch(id_type, "cas"="CAS", "tox21"="Tox21.ID", "ncgc"="Sample.ID")
  result <- data.frame(input)
  colnames(result) <- id_type
  result <- join(result, mapping)
  #result <- subset(result, select=c(CAS, Tox21.ID, Sample.ID, Chemical.Name,StructureID))
  #result <- subset(result, select=c(CAS, Tox21.ID, Chemical.Name,StructureID))
  return(result)
}

get_plot <- function (qhts_melt, mode=c('parallel', 'overlay'), plot_options=plot_options, fontsize=20, pointsize=3)
{
  if (mode == 'parallel')
  {
    p <- ggplot(qhts_melt, aes(x=x, y=raw, color=readout)) + geom_point(size=pointsize)+
      theme(text = element_text(size=fontsize) ) + scale_x_continuous('log10(conc(M))') + 
      scale_y_continuous('resp (%)')
    if (sum (plot_options %in% 'raw') == 1) p <- p + geom_line()
    if (sum (plot_options %in% 'curvep') == 1) p <- p + geom_line(aes(x=x, y=curvep, color=readout),  linetype=5)
    if (sum (plot_options %in% 'hill') == 1) p <- p + geom_line(aes(x=x, y=hill, color=readout), linetype=6)
    if (sum (plot_options %in% 'mask') == 1) p <- p + geom_point(aes(x=x, y=mask, color=readout),shape = 4, size=pointsize*2)
    #p <- p + facet_grid(display_name ~ pathway)
    
  } else if (mode == 'overlay')
  {
    qhts_melt$path_readout <- paste(qhts_melt$pathway, "|\n", qhts_melt$readout, sep="")
    p <- ggplot(qhts_melt, aes(x=x, y=raw, color=path_readout)) +  geom_point(size=pointsize) + 
      theme(text = element_text(size=fontsize) ) + scale_x_continuous('log10(conc(M))') + 
      scale_y_continuous('resp (%)')
    if (sum (plot_options %in% 'raw') == 1) p <- p + geom_line()
    if (sum (plot_options %in% 'curvep') == 1) p <- p + geom_line(aes(x=x, y=curvep, color=path_readout),  linetype=5)
    if (sum (plot_options %in% 'hill') == 1) p <- p + geom_line(aes(x=x, y=hill, color=path_readout), linetype=6)
    if (sum (plot_options %in% 'mask') == 1) p <- p + geom_point(aes(x=x, y=mask, color=path_readout),shape = 4, size=pointsize*2)
    #p <- p  + facet_wrap(~ display_name  , ncol=2)
    #p <- p + facet_grid(display_name ~. )
  }
  return(p)
}


get_blank_data <- function (qhts_melt, n_page)
{
  result <- qhts_melt
  nn <- unique(qhts_melt$display_name)
  base <- subset(qhts_melt, display_name == nn[1])
  if (! is.null(base$raw) ) base$raw <- 0
  if (! is.null(base$curvep) ) base$curvep <- 0
  if (! is.null(base$hill) ) base$hill <- 0
  if (! is.null(base$mask) ) base$mask <- 0
  
  blank_n <- n_page - (length(nn) %% n_page)
  for (i in 1:blank_n)
  {
    base$display_name <- paste('zzzzz',".", i, sep="")
    result <- rbind(result, base)
  }
  return(result)
  
}

get_relevant_cmpd_library <- function (qhts)
{
  autod <- qhts[grep("cell|medi", qhts$readout), ]
  otherd <- qhts[grep("cell|medi", qhts$readout, invert=TRUE), ]
  result <- data.frame()
  
  if (nrow(autod) > 0)
  {
    # the situation where only auto data are available ...
    if (nrow(otherd) > 0)
    {
      result <- lapply(unique(otherd$Tox21AgencyID), function (x)
        {
          v <- strsplit(x, "@")[[1]]
          id <- v[1]
          cmpd_library <- v[2]
          ind <- grep(id, autod$Tox21.ID, fixed=TRUE)
          if (length(ind) > 1) ind <- sample(ind,1)
        
          s <- split(autod, autod$readout)
          s <- lapply(s, function (x)
            { 
              ind <- grep(id, x$Tox21.ID, fixed=TRUE)
              if (length(ind) > 1) ind <- sample(ind,1)
              return(x[ind,])
            }
          )
        
          result <- as.data.frame(do.call("rbind", s))
          result[, "Tox21AgencyID"] <- x
          #result[, "Cmpd_Library"] <- cmpd_library
          return(result)
        }
      )
      
      result <- as.data.frame(do.call("rbind.fill", result))
    } else
    {
      result <- autod
    }
   
  }
  result <- rbind.fill(otherd, result)
  return(result)
}

## different from the CurveVisualizationGUI
get_mask_column <- function (qhts)
{
  result <- 'Mask.Flags'
#   if ( nrow(qhts) == sum(qhts$mask == '') )
#   {
#     if (! is.null(qhts$curvep_mask))
#     {
#       result <- 'curvep_mask'
#     } 
#   }
  return(result)
}