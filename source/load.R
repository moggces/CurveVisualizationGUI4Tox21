
# same as the new_qhtsdata2rdata
####### dependency: plyr:join
# mapping_file: file name
load_mapping_file <- function (mapping_file)
{
  mapping_df <- read.delim(mapping_file, stringsAsFactors=FALSE)
  #s <- strsplit(mapping_df$Sample.ID, split = ";")
  #t <- strsplit(mapping_df$Tox21.ID, split = ";")
  #result <- data.frame(CAS = rep(mapping_df$CAS, sapply(s, length)), Sample.ID = unlist(s), Tox21.ID = unlist(t), stringsAsFactors=FALSE)
  #result <- join(result, subset(mapping_df, select=c(CAS,Chemical.Name,StructureID)))
  result <- subset(mapping_df, select=c(CAS, Tox21.ID, Chemical.Name,StructureID,Chemical.ID))
  return(result)
}

# a wrapper function
load_profile <- function (profile_file)
{
  #return ( read.table(paste(getwd(),  profile_file, sep=""), header = TRUE, sep = "\t", quote = '"', comment.char = "") )
  return ( read.table( profile_file, header = TRUE, sep = "\t", quote = '', check.names=FALSE, comment.char = "") )
}

####### dependency: plyr:join, load_mapping_file
# cebs_file: file name
# mapping: output of load_mapping_file
# pathway: the name you want to assign ; NULL to disable
# readout: type. limited to the one listed ; NULL to disable
load_cebs_file <- function (cebs_file, mapping, pathway="", readout='')
{
  #skip the lines staring with #
  skipl <- length(grep("^#", readLines(cebs_file)))
  cebs_df <- read.delim(cebs_file, quote = "", skip=skipl, stringsAsFactors=FALSE)
  
  base_cols <- c('CAS', 'Chemical.Name','StructureID','Chemical.ID')
  result <- join(cebs_df[, ! colnames(cebs_df) %in% base_cols, drop=FALSE], 
                 mapping[, colnames(mapping) %in% c('Tox21.ID', base_cols), drop=FALSE],  by="Tox21.ID",  match = "first")
  
  #create pathway column 
  if (! is.null(pathway)) result$pathway <- pathway
  
  #create Tox21AgencyID column
  result$Tox21AgencyID <- paste(result$Tox21.ID, "@", result$Cmpd_Library, sep="")
  
  #create readout column
  if (! is.null(readout)) 
  {
    result[, 'readout'] <- paste(readout, substr(result$Library, nchar(result$Library), nchar(result$Library) ), sep=".")
    result[, 'tempid'] <- paste(result$Tox21.ID, result$Library)
    result <- result[order(result$tempid, result$readout),]
    replicate <- names(table(result$tempid)[table(result$tempid) > 1])
    times <- as.vector(t(table(result[result$tempid %in% replicate, "tempid"])))
    result[result$tempid %in% replicate,"readout"] <- paste(result[result$tempid %in% replicate,]$readout,  "_", unlist(sapply(times, seq)), sep="")
    result[, "tempid"] <- NULL
    
    #    result$readout <- paste(readout, substr(result$Library, nchar(result$Library), nchar(result$Library) ), sep=".")
    #    result <- result[order(result$Tox21AgencyID, result$readout),]
    #    replicate <- names(table(result$Tox21AgencyID)[table(result$Tox21AgencyID) > 3])
    #    times <- as.vector(t(table(result[result$Tox21AgencyID %in% replicate, c("Tox21AgencyID", "readout")])))
    #    result[result$Tox21AgencyID %in% replicate,]$readout <- paste(result[result$Tox21AgencyID %in% replicate,]$readout,  "_", unlist(sapply(times, seq)), sep="")
  }
  
  # get only the readout name (remove the batch info)
  readout_n <- strsplit(result$readout, ".", fixed=TRUE)[[1]][1]
  
  ### add more lines for the autofluro data
  if ( grepl("cell|medi", readout_n) ) 
  {
    app <- lapply(grep("_[0-9]$", mapping$Tox21.ID, value=TRUE), function (x)
    {
      x1 <- sub("_[0-9]$", "", x)
      l <- grep(x1, result$Tox21.ID, fixed=TRUE)
      if (length(x) > 1) l <- sample(l,1)
      l <- result[l, ]
      l[, "Tox21.ID"] <- x
      return(l)
    }
    )
    app <- as.data.frame(do.call("rbind", app))
    result <- rbind(result, app)
  }
  
  return(list(plate_data=result, readout=readout_n, pathway=unique(result$pathway)))
  
}