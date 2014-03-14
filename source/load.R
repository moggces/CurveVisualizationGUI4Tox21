
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

####### dependency: plyr:join, load_mapping_file
# cebs_file: file name
# mapping: output of load_mapping_file
# pathway: the name you want to assign
# readout: type. limited to the one listed
load_cebs_file <- function (cebs_file, mapping, pathway="", readout=c('ratio', 'ch2', 'ch1', 'luc', 'via', '653', '657', '100'))
{
  cebs_df <- read.delim(cebs_file, quote = "", stringsAsFactors=FALSE)
  result <- join(subset(cebs_df, select=-CAS), subset(mapping, select=c(CAS, Tox21.ID, Chemical.Name,StructureID,Chemical.ID)), by="Tox21.ID",  match = "first")
  result$readout <- paste(readout, substr(result$Library, nchar(result$Library), nchar(result$Library) ), sep=".")
  result$pathway <- pathway
  result$Tox21AgencyID <- paste(result$Tox21.ID, "@", result$Cmpd_Library, sep="")

#   replicate <- names(table(result$Tox21AgencyID)[table(result$Tox21AgencyID) > 3])
#   for (name in replicate)
#   {
#     for (name2 in paste(readout, seq(1,3), sep="."))
#     {
#       times <- sum(result$Tox21AgencyID == name & result$readout == name2)
#       result[result$Tox21AgencyID == name & result$readout == name2, ]$readout <-  paste(name2, seq(1,times), sep="_")
#     }
#   }
   result <- result[order(result$Tox21AgencyID, result$readout),]
   replicate <- names(table(result$Tox21AgencyID)[table(result$Tox21AgencyID) > 3])
   times <- as.vector(t(table(result[result$Tox21AgencyID %in% replicate, c("Tox21AgencyID", "readout")])))
   result[result$Tox21AgencyID %in% replicate,]$readout <- paste(result[result$Tox21AgencyID %in% replicate,]$readout,  "_", unlist(sapply(times, seq)), sep="")
  return(result)
  
}


