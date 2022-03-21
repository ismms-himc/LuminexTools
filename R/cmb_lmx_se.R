#'@description combind a list of luminex summarizedexperiment
#'@param se_list list of npx summarizedexperiment obj
#'rowData will be cbind, common colData
#'@return combined summarizedexperiment obj
#'@export
#'@md
#'
cmb_lmx_se <- function(se_list){

  rowdata_list <- lapply(se_list, function(x){
    temp <- data.frame(x@elementMetadata@listData)
    Plate.ID <- gsub("(^.*/|\\.Detail.xls)",  "",x$File[1])
    colnames(temp) <- paste(colnames(temp), Plate.ID, sep = "_")
    temp
  })

  row_data <- rowdata_list[[1]]
  #colnames(row_data)[1] <- "Analyt"
  for (i in 2 : length(rowdata_list)) {
    row_data <- cbind(row_data, rowdata_list[[i]][ ,-c(1,2)]) #Analyt Assay
  }

  colnames(row_data) <- make.names(colnames(row_data), unique = T)
  row_data$LOD <- rowMeans(row_data[ , grepl("LOD_.*$", colnames(row_data))], na.rm = T)
  row_data$HOD <- rowMeans(row_data[ , grepl("HOD_.*$", colnames(row_data))], na.rm = T)

  # remove LOD from combined se object does not make senes to keep it
  se_list <- lapply(se_list, function(x){
    x@elementMetadata@listData$LOD <- x@elementMetadata@listData$HOD <- NULL
    x@elementMetadata@listData <- data.frame()
    x
  })
  temp <- se_list[[1]]
  for (i in 2 : length(se_list)) {
    temp <- cbind(temp, se_list[[i]])
  }
  rowData(temp) <- row_data
  return(temp)
}
