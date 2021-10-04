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
    colnames(temp) <- paste(colnames(temp), x$Plate.ID[1], sep = "_")
    temp
  })

  row_data <- rowdata_list[[1]]
  #colnames(row_data)[1] <- "Analyt"
  for (i in 2 : length(rowdata_list)) {
    row_data <- cbind(row_data, rowdata_list[[i]][ ,-c(1,2)]) #Analyt Assay
  }

  colnames(row_data) <- make.names(colnames(row_data), unique = T)

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
