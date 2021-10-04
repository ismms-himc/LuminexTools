#'bdg_norm_lmx (Reference sample based normalization)
#'@description normalize by adjust samples' value per plate to make that
#'             the reference sample among all plates are the same. option to
#'             pick using median, max, mean
#'
#'@param bridge.ls list of bridge.summarizedexperiment objs using pull_bdg()
#'@param data.ls list of summarizedexperiment objs using read_npx().
#'               normalized assay slot will be added to each of the object in the list.
#'@param between.plate.method method to set the inter plate reference value of bridging samples
#'       using max can garantee no negative value of the normalized data
#'@param round_digits digit kept 
#'@param from_assay select assay slot to be normalized
#'@param save_assay name the assay slot for normalized data
#'@export
#'@md
#'
bdg_norm_lmx <- function(bridge.ls, data.ls, between.plate.method = "max", round_digits = 3,
                         from_assay = "default", save_assay = "normed"){
  
  if(is.null(names(data.ls)) | (sum(names(bridge.ls) != names(data.ls)) != 0)){
    stop("list of bridges and data has to be named list!\n")
  }
  
  assay.from <- paste("data", from_assay, sep = "_")
  #assay.from.mfi <- paste("mfi", from_assay, sep = "_")
  
  bridge <- cmb_lmx_se(bridge.ls)
  
  #  count
  bridge.plate.mean <- cbind.data.frame(File = bridge$File,
                                        t(bridge@assays@data[[assay.from]] %>% log10))%>% 
    group_by(File)%>%
    summarize_all(.fun = mean, na.rm = T)
  
  
  bridge.median <- bridge.plate.mean%>%
    dplyr::select(-File)%>%
    summarize_all(.funs = between.plate.method, na.rm = T)
  
  #  mfi
  bridge.plate.mean.mfi <- cbind.data.frame(File = bridge$File,
                                            t(bridge@assays@data$mfi_default %>% log10))%>%
    group_by(File)%>%
    summarize_all(.fun = mean, na.rm = T)
  
  
  bridge.median.mfi <- bridge.plate.mean.mfi%>%
    dplyr::select(-File)%>%
    summarize_all(.funs = between.plate.method, na.rm = T)  
  
  # update adjust factor   
  query <- bridge.plate.mean$File
  names(query) <- query
  
  # adjust factor count
  bridge.adj <- lapply(query, function(x){
    bridge.plate.mean[bridge.plate.mean$File == x ,-1] - unlist(bridge.median)
  })
  
  # adjust factor mfi
  bridge.adj.mfi <- lapply(query, function(x){
    bridge.plate.mean.mfi[bridge.plate.mean.mfi$File == x ,-1] - unlist(bridge.median.mfi)
  })
  
  
  save.assay <- paste("data", save_assay, sep = "_")
  #save.assay.mfi <- paste("mfi", save_assay, sep = "_")
  data.ls <- lapply(query, function(x){
    temp <- data.ls[[x]]
    temp@assays@data[[save.assay]] <- 10^(log10(temp@assays@data[[assay.from]]) - unlist(bridge.adj[[x]]))%>%round(digits = round_digits)
    temp@assays@data$mfi_normed <- 10^(log10(temp@assays@data$mfi_default) - unlist(bridge.adj.mfi[[x]]))%>%round(digits = round_digits)
    #temp@elementMetadata$LOD_normed <- 10^(temp@elementMetadata$LOD - unlist(bridge.adj[[x]]))
    #temp@elementMetadata$HOD_normed <- 10^(temp@elementMetadata$HOD - unlist(bridge.adj[[x]]))
    temp
  })
  
  data.ls
}
