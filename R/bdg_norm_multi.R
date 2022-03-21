#'bdg_norm_multi (Reference sample based normalization)
#'@description normalize by adjust samples' value per plate to make that
#'             the reference sample among all plates are the same. option to
#'             pick using median, max, mean. Adjusted factor is computed on log10
#'             linear scale, normalization is done on log10 scale then transfered back to
#'             the original scale as raw.
#'
#'@param bridge.str vector of common strings to identify sets of bridging samples
#'@param data.ls list of summarizedexperiment objs using read_lmx().
#'               normalized assay slot will be added to each of the object in the list.
#'@param between.plate.method method to set the inter plate reference value of bridging samples
#'       using max can garantee no negative value of the normalized data
#'@param from_assay select assay slot to be normalized
#'@param save_assay name the assay slot for normalized data
#'@export
#'@md
#'
bdg_norm_multi <- function(bridge.str, data.ls, between.plate.method = "median",
                           from_assay = "data_default", save_assay = "data_default_normed"){
  if(
    sapply(bridge.str, function(x){
      lapply(data.ls, function(y){
        sum(is.na(grepl(x, y$Sample)))
      })
    })%>%unlist()%>%sum() != 0
  ){
    stop("not all bridging samples exist in each plate!")
  }

  names(bridge.str) <- bridge.str
  adj.ls <- lapply(bridge.str, function(x){
    bridge <- pull_bdg(data.ls, pattern = x, fields = "Sample")%>%
      cmb_lmx_se()

    if(is.na(bridge@assays@data[[from_assay]]) %>% sum() != 0){
      stop("one or more analyte has NA values in bridging sample, can not perform normalization! Select other Bridging sample to try again.")
    }
    #  count
    bridge.plate.mean <- cbind.data.frame(File = bridge$File,
                                          t(bridge@assays@data[[from_assay]]))%>%
      group_by(File)%>%
      summarize_all(.fun = mean, na.rm = T)


    bridge.median <- bridge.plate.mean%>%
      dplyr::select(-File)%>%
      summarize_all(.funs = between.plate.method, na.rm = T)


    # update adjust factor
    query <- bridge.plate.mean$File
    names(query) <- query

    # adjust factor count
    bridge.adj <- lapply(query, function(x){
      bridge.plate.mean[bridge.plate.mean$File == x ,-1] - unlist(bridge.median)
    })

    return(bridge.adj)
  })

  adj.mean <- adj.ls[[names(adj.ls)[1]]]
  if(length(adj.ls) > 1){
    for (i in names(adj.mean)) {
      for (j in 2 : length(adj.ls)) {
        adj.mean[[i]] <- adj.mean[[i]] + adj.ls[[j]][[i]]
      }
      adj.mean[[i]] <- adj.mean[[i]]/length(adj.ls)
    }
  }

  query <- names(adj.mean)
  names(query) <- query

  data.ls <- lapply(query, function(x){
    temp <- data.ls[[x]]
    temp@assays@data[[save_assay]] <- 10^(log10(temp@assays@data[[from_assay]]) - unlist(adj.mean[[x]]))%>%round(5)
    temp
  })

  data.ls

}
