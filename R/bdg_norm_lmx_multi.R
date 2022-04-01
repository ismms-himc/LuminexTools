#'bdg_norm_lmx_multi (Reference sample based normalization)
#'@description normalize by adjust samples' value per plate to make that
#'             the reference sample among all plates are the same. option to
#'             pick using median, max, mean
#'
#'@param bridge.str vector of common strings to identify sets of bridging samples
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
bdg_norm_lmx_multi <- function(bridge.str, data.ls, between.plate.method = "median", round_digits = 3,
                         from_assay = "default", save_assay = "normed"){

  if(is.null(names(data.ls))){
    stop("list of data has to be named list!\n")
  }

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
    #  count
    bridge.plate.mean <- cbind.data.frame(File = bridge$File,
                                          t(bridge@assays@data[[from_assay]] %>% log10))%>%
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
        #adj.mean[[i]] <- adj.mean[[i]] + adj.ls[[j]][[i]]
        adj.mean[[i]] <- mapply(function(x, y) sum(x, y, na.rm = T), x = adj.mean[[i]], y = adj.ls[[j]][[i]])
      }
      adj.mean[[i]] <- adj.mean[[i]]/length(adj.ls)
    }
  }


  adj.ls.mfi <- lapply(bridge.str, function(x){
    bridge <- pull_bdg(data.ls, pattern = x, fields = "Sample")%>%
      cmb_lmx_se()
    #  mfi
    bridge.plate.mean.mfi <- cbind.data.frame(File = bridge$File,
                                              t(bridge@assays@data$mfi_default %>% log10))%>%
      group_by(File)%>%
      summarize_all(.fun = mean, na.rm = T)


    bridge.median.mfi <- bridge.plate.mean.mfi%>%
      dplyr::select(-File)%>%
      summarize_all(.funs = between.plate.method, na.rm = T)


    # update adjust factor
    query <- bridge.plate.mean.mfi$File
    names(query) <- query

    # adjust factor mfi
    bridge.adj.mfi <- lapply(query, function(x){
      bridge.plate.mean.mfi[bridge.plate.mean.mfi$File == x ,-1] - unlist(bridge.median.mfi)
    })

    return(bridge.adj.mfi)
  })

  adj.mean.mfi <- adj.ls.mfi[[names(adj.ls.mfi)[1]]]
  if(length(adj.ls.mfi) > 1){
    for (i in names(adj.mean.mfi)) {
      for (j in 2 : length(adj.ls.mfi)) {
        #adj.mean.mfi[[i]] <- adj.mean.mfi[[i]] + adj.ls.mfi[[j]][[i]]
        adj.mean.mfi[[i]] <- mapply(function(x, y) sum(x, y, na.rm = T), x = adj.mean.mfi[[i]], y = adj.ls.mfi[[j]][[i]])
      }
      adj.mean.mfi[[i]] <- adj.mean.mfi[[i]]/length(adj.ls.mfi)
    }
  }


  query <- names(adj.mean.mfi)
  names(query) <- query
  #save.assay.mfi <- paste("mfi", save_assay, sep = "_")
  data.ls <- lapply(query, function(x){
    temp <- data.ls[[x]]
    temp@assays@data[[save_assay]] <- 10^(log10(temp@assays@data[[from_assay]]) - unlist(adj.mean[[x]]))%>%round(digits = round_digits)
    temp@assays@data$mfi_normed <- 10^(log10(temp@assays@data$mfi_default) - unlist(adj.mean.mfi[[x]]))%>%round(digits = round_digits)
    #temp@elementMetadata$LOD_normed <- 10^(temp@elementMetadata$LOD - unlist(bridge.adj[[x]]))
    #temp@elementMetadata$HOD_normed <- 10^(temp@elementMetadata$HOD - unlist(bridge.adj[[x]]))
    temp
  })

  data.ls
}
