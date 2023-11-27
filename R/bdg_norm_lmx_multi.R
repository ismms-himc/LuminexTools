#'bdg_norm_lmx_multi (Reference sample based normalization)
#'@description normalize by adjust samples' value per plate to make that
#'             the reference sample among all plates are the same. option to
#'             pick using median, max, mean
#'
#'@param bridge.str vector of common strings to identify sets of bridging samples
#'@param data.ls list of summarizedexperiment objs using read_npx().
#'               normalized assay slot will be added to each of the object in the list.
#'@param round_digits digit kept
#'@param from_assay select assay slot to be normalized
#'@param save_assay name the assay slot for normalized data
#'@export
#'@md
#'
bdg_norm_lmx_multi <- function(bridge.str, data.ls, round_digits = 3, between.plate.method = "mean", ref_batch = NULL,
                         from_assay = "default", save_assay = "normed"){

  #from_assay <- paste("data", from.assay, sep = "_")

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

  if(!is_null(ref_batch) & ( !isTRUE(ref_batch %in% names(data.ls)) | (length(ref_batch) > 1))){
    stop("set a single ref_batch from names(data.ls)!")
  }

  #check avaliable ref
  actual_n_ref <- sapply(rownames(data.ls[[1]]), function(x){
    analyte_eval <- lapply(data.ls, function(y){
      y@assays@data[[from_assay]][x, match(bridge.str, y$Sample)]
    })%>%
      do.call(what = "rbind")%>%
      colSums()
    sum(is.finite(analyte_eval))
  })

  names(bridge.str) <- bridge.str

  adj.ls <- lapply(bridge.str, function(x){
    bridge <- pull_bdg(data.ls, pattern = x, fields = "Sample")%>%
      cmb_lmx_se()
    #  count

    bridge.plate.mean <- cbind.data.frame(File = bridge$File,
                                          t(bridge@assays@data[[from_assay]] %>% log10))%>%
      group_by(File)%>%
      summarize_all(.fun = mean, na.rm = T)

    if(is_null(ref_batch)){
      bridge.median <- bridge.plate.mean%>%
        dplyr::select(-File)%>%
        summarize_all(.funs = between.plate.method, na.rm = T)
    }else{
      bridge.median <- bridge.plate.mean[bridge.plate.mean$File == ref_batch, -1]
    }


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
      adj.mean[[i]] <- adj.mean[[i]]/actual_n_ref
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

    if(is_null(ref_batch)){
      bridge.median.mfi <- bridge.plate.mean.mfi%>%
        dplyr::select(-File)%>%
        summarize_all(.funs = between.plate.method, na.rm = T)
    }else{
      bridge.median.mfi <- bridge.plate.mean.mfi[bridge.plate.mean.mfi$File == ref_batch, -1]
    }


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
    for (i in 1 : ncol(temp)) {
      temp@assays@data[[save_assay]][ , i] <- ifelse(is.finite(temp@assays@data[[from_assay]][ , i]),
                                                    yes = temp@assays@data[[save_assay]][ , i],
                                                    no = temp@assays@data[[from_assay]][ , i])
    }
    temp@assays@data$mfi_normed <- 10^(log10(temp@assays@data$mfi_default) - unlist(adj.mean.mfi[[x]]))%>%round(digits = round_digits)
    temp@elementMetadata$LOD_normed <- 10^(log10(temp@elementMetadata$LOD) - unlist(adj.mean[[x]]))
    temp@elementMetadata$HOD_normed <- 10^(log10(temp@elementMetadata$HOD) - unlist(adj.mean[[x]]))
    temp
  })

  data.ls
}

