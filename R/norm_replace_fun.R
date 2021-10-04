#'@description  normalized MFI component by selected reference sample,
#'              position-wise average the STD from each plates,
#'              only the 1st plate will keep STD as the averaged value,
#'              (this is adapt to the luminex analyser's format).
#'
#'
#'@param f_list a named list created from read_lmx_csv()
#'@param morm_method normalization method, such as mean, median...
#'@param map_csv csv files stored in inst/extra, or user
#'@param ref_sample reference sample for between plates normalization
#'@return a list with the MFI component ref-sample normalized.
#'@md
norm_replace_fun <- function(f_list, morm_method = "median", ref_sample = "QC", map_csv){
  if(is.null(names(f_list))){
    stop("list must be named!\n")
  }
  query <- names(f_list)
  names(query) <- query

  f_list_median <- lapply(query, function(x){
    temp <- f_list[[x]]
    temp <- temp$unit[[grep("Median", names(temp$unit))]]
    temp[ ,1 : which(colnames(temp) == "Total Events")]
  })



  # std3 as reference normalize

  f_list_std3 <- lapply(f_list_median, function(x){
    x%>%
      filter(grepl(ref_sample, Sample))%>%
      select(-Location, -Sample) %>%
      mutate_all(as.numeric)%>%
      summarize_all(.funs = morm_method, na.rm = T)%>%
      mutate_all(log10)
  })

  all_f_list_std3_mean <- colMaxs(as.matrix(do.call(rbind, f_list_std3)))

  adj_factor <- lapply(f_list_std3, function(x){
    x - all_f_list_std3_mean
  })


  f_list_median <- lapply(query, function(x){
    temp_loc_samp <- f_list_median[[x]]%>%
      select(Location, Sample)
    temp_last <- f_list_median[[x]]%>%
      select("Total Events")
    temp <- f_list_median[[x]]%>%
      select(-Location, -Sample)%>%
      mutate_all(as.numeric)
    cbind(temp_loc_samp,
          t(apply(temp, 1, function(y){
            round(10^(log10(y+1) - unlist(adj_factor[[x]]))-1, 1)
          }))[ ,-ncol(temp)],
          temp_last
    )
  })


  # extract std and sample

  f_list_median_std <- lapply(f_list_median, function(x){
    x%>%
      filter(grepl("(Blank|STD)", ignore.case = T, Sample))
  })

  f_list_median_sample <- lapply(f_list_median, function(x){
    x%>%
      filter(!grepl("(Blank|STD)", ignore.case = T, Sample))
  })

  loc_samp <- f_list_median_std[[1]]%>%
    select(Location, Sample)


  #- averaged std to replace back

  ave_std <- loc_samp%>%
    select(Location, Sample)%>%
    left_join(
      do.call(rbind, f_list_median_std) %>%
        #rownames_to_column( var = "f_name")%>%
        group_by(Location, Sample) %>%
        mutate_all(as.numeric)%>%
        summarize_all(.funs = function(x) log10(x+1))%>%
        summarize_all(mean, na.rm = T)%>%
        group_by(Location, Sample) %>%
        summarize_all(.funs = function(x) round(10^x-1, 1))
    )



  # normalized to replace median
  f_list_median <- lapply(f_list_median_sample, function(x){
    map_csv %>%
      left_join(
        rbind(ave_std, x)
      )

  })



  f_list <- lapply(query, function(x){
    temp <- f_list[[x]]
    temp$unit[[grep("Median", names(temp$unit))]] <-  f_list_median[[x]]
    temp
  })

  return(f_list)
}
