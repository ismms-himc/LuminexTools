#'@description  read lmx raw file to list of dataframes of each components
#'
#'
#'@param f a string file path
#'@param plex "384" or "96"
#'@return a list each elemment is a component from the raw csv.
#'@export
#'@md

lmx_norm_from_raw <- function(f_list, plex, pre_fix, map_csv, morm_method = "median", ref_sample = "QC"){
  if(rlang::is_null(names(f_list))){
    names(f_list) <- f_list
  }

  #readin csv to components list
  f_list <- lapply(f_list, read_lmx_csv, plex = plex)

  #plot beads count heatmap
  count_qc_plot(f_list, pre_fix = pre_fix)


  f_list <- norm_replace_fun(f_list, morm_method = morm_method, ref_sample = ref_sample, map_csv)

  save_comp2csv(f_list,pre_fix = pre_fix)
}
