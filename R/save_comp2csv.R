#'@description  save normalized component list to csv
#'
#'
#'@param f_list a named list created from read_lmx_csv()
#'@param pre_fix pre_fix for the csv to be saved.
#'@md
save_comp2csv <- function(f_list, pre_fix = "~/Downloads/"){
  if(is.null(names(f_list))){
    stop("list must be named!\n")
  }

  for(y in 1 : length(f_list)){

    temp <- f_list[[y]]$comp_fix
    temp <- do.call(rbind, temp)

    if(y == 1){
      unit <- f_list[[y]]$unit
    }else{
      unit <- f_list[[y]]$unit
      unit[[grep("Median", names(unit))]][1 :16, 2 : ncol(unit[[grep("Median", names(unit))]])] <- NA
    }

    for (i in 1 : length(unit)) {
      temp_unit <- unit[[i]]

      if(ncol(temp) > ncol(temp_unit)){
        temp_unit <- cbind(temp_unit,
                           data.frame(matrix(nrow = nrow(temp_unit), ncol = (ncol(temp) - ncol(temp_unit))))%>%
                             set_colnames(value = NA))
      }

      colnames(temp_unit)[grep("^Var", colnames(temp_unit))] <- NA
      emp_rows <- data.frame(matrix(nrow = 2 ,ncol = ncol(temp_unit)))%>%
        set_colnames(value = colnames(temp_unit))
      emp_rows[1, 1 : 2] <- str_split(names(unit)[[i]], pattern = ",")[[1]]
      emp_rows[2, ] <- colnames(emp_rows)

      emp_row <- data.frame(matrix(ncol = ncol(temp_unit)))%>%
        set_colnames(value = colnames(temp_unit))

      colnames(temp) <- colnames(temp_unit)
      temp_unit <- rbind(emp_rows, temp_unit)

      temp <- rbind(temp, temp_unit)
      temp <- rbind(temp, emp_row)
    }

    f_list[[y]]$comp_last <- f_list[[y]]$comp_last[ ,1 : ncol(temp)]
    colnames(temp) <- colnames(f_list[[y]]$comp_last)
    temp <- rbind(temp, f_list[[y]]$comp_last)
    #colnames(temp) <- NULL

    f_name <- paste0(pre_fix,
                     gsub("(^.*/|.csv)", "", names(f_list)[[y]]),
                     "std_aved_med_normed.csv")
    write_csv(temp, f_name, col_names = F, na = "") #2 need colname
  }
}
