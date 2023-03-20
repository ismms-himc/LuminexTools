#'read_lmx_csv
#'@description  read lmx raw file to list of dataframes of each components
#'
#'
#'@param f a string file path
#'@param plex "384" or "96"
#'@return a list each elemment is a component from the raw csv.
#'@md

read_lmx_csv <- function(f, plex){

  check_row <- c(
    "DataType:,Median",
    "DataType:,Net MFI",
    "DataType:,Count",
    "DataType:,Avg Net MFI",
    "DataType:,Mean",
    "DataType:,%CV",
    "DataType:,Peak",
    "DataType:,Std Dev",
    "DataType:,Trimmed Count",
    "DataType:,Trimmed Mean",
    "DataType:,Trimmed % CV of Microspheres",
    "DataType:,Trimmed Peak",
    "DataType:,Trimmed Standard Deviation",
    "DataType:,Units",
    "DataType:,Per Bead Count",
    "DataType:,Acquisition Time",
    "DataType:,Dilution Factor",
    "DataType:,Analysis Types",
    "DataType:,Audit Logs",
    "DataType:,Warnings/Errors",
    "CRC"
  )

  unit <- scan(f, what = "character", sep = "\t")

  #-----
  comp_idx <- numeric()
  for(i in 1:length(check_row)){
    try(comp_idx[i] <- grep(check_row[i], unit))
  }
  comp_idx <- comp_idx[!is.na(comp_idx)]

  comp_last <- unit[comp_idx[length(comp_idx)]: (comp_idx[length(comp_idx)] + 1)]

  #---
  comp_name <- sapply(comp_idx[1 : (length(comp_idx) - 1)], function(x){
    unit[x]
  })

  comp_idx <- sapply(1 : (length(comp_idx) - 1), function(x){
    seq((comp_idx[x] + 1), (comp_idx[x + 1] - 1))
  })
  #-----

  #fix_comp_idx <- check_fix_row[length(check_fix_row)]
  #comp_fix <- unit[1 : grep(fix_comp_idx, unit)]
  if(plex == "384"){
    comp_fix <- list(
      unit[1 : 3],
      ",,,",
      unit[5 : 26],
      ",,,",
      unit[28 : 34],
      ",,,",
      unit[36 : 45],
      ",,,",
      ",,,",
      unit[48],
      ",,,",
      unit[50],
      ",,,"
    )
  }else{
    comp_fix <- list(
      unit[1 : 3],
      ",,,",
      unit[4 : 24],
      ",,,",
      unit[25 : 29],
      ",,,",
      unit[30 : 36],
      ",,,",
      ",,,",
      unit[37],
      ",,,",
      unit[38],
      ",,,"
    )
  }


  fix_c_name <- max(sapply(comp_fix, function(f){
    str_count(f, pattern = ",")
  })%>% unlist())

  unit <- lapply(comp_idx, function(x){
    unit[x]
  })

  names(unit) <- comp_name


  c_name <- max(sapply(unit, function(x){
    str_count(x[1], pattern = ",")
  }))
  c_name <- paste("c", 1 : (max(c(c_name, fix_c_name)) + 2), sep = "_")

  # -to dataframe

  comp_last <- data.frame(dat = comp_last)%>%
    separate(col = dat, into = c_name, sep = ",")

  for (i in 1 : length(comp_fix)) {
    comp_fix[[i]] <-  data.frame(dat = comp_fix[[i]])%>%
      separate(col = dat, into = c_name, sep = ",")
  }

  unit <- lapply(unit, function(x){
    cname <- unlist(strsplit(x[1], split = ","))
    temp <- data.frame(dat = x[-1])%>%
      separate(col = dat, into = c_name, sep = ",")%>%
      mutate(c_1 = paste(c_1, c_2, sep = ","))

    if(temp$c_2[1] == "A1)" & (nrow(temp) != 0)){
      temp <- temp%>%select(-c_2)
    }

    temp <- temp[temp[ ,1] != ",", ] # random string from readin

    colnames(temp) <- cname
    temp
  })


  list("comp_fix" = comp_fix,
       "unit" = unit,
       "comp_last" = comp_last)
}
