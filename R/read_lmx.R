#'read_lmx
#'@description  read lmx raw file to create a summarizedexperiment
#'              >HOD is replaced by Inf, <LOD is replaced by -Inf;
#'              retrieve STD3, and the lowest, the hightest STD's MFI, Conc.
#'              (correlated to LLOD and HLOD and matching MFI)
#'@import dplyr
#'@export
#'@param f a string file path
#'@param lot: optional parameter to record lot information
#'@param analyt_start_row 8 for npx, 9 for quant
#'@return a summarizedexperimet object.
#'        rowData is the feature and HOD, and LOD,
#'        colData is the unique_id and subject id given by the operator, and
#'        assay file name where the sample was assayed
#'@md

read_lmx <- function(f, lot = "default", analyt_start_row){

  if(grepl("xls$", f)){
    fdata <- readxl::read_xls(f, sheet = 1, col_names = T)
  }
  if(grepl("xlsx$", f)){
    fdata <- read.xlsx(f, sheetName = "Summary", header = F)
  }

  fdata <- fdata[1 : (which(fdata[ , 1] == "Notes:")-1), ]

  coldata <- data.frame(Location = unlist(fdata[-c(1 : 4), 1]),
                        Sample = unlist(fdata[-c(1 : 4), 2]))%>%
    mutate(unique_id = paste(Location, Sample, sep = "_"),
           File = f,
           Lot = lot)%>%
    magrittr::set_rownames(value = .$unique_id)

  analyt <- fdata[3, -c(1:2)]
  analyt <- analyt[!is.na(analyt)]
  rowdata <- data.frame(Analyt = make.names(analyt),
                        Unit = as.character(fdata[4, -c(1:2)][1: length(analyt)]))%>%
    magrittr::set_rownames(value = .$Analyt)

  fdata <- data.frame(fdata[-c(1:4), -c(1:2)])
  fdata <- data.frame(fdata[, 1:length(analyt)])

  # get lod
  rowdata$LOD <- sapply(fdata, function(x){
    as.numeric(gsub("(↓|↑|<|>)", "", x[grep("<", x)][1]))
  })
  # get hod
  rowdata$HOD <- sapply(fdata, function(x){
    as.numeric(gsub("(↓|↑|<|>)", "", x[grep(">", x)][1]))
  })

  # replace  lod to -Inf Inf
  fdata_na <- apply(fdata, 2, function(x){
    case_when(
      grepl("(<|↓)", x) ~ -Inf,
      grepl("(>|↑)", x) ~ Inf,
      TRUE ~ as.numeric(x)
    )
  })

  # get missing lod pct
  rowdata$oor_pct <- sapply(data.frame(is.infinite(fdata_na)), function(x){
    round(sum(x) / length(x) * 100, 2)
  })


  # replace below lod to facevalue
  fdata <- apply(fdata, 2, function(x){
    as.numeric(gsub("(<|↓|>|↑)", "", x))
  })

  colnames(fdata_na) <- colnames(fdata) <-rowdata$Analyt
  rownames(fdata_na) <- rownames(fdata) <- coldata$unique_id

  #---flour intensity raw
  #sheet.list <- make.names(readxl::excel_sheets(f))
  sheet.list <- rowdata$Analyt

  rdata <- lapply(rowdata$Analyt, function(x){
    #temp <- read.xlsx(f, sheetIndex = which(sheet.list == x), header = T, startRow = analyt_start_row)
    if(grepl("xls$", f)){
      temp <- readxl::read_xls(f, sheet = (2 + which(sheet.list == x)), col_names = T, skip = (analyt_start_row - 1))
      temp <- temp[1: (grep("Note", temp$Location)-2), ]
    }
    if(grepl("xlsx$", f)){
      temp <- read.xlsx(f, sheetIndex = (2 + which(sheet.list == x)), header = T, startRow = analyt_start_row)
      temp <- temp[1: (grep("Note", temp$Location)-1), ]
    }


    temp <- temp%>%
      dplyr::select(Location, Sample, MFI, CV)%>%
      dplyr::filter(!is.na(Sample))%>%
      magrittr::set_colnames(value = c("Location", "Sample",
                             paste0(x, "_MFI"),
                             paste0(x, "_CV")))
    temp
  })

  rstd <- lapply(rowdata$Analyt, function(x){
    #temp <- read.xlsx(f, sheetIndex = which(sheet.list == x), header = T, startRow = analyt_start_row)
    if(grepl("xls$", f)){
      temp_std <- readxl::read_xls(f, sheet = (2 + which(sheet.list == x)), col_names = T, skip = 6)
      temp_std <- temp_std[1: (grep("Note", unlist(temp_std[ , 1]))[1]-2), ]
    }
    if(grepl("xlsx$", f)){
      temp_std <- read.xlsx(f, sheetIndex = (2 + which(sheet.list == x)), header = T, startRow = 7)
      temp_std <- temp_std[1: (grep("Note", unlist(temp_std[ , 1]))[1]-1), ]
    }

    idx_l <- which(!is.na(temp_std$MFI))[2]
    idx_h <- which(!is.na(temp_std$MFI))
    idx_h <- idx_h[length(idx_h)]
    c(temp_std[idx_l, "MFI"], temp_std[idx_l, 11], # lowest std gradient's mfi and conc.
      temp_std[idx_h, "MFI"], temp_std[idx_h, 11], # highest std gradient's mfi and conc.
      temp_std[temp_std$Location == "1E1", "MFI"], temp_std[temp_std$Location == "1E1", 11], temp_std[temp_std$Location == "1E1", 13],
      temp_std[temp_std$Location == "1G1", "MFI"], temp_std[temp_std$Location == "1G1", 11], temp_std[temp_std$Location == "1G1", 13])%>%
      as.numeric()# std3 conc.
  })%>%
    do.call(what = rbind)%>%
    magrittr::set_colnames(value = c("low_MFI", "low_conc.", "high_MFI", "high_conc.", "std3_MFI", "std3_conc.", "std3_cv",
                           "std4_MFI", "std4_conc.", "std4_cv"))

  rowdata <- rowdata%>%
    cbind(data.frame(rstd))

  merge.rdata <- rdata[[1]]
  #merge.rdata$Sample[which(merge.rdata$Sample == "MM-0535-96-H_")] <- "MM-0535-96-H_2"
  if(length(rdata) > 1){
    for (i in 2 : length(rdata)) {
      merge.rdata <- merge(merge.rdata, rdata[[i]])
    }
  }
  rownames(merge.rdata) <- paste(merge.rdata$Location, merge.rdata$Sample, sep = "_")

  rdata <- data.frame(merge.rdata[match(rownames(fdata), rownames(merge.rdata)) , grep("MFI$", colnames(merge.rdata))])

  cvs <- data.frame(merge.rdata[match(rownames(fdata), rownames(merge.rdata)) , grep("CV$", colnames(merge.rdata))])
  cvs <- sapply(cvs, function(x){
    as.numeric(gsub("%", "", x))/100
  })
  colnames(rdata) <- colnames(cvs) <- colnames(fdata)
  rownames(rdata) <- rownames(cvs) <- rownames(fdata)

  re <- SummarizedExperiment(colData = coldata,
                             rowData = rowdata,
                             assays = list(data_default = t(fdata_na),
                                           data_imputed = t(fdata),
                                           mfi_default = t(rdata),
                                           cv = t(cvs)),
                             metadata = list("file_name" = toString(f)))
  return(re)
}

