#'count_qc_plot
#'@description  Heatmap plot count location to visualize the
#'             location with beads count < 35. red color indicate
#'             low beads count
#'
#'
#'@param x component object created by read_lmx_csv()
#'@param pre_fix pre_fix for the csv to be saved.
#'@md
#'
count_qc_plot <- function(f_list, pre_fix){
  for (x in 1 : length(f_list)) {
    temp <- f_list[[x]]$unit[[grep("DataType:,Count", names(f_list[[x]]$unit))]]

    temp <- temp[ , 1 : (which(colnames(temp) == "Total Events") - 1)]

    mat <- temp[ ,-ncol(temp)]%>%
      select(-Location)%>%
      group_by(Sample)%>%
      mutate_all(as.numeric)%>%
      mutate_all(.funs = function(y) {
        case_when(
          y >= 50 ~ y * 0.01,
          (y < 50) & (y >= 35) ~ 5,
          y < 35 ~ 10
        )
      })

    col_dat <- data.frame(temp[, 1:2])
    mat <- t(mat[, -1])

    colnames(mat) <- rownames(col_dat) <- paste(col_dat$Location, col_dat$Sample)

    f_name <- paste0(pre_fix, gsub("(^.*/|.csv)", "", names(f_list)[x]), ".png")
    png(filename = f_name, width = 600 + ncol(mat) * 10, height = 600 + nrow(mat) * 10)
    pheatmap::pheatmap(mat, annotation_col = col_dat, cluster_cols = F, cluster_rows = F, scale = "row")
    dev.off()
  }
}
