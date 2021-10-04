#' update samples with different dilution factors
#'@rdname update_dilute
#'@description update samples with different dilution factors
#'@param se npx summarizedexpermient object
#'@param sample_pattern string pattern to sample
#'@param fields colname in metadata to filter sample
#'@param dilute.f numeric dilution factor
#'@export
#'@md
#'
#'
update_dilute <- function(se, sample_pattern = "QC", fields = "Sample", dilute.f = 1){
  se@assays@data$data_default[ , grep(sample_pattern, se[[fields]])] <- se@assays@data$data_default[ , grep(sample_pattern, se[[fields]])] * dilute.f
  se
}
