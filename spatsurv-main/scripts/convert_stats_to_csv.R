library(spatstat)

#' Create dataframe row names for K functions
#'
#' @param sample_ids list of sample/region acquisition ids
#' @param extension string value indicating the end substring of the RDS files
#' @return Returns a string vector of row names
create_rownames <- function(sample_ids, extension, markers) {
  row_names <- vector(mode='list', length=length(sample_ids)*length(markers))
  for (i in 1:length(sample_ids)) {
    for (j in 1:length(markers)) {
      row_names[(i-1)*length(sample_ids) + j] <- sprintf('%s_%s_%s', sample_ids[i], markers[j], extension)
    }
  }
  return(unlist(row_names))
}

#' Read K function files into a dataframe
#'
#' @param file_dir string of path to directory with the RDS files
#' @param title string indicating end substring of RDS files
#' @param sample_ids list of sample/region acquisition ids
#' @param markers list of biomarkers
#' @param iso optional boolean indicating if RDS object requires $iso to get the values
#' @return dataframe with columns: r1-->r401 for the K fn values for each sample
read_files_to_df <- function(file_dir, title, sample_ids, markers, iso=FALSE) {
  row_names <- create_rownames(sample_ids, title, markers)
  fn_data <- data.frame(matrix(ncol=401, nrow=length(sample_ids)*length(markers)), row.names=row_names) 
  for (sample_id in sample_ids){
    for (marker in markers) {
      file_name <- sprintf('%s_%s_%s', sample_id, marker, title)
      fn_vals <- readRDS(sprintf('%s/%s.RDS', file_dir, file_name))
      if (iso) {
        fn_data[file_name, ] <- fn_vals$iso
      } else {
        fn_data[file_name, ] <- fn_vals
      }
    }
  }
  saveRDS(fn_data, sprintf('%s_df.RDS', title))
  return(fn_data)
}

# run get_study_metadata.R first
setwd('./')

norm_marker_data <- readRDS('data/all_arcsinh_norm_data.RDS')
markers <- colnames(norm_marker_data)[1:40]

Lmark_norm_df <- read_files_to_df('k_fns_norm_by_uw', 
                                  'Lmark_normalized', 
                                  qc_acq_ids_labeled, 
                                  markers, 
                                  iso=FALSE)
write.csv(Lmark_norm_df, 'data/k_fns_norm_by_uw_qc_labeled.csv', quote=FALSE)

