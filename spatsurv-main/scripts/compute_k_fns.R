library(spatstat)
library(parallel)

compute_weighted_k <- function(sample_name, coords_df, marks_df, output_dir, markers) {
  dir.create(output_dir, showWarnings = FALSE)
  sub_coords_df <- subset(coords_df, sample==sample_name)
  sub_marks_df <- subset(marks_df, sample==sample_name)
  x <- sub_coords_df$x
  y <- sub_coords_df$y
  for (marker in markers) {
    if (file.exists(sprintf('%s/%s_%s_Lmark.RDS', output_dir, sample_name, marker))) {
      next
    }
    pp <- ppp(x,y, c(0, max(x)+1), c(0, max(y)+1), marks=sub_marks_df[[marker]])
    Lm <- Kmark(pp, correction='isotropic', returnL=TRUE, r=seq(0,400,length.out=401))
    saveRDS(Lm, file=sprintf('%s/%s_%s_Lmark.RDS', output_dir, sample_name, marker))
  }
  return()
}
compute_unweighted_k <- function(sample_name, coords_df, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  if (file.exists(sprintf('%s/%s_L_unweighted.RDS', output_dir, sample_name))) {
    return()
  }
  sub_coords_df <- subset(coords_df, sample==sample_name)
  x <- sub_coords_df$x
  y <- sub_coords_df$y
  pp <- ppp(x,y, c(0, max(x)+1), c(0, max(y)+1))
  L <- Lest(pp, correction='isotropic', returnL=TRUE, r=seq(0,400,length.out=401))
  saveRDS(L, file=sprintf('%s/%s_L_unweighted.RDS', output_dir, sample_name))
  return()
}
compute_norm_by_uw_k <- function(sample_name, output_dir, markers, weighted_dir, unweighted_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  for (marker in markers) {
    if (file.exists(sprintf('%s/%s_%s_Lmark_normalized.RDS', output_dir, sample_name, marker))) {
      next
    }
    weighted_k <- readRDS(sprintf('%s/%s_%s_Lmark.RDS', weighted_dir, sample_name, marker))
    unweighted_k <- readRDS(sprintf('%s/%s_L_unweighted.RDS', unweighted_dir, sample_name))
    normalized_k <- weighted_k$iso - unweighted_k$iso
    
    saveRDS(normalized_k, sprintf('%s/%s_%s_Lmark_normalized.RDS', output_dir, sample_name, marker))
  }
}

setwd('./')
xy_data <- readRDS('data/all_sm_xy_data.RDS')
norm_marker_data <- readRDS('data/all_arcsinh_norm_data.RDS')

sample_list <- unique(xy_data[['sample']])
markers <- colnames(norm_marker_data)[1:40]

print('Computing unweighted k functions')
mclapply(sample_list, compute_unweighted_k, xy_data, 'k_unweight_fns', mc.cores = detectCores())
print('Computing weighted k functions')
mclapply(sample_list, compute_weighted_k, xy_data, norm_marker_data, 'kmark_fns', markers, mc.cores = detectCores())
print('Computing weighted minus unweighted k functions')
mclapply(sample_list, compute_norm_by_uw_k, 'k_fns_norm_by_uw', markers, 'k_unweight_fns', 'k_unweight_fns', mc.cores = detectCores())


