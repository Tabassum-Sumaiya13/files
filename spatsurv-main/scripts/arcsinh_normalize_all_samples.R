library(emconnect)
library(SpatialMap)

con <- initialize_DB_connection()
con <- emconnect()

output_dir = 'data'
dir.create(output_dir, showWarnings = FALSE)

setwd(output_dir)

# Getting region ids and list of markers from the UPMC HNC1 dataset
acquisition_ids <- list_acquisition_ids('UPMC')$ACQUISITION_ID
markers <- features(spatialmap_from_db(con, schema="CORE_DATA", acquisition_ids=acquisition_ids[1], study_id=240))

# create dataframe with all cells x markers. columns include marker values, cell id, and sample
all_normalized_data <- data.frame(matrix(ncol=length(markers)+2, nrow=0))
colnames(all_normalized_data) <- c(markers, 'sample', 'cell_id')
idx=1
for (i in 1:length(acquisition_ids)) {
  sm <- spatialmap_from_db(con, schema="CORE_DATA", acquisition_ids=acquisition_ids[i], study_id=240)
  sm <- sm %>% Normalize('asinh')
  num_cells <- dim(sm[1]@regions[[1]]@NormalizedData)[2]
  if (num_cells > 0) {
    all_normalized_data[idx:(idx+num_cells-1),1:length(markers)] <- t(sm[1]@regions[[1]]@NormalizedData)
    all_normalized_data[idx:(idx+num_cells-1),(length(markers)+1):(length(markers)+2)] <- sm[1]@regions[[1]]@cellMetadata
    idx <- idx + num_cells
  }
}

# Save arcsinh-normalized data to a local RDS and csv file
saveRDS(all_normalized_data, 'all_arcsinh_norm_data.RDS')
write.csv(all_normalized_data,'all_arcsinh_norm_data.csv')

all_normalized_data <- readRDS('all_arcsinh_norm_data.RDS')

# create dataframe with all cells x coordinates. columns include coordinates, cell id, and sample
all_coords_data <- data.frame(matrix(ncol=4, nrow=0))
colnames(all_coords_data) <- c('x', 'y', 'sample', 'cell_id')
idx=1
for (i in 1:length(acquisition_ids)) {
  print(i)
  sm <- spatialmap_from_db(con, schema="CORE_DATA", acquisition_ids=acquisition_ids[i], study_id=240)
  sm <- sm %>% Normalize('asinh')
  num_cells <- dim(sm[1]@regions[[1]]@NormalizedData)[2]
  if (num_cells > 0) {
    all_coords_data[idx:(idx+num_cells-1),1:2] <- sm[1]@regions[[1]]@coordinates[,6:7]
    all_coords_data[idx:(idx+num_cells-1),3:4] <- sm[1]@regions[[1]]@cellMetadata
    idx <- idx + num_cells
  }
}

saveRDS(all_coords_data, 'all_sm_xy_data.RDS')
write.csv(all_coords_data,'all_sm_xy_data.csv')

all_coords_data <- readRDS('all_sm_xy_data.RDS')

