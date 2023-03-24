
# Raster Processing -------------------------------------------------------
library(terra)
library(tictoc)
library(exactextractr)
library(tidyverse)
library(sf)
library(feather)


#Specify Drive Path
drive_path <- "//worldpop.files.soton.ac.uk/worldpop/Projects/WP517763_GRID3/"
shapefile_path <-  paste0(drive_path, "Working/GHA/Ortis/Shapefile/")
raster_path <- paste0(drive_path, "Working/GHA/Ortis/Final_covariates/")
raster_path2 <- paste0(drive_path, "Working/GHA/Ortis/Other_covariates/")
output_path <- paste0(drive_path, "Working/GHA/Ortis/Output/")



#Read dataset
GHA_Data_shp <- st_read(paste0(shapefile_path, "Ghana_District.shp"))


# Extract Covariates ----------------------------------------------

#import all covariates and stack them

raster_list <-list.files(path=raster_path, pattern=".tif$", all.files=TRUE, full.names=FALSE)
raster_list

#Stack all covariates
raster_covariates <- rast(paste0(raster_path, c(raster_list)))

#Project CMR_Data to same spatial reference as raster and select variables
GHA_Data_shp <- st_transform(GHA_Data_shp, crs = st_crs(raster_covariates))


#Extract rasters using their mean values
tic()

raster_extract <- exactextractr::exact_extract(raster_covariates, GHA_Data_shp, fun = 'mean')

toc()

#Extract variable names
var_names <- names(raster_extract)

#Change names
colnames(raster_extract) <- c(paste0('x', 1:24))

#Extract names of raster
var_names2<- names(raster_extract)

#cbind names
var_names <- cbind(var_names, var_names2) %>% 
  as_tibble()

#Export names
write.csv(var_names, paste0(output_path, "var_names.csv"))

#Cbind raster_extract to data

GHA_Data_shp <- GHA_Data_shp %>% 
  cbind(raster_extract)


#write_sf(GHA_Data_shp, paste0(output_path,"GHA_Data.shp"))

#convert shapefile to dataframe
GHA_Data_df<- as.data.frame(GHA_Data_shp)

GHA_Data_df <- GHA_Data_df %>% 
  dplyr::select(-geometry)

write.csv(GHA_Data_df, paste0(output_path, "GHA_Data_df.csv"))

##############################################################################################


# Extract dataset for prediction ------------------------------------------

masterGrid<- rast(paste0(raster_path2, "GHA_Mastergrid.tif"))
Grid_ID <-  rast(paste0(raster_path2, "Grid_ID.tif"))

#stack rasters

raster_stack <- c(raster_covariates, masterGrid, Grid_ID)

#get raster values
raster_df <- terra::values(raster_stack, dataframe = T) %>% 
  rename(Dist_ID = GHA_Mastergrid )

#filter only pixels with building count
settled_df <- raster_df %>% 
  filter(Dist_ID >0, GHA_building_count >0) 

#Check for duplicate
any(duplicated(settled_df$Grid_ID))

##Change names of covariate
colnames(settled_df) <- c(paste0('x', 1:24), "Dist_ID", "Grid_ID")

#Replace NA values with 0
#settled_df <- settled_df %>% 
 # mutate_all(~replace_na(.x, 0))

tic()

write_feather(settled_df, paste0(output_path, "settled_df.feather"))
#write.csv(settled_df, paste0(output_path, "settled_df.csv"))
toc()












