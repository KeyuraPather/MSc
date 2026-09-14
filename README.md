#1. CLIMATE:

# Install and Load Required Libraries
#install.packages(c("terra", "dplyr"))
library(terra)
library(dplyr)

##downloaded and extracted relevant layers (precip and temp) from CRU

# Step 1: Load CRU Data
precip <- rast("precip")
temp <- rast("temp")

# Subset data for 1940-1975 
time_index <- 12 * (1940 - 1901) + 1:(12 * (1975 - 1940 + 1))
precip <- precip[[time_index]]
temp <- temp[[time_index]]

# Generate Year-Month Labels for Time Period
start_date <- as.Date("1940-01-01")
dates <- seq(start_date, by = "month", length.out = length(time_index))
date_labels <- format(dates, "%Y-%m")  # Format as "YYYY-MM"

# Step 2: Load and Process Species Distribution Data
species_data <- read.csv("Final_occurrences_clean.csv")  # Load CSV with longitude & latitude
species_points <- vect(species_data, geom = c("DDE", "DDS"), crs = crs(precip))  # Convert to spatial points

# Step 3: Extract Climate Data for Species Locations
precip_values <- extract(precip, species_points)  # Extract precipitation values
temp_values <- extract(temp, species_points)      # Extract temperature values

# Rename Columns with Year-Month Labels
colnames(precip_values)[-1] <- paste0("precip_", date_labels)  # Exclude ID column
colnames(temp_values)[-1] <- paste0("temp_", date_labels)      # Exclude ID column

# Step 4: Combine Climate Data with Species Data
species_clipped <- cbind(species_data, precip_values[, -1], temp_values[, -1])  # Remove ID columns

summary(species_clipped)

#saw was not filtered properly by date so filtered again
species_clipped <- species_clipped[species_clipped$Coll_Year >= 1940 & species_clipped$Coll_Year <= 1975, ]

species_clipped <- na.omit(species_clipped)


# Step 5: Save Results to CSV
write.csv(species_clipped, "Final/New_cru.csv", row.names = FALSE)



##AAD CALCULATIONS
# Load libraries
#install.packages("dplyr")
library(dplyr)

# Step 1: Load Data
data <- read.csv("Final/New_cru.csv")  # Replace with your file path

# Step 2: Extract Monthly Data
# Select columns with monthly precipitation and temperature
monthly_precip <- data %>% select(starts_with("precip_"))
monthly_temp <- data %>% select(starts_with("temp_"))

# Step 3: Calculate AAD for Precipitation
aad_precip <- data %>%
  mutate(
    mean_precip = rowMeans(monthly_precip, na.rm = TRUE),  # Mean monthly precip per species
    n_precip = ncol(monthly_precip),  # Total number of months (data points)
    aad_precip = rowSums(abs(monthly_precip - mean_precip), na.rm = TRUE)  # Sum of absolute deviations
  ) %>%
  mutate(aad_precip = aad_precip / n_precip) %>%  # Divide by number of months to get average
  group_by(X.Species) %>%  # Ensure one row per species
  summarise(aad_precip = unique(aad_precip)[1], .groups = "drop")  # Aggregate to one row per species

# Step 4: Calculate AAD for Temperature
aad_temp <- data %>%
  mutate(
    mean_temp = rowMeans(monthly_temp, na.rm = TRUE),  # Mean monthly temp per species
    n_temp = ncol(monthly_temp),  # Total number of months (data points)
    aad_temp = rowSums(abs(monthly_temp - mean_temp), na.rm = TRUE)  # Sum of absolute deviations
  ) %>%
  mutate(aad_temp = aad_temp / n_temp) %>%  # Divide by number of months to get average
  group_by(X.Species) %>%  # Ensure one row per species
  summarise(aad_temp = unique(aad_temp)[1], .groups = "drop")  # Aggregate to one row per species

# Step 5: Combine Results
aad_results <- left_join(aad_precip, aad_temp, by = "X.Species")

# Step 6: Save Results to CSV
write.csv(aad_results, "Final/species_aad.csv", row.names = FALSE)


#finding the QUANTILES
# Load necessary library
#install.packages("dplyr")
library(dplyr)

# Step 1: Load the AAD results CSV file
aad_results <- read.csv("Final/species_aad.csv")  # Make sure this is the correct path

# Step 2: Calculate the quantiles for precipitation and temperature separately
quantiles_precip <- quantile(aad_results$aad_precip, probs = c(0.25, 0.75), na.rm = TRUE)
quantiles_temp <- quantile(aad_results$aad_temp, probs = c(0.25, 0.75), na.rm = TRUE)

# Step 3: Categorize species based on AAD values for Precipitation
aad_results <- aad_results %>%
  mutate(
    category_precip = case_when(
      aad_precip <= quantiles_precip[1] ~ "Lowest 25%",  # AAD for precipitation <= 25th percentile
      aad_precip > quantiles_precip[1] & aad_precip <= quantiles_precip[2] ~ "25-75%",  # AAD for precipitation between 25th and 75th percentiles
      aad_precip > quantiles_precip[2] ~ "Highest 75%"  # AAD for precipitation >= 75th percentile
    ),
    
    category_temp = case_when(
      aad_temp <= quantiles_temp[1] ~ "Lowest 25%",  # AAD for temperature <= 25th percentile
      aad_temp > quantiles_temp[1] & aad_temp <= quantiles_temp[2] ~ "25-75%",  # AAD for temperature between 25th and 75th percentiles
      aad_temp > quantiles_temp[2] ~ "Highest 75%"  # AAD for temperature >= 75th percentile
    )
  )

# Step 4: Save the categorized results to a single CSV file
write.csv(aad_results, "Final/categorized_aad.csv", row.names = FALSE)



# 2. EXPOSURE
#create climate stacks
library(terra)

# Define paths to your data, changed folder name for each model
current_path <- "MIROC6/current/" # Folder containing historical BIO layers
future_126_path <- "MIROC6/rcp_126.tif"  # Multiband file for RCP 1.26
future_585_path <- "MIROC6/rcp_585.tif"  # Multiband file for RCP 5.85

# Define variables to include (BIO1, BIO5, BIO10, BIO11, BIO12, BIO15)
variables <- c(1, 5, 10, 11, 12, 15)

# Load and stack historical (current) data using terra
current_stack <- rast(paste0(current_path, "bio", variables, ".tif"))

# Load future scenarios and extract selected variables using terra
future_126_stack <- rast(future_126_path)[[variables]]
future_585_stack <- rast(future_585_path)[[variables]]

# Define the target CRS (WGS84 as an example)
target_crs <- "+proj=longlat +datum=WGS84 +no_defs"

# Project each stack to the target CRS
current_stack <- project(current_stack, target_crs)
future_126_stack <- project(future_126_stack, target_crs)
future_585_stack <- project(future_585_stack, target_crs)

# Optional: Align resolution and extent (if necessary)
future_126_stack <- resample(future_126_stack, current_stack)
future_585_stack <- resample(future_585_stack, current_stack)

# Check the CRS of each stack to ensure they match
print(crs(current_stack))  # Should be WGS84
print(crs(future_126_stack))  # Should match current stack's CRS
print(crs(future_585_stack))  # Should match current stack's CRS

# Load the terra package
library(terra)

# Assuming you have a RasterLayer object (e.g., 'current_maxent_pred')
writeRaster(current_stack, "MIROC6/current.tif", overwrite=TRUE)
writeRaster(future_126_stack, "MIROC6/future126.tif", overwrite=TRUE)
writeRaster(future_585_stack, "MIROC6/future585.tif", overwrite=TRUE)

# Load necessary libraries
library(raster)
library(bioclim)
library(sp)
library(dismo)
library(pROC)

# Define Albers Equal Area Projection CRS for Southern Africa
aea_crs <- "+proj=aea +lat_1=-30 +lat_2=-20 +lon_0=23 +datum=WGS84 +units=m +no_defs"

# Load bioclim and CMIP data for current, RCP 126, and RCP 585
bioclim_current <- stack("MIROC6/current.tif")
bioclim_rcp126 <- stack("MIROC6/future126.tif")
bioclim_rcp585 <- stack("MIROC6/future585.tif")

# Define the extent for Southern Africa (Longitude: 10°E to 32°E, Latitude: -35°S to -10°S)
southern_africa_extent <- extent(10, 32, -35, -10)

# Crop the original bioclimatic raster to this extent
bioclim_current <- crop(bioclim_current, southern_africa_extent)
bioclim_rcp126 <- crop(bioclim_rcp126, southern_africa_extent)
bioclim_rcp585 <- crop(bioclim_rcp585, southern_africa_extent)

# Reproject bioclimatic data to Albers Equal Area (AEA)
bioclim_current_aea <- projectRaster(bioclim_current, crs = aea_crs, method = "bilinear")
bioclim_rcp126_aea <- projectRaster(bioclim_rcp126, crs = aea_crs, method = "bilinear")
bioclim_rcp585_aea <- projectRaster(bioclim_rcp585, crs = aea_crs, method = "bilinear")

# Define the desired resolution in meters (e.g., 1000 meters = 1 km)
target_resolution <- 1000

# Create a new raster with the desired resolution
new_raster <- raster(crs = crs(bioclim_current_aea))  # Match CRS to the original raster
extent(new_raster) <- extent(bioclim_current_aea)  # Match extent to the original raster
res(new_raster) <- target_resolution  # Set the resolution

# Resample the original raster to the new resolution
bioclim_current_aea <- resample(bioclim_current_aea, new_raster, method = "bilinear")
bioclim_rcp126_aea <- resample(bioclim_rcp126_aea , new_raster, method = "bilinear")
bioclim_rcp585_aea <- resample(bioclim_rcp585_aea, new_raster, method = "bilinear")

# Occurrence data (latitude, longitude)
# Assuming 'occurrence_data' is a data.frame with columns: longitude, latitude
#ONLY RUN FROM HERE AFTER THE FIRST TIME WITH NEXT SPECIES LISTS 
occurrence_data <- read.csv(("vach.csv"))
coordinates(occurrence_data) <- ~longitude + latitude
proj4string(occurrence_data) <- CRS("+proj=longlat +datum=WGS84")

# Reproject occurrence data to AEA
occurrence_data_aea <- spTransform(occurrence_data, CRS(aea_crs))

# Train the species models using bioclim package
model_current <- bioclim(bioclim_current_aea, occurrence_data_aea)
model_rcp126 <- bioclim(bioclim_rcp126_aea, occurrence_data_aea)
model_rcp585 <- bioclim(bioclim_rcp585_aea, occurrence_data_aea)

#Get suitability scores from the bioclim model
model_current <- predict(model_current, bioclim_current_aea)
model_rcp126 <- predict(model_rcp126, bioclim_rcp126_aea)
model_rcp585 <- predict(model_rcp585, bioclim_rcp585_aea)

# Thresholding to create binary presence/absence maps for each model
# Use 10th percentile for threshold (or adjust based on your needs)
threshold_current <- quantile(values(model_current), 0.1, na.rm = TRUE)
threshold_rcp126 <- quantile(values(model_rcp126), 0.1, na.rm = TRUE)
threshold_rcp585 <- quantile(values(model_rcp585), 0.1, na.rm = TRUE)

# Create binary presence/absence maps
binary_model_current <- model_current > threshold_current
binary_model_rcp126 <- model_rcp126 > threshold_rcp126
binary_model_rcp585 <- model_rcp585 > threshold_rcp585

# Calculate total area in m² (for each model)
pixel_area_m2 <- res(bioclim_current_aea)[1] * res(bioclim_current_aea)[2]  # pixel area in square meters

#Calculate total area in km² (for each model)
pixel_area <- (pixel_area_m2) / 1e6

pixel_area

total_area_current <- cellStats(binary_model_current, sum) * pixel_area
total_area_rcp126 <- cellStats(binary_model_rcp126, sum) * pixel_area
total_area_rcp585 <- cellStats(binary_model_rcp585, sum) * pixel_area

# Calculate range overlap (in m²) between current and future scenarios
overlap_rcp126 <- binary_model_current & binary_model_rcp126
overlap_rcp585 <- binary_model_current & binary_model_rcp585

overlap_area_rcp126 <- cellStats(overlap_rcp126, sum) * pixel_area
overlap_area_rcp585<- cellStats(overlap_rcp585, sum) * pixel_area 


# Calculate total range change (in m²) between current and future scenarios
range_change_rcp126<- total_area_rcp126 - total_area_current 
range_change_rcp585 <- total_area_rcp585 - total_area_current 

#AUC
evaluation_data <- data.frame(
  presence = c(1, 0, 1),  #(1 = presence, 0 = absence)
  suitability = values(model_current))  # Use the suitability values for the model predictions

# Calculate AUC
auc_result <- roc(evaluation_data$suitability, evaluation_data$presence)

cat("Total suitable area (current):", total_area_current)
cat("Total suitable area (RCP 126):", total_area_rcp126)
cat("Total suitable area (RCP 585):", total_area_rcp585)
cat("Overlap area (current vs RCP 126):", overlap_area_rcp126)
cat("Overlap area (current vs RCP 585):", overlap_area_rcp585)
cat("Total range change (current vs RCP 126):", range_change_rcp126)
cat("Total range change (current vs RCP 585):", range_change_rcp585)
cat("AUC:", auc_result$auc)


# 3.  VEGETATION TYPES
# Load required libraries
#install.packages(c("sf", "dplyr", "terra")) # Install packages if not already installed
library(sf)
library(dplyr)
library(terra)

# Step 1: Load species occurrence data
species_data <- read.csv("Final_occurrences_clean.csv")  # Replace with your CSV file
species_sf <- st_as_sf(
  species_data,
  coords = c("DDE", "DDS"),  # Replace with your coordinate column names
  crs = 4326                           # WGS 84 coordinate system
)

##loading vegetation type shapefiles
sa<-st_read("vegtypes.shp")
na<- st_read("Namibia_vegtypes.shp")

#SOUTH AFRICA
# Step 2: Load vegetation map
vegetation_map <-sa  # Replace with your shapefile path
st_crs(vegetation_map)
head(st_coordinates(vegetation_map))


# Step 3: Ensure CRS alignment
species_sf <- st_transform(species_sf, st_crs(vegetation_map))

# Step 4: Spatial join to link species points with vegetation polygons
joined_data <- st_join(species_sf, vegetation_map)

# Step 5: Count unique vegetation types for each species
veg_types_per_species <- joined_data %>%
  group_by(X.Species) %>%          # Replace 'species_name' with your species identifier column
  summarise(num_veg_types = n_distinct(T_Name))  # Replace 'vegetation_type' with the appropriate column name

# Step 6: View and export results
print(veg_types_per_species)  # Display the result
write.csv(veg_types_per_species, "veg_types_per_species.csv")  # Export to CSV



#Namibia
# Step 2: Load vegetation map
vegetation_map <- na # Replace with your shapeless path

# Check for invalid geometries
invalid_geom <- st_is_valid(vegetation_map)

# Print summary
summary(invalid_geom)

# Fix invalid geometries
vegetation_map <- st_make_valid(vegetation_map)

# Verify again
summary(st_is_valid(vegetation_map))

# Step 3: Ensure CRS alignment
species_sf <- st_transform(species_sf, st_crs(vegetation_map))

# Step 4: Spatial join to link species points with vegetation polygons
joined_data <- st_join(species_sf, vegetation_map)

# Step 5: Count unique vegetation types for each species
veg_types_per_species <- joined_data %>%
  group_by(X.Species) %>%          # Replace 'species_name' with your species identifier column
  summarise(num_veg_types = n_distinct(VEG_TYPE))  # Replace 'vegetation_type' with the appropriate column name

# Step 6: View and export results
print(veg_types_per_species)  # Display the result
write.csv(veg_types_per_species, "Namibia_veg_types_per_species.csv")  # Export to CSV



# 4. VEGETATION TYPES STATS

library(sf)
library(dplyr)
library(terra)


###checking the resolution of each file and chnaging where needed

sa<-st_read("vegtypes.shp")
na<- st_read("Namibia_vegtypes.shp")

#1) crs
st_crs(sa)
st_crs(na)

na <- st_transform(na, crs = st_crs(sa)) #now both have same crs

#2) precision
head(st_coordinates(sa))
head(st_coordinates(na)) #same decimals so equal precision for each

summary(st_coordinates(st_geometry(sa)))
summary(st_coordinates(st_geometry(na))) #same decimals so equal precision for each

#3) density

# Calculate total area (in km2)
area_sa <- sum(st_area(sa)) / 10^6
area_na <- sum(st_area(na)) / 10^6 

# Calculate number of polygons per 1000 km²
density_sa <- nrow(sa) / (area_sa / 1000)
density_na <- nrow(na) / (area_na / 1000)

# Print results
cat("Polygon density per 1000 km²:\n")
cat("South Africa:", round(density_sa, 2), "\n")
cat("Namibia:", round(density_na, 2), "\n")  #south africa has higher density per 1000km so need to simplify it - recomended if doing point in polygon

# Manually simplify - try with 0.01 first
sa_simp <- st_simplify(sa, dTolerance =  0.01, preserveTopology = TRUE)

# Check new density
area_sa <- sum(st_area(sa_simp)) / 10^6
density <- nrow(sa) / (area_sa / 1000)

cat("Density:", density)

# 5. HEAT MAPS 
#ADDITIVE
###BEST CASE
# Load necessary libraries
library(sf)
library(raster)
library(ggplot2)
library(viridis)

#make sure shapefile is correct, I cut mine to make a new shape file (only the first time) - downloaded from The South African Protected Area Database (SAPAD)
# Load the Richtersveld shapefile(only run this first time, dont need otherwise)
#PA <- st_read("Basemap.shp")  # Adjust path to the Richtersveld shapefile

#richtersveld <- PA %>%
#filter(grepl("Richtersveld National Park", CUR_NME, ignore.case = TRUE))  # Adjust "NAME" to your column name
# Step 3: Ensure the CRS is EPSG:4326 (WGS 84)
#richtersveld <- st_transform(richtersveld, crs = 4326)  # Ensure it's in WGS 84

# Step 4: Create a bounding box (for example, for a region of interest)
# Define the coordinates of the bounding box (you can adjust these as needed)
#bbox <- st_sfc(st_polygon(list(matrix(c(16.6, -28.6, 17.8, -28.6,                       17.8, -27.8, 16.6, -27.8, 16.6, -28.6), ncol = 2, byrow = TRUE))))

# Step 5: Assign the same CRS as the Richtersveld shapefile to the bounding box
#st_crs(bbox) <- st_crs(richtersveld)  # Assign CRS to the bbox to match richtersveld's CRS

# Step 6: Clip the Richtersveld shapefile with the bounding box
#richtersveld <- st_intersection(richtersveld, bbox)


# Step 7: Plot the clipped Richtersveld
#plot(richtersveld)

# Save the clipped shapefile to your desired path
#st_write(richtersveld, "richtersveld.shp")


# Load the species occurrence data (assuming it has 'longitude', 'latitude' and 'species_name')
species_data <- read.csv("Final_occurrences_clean.csv")  # Adjust path to your species data CSV
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Load the Richtersveld shapefile(only run this first time, dont need otherwise)
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile


# Ensure the CRS of species data and Richtersveld match
species_sf <- st_transform(species_sf, crs = st_crs(richtersveld))

# Define the species you want to include (example: "Species_A", "Species_B", "Species_C")
selected_species <- c("Albuca etesiogaripensis U.Müll.-Doblies", "Albuca longipes Baker", "Anacampseros albissima Marloth", "Astridia citrina (L.Bolus) L.Bolus", "Bromus pectinatus Thunb.", "Codon schenckii Schinz","Crotalaria humilis Eckl. & Zeyh.", "Cyrtanthus herrei (F.M.Leight.) R.A.Dyer",  "Codon schenckii Schinz", "Crotalaria humilis Eckl. & Zeyh.", "Cyrtanthus herrei (F.M.Leight.) R.A.Dyer", "Dimorphotheca pinnata (Thunb.) Harv.", "Drosanthemum salicola L.Bolus", "Euphorbia rhombifolia Boiss.", "Helichrysum leontonyx DC.", "Isolepis hemiuncialis (C.B.Clarke) J.Raynal", "Mesembryanthemum pellitum Friedrich", "Ozoroa crassinervia (Engl.) R.Fern. & A.Fern.", "Pelargonium antidysentericum (Eckl. & Zeyh.) Kostel. ", "Pelargonium echinatum Curtis",  "Ruschia glauca L.Bolus",  "Stipagrostis geminifolia Nees", "Trachyandra aridimontana J.C.Manning", "Ursinia nana DC.") # Modify with the species of your choice

# Filter species data to include only the selected species
species_filtered <- species_sf[species_sf$X.Species %in% selected_species, ]

# Clip the filtered species data to the Richtersveld boundary
species_clipped <- st_intersection(species_filtered, richtersveld)

# Extract coordinates from species_clipped
coordinates <- st_coordinates(species_clipped)

# Convert the sf object to a data frame, including the coordinates
species_df <- as.data.frame(species_clipped)

# Add the coordinates (longitude and latitude) to the data frame
species_df$DDE <- coordinates[, 1]
species_df$DDS <- coordinates[, 2]

species_df

# Save the data frame as a CSV
write.csv(species_df, "species_clipped.csv", row.names = FALSE)

library(sf)
library(raster)
library(ggplot2)

# Step 1: Read the CSV file containing species data with coordinates
species_data <- read.csv("species_clipped.csv", row.names = NULL)
str(species_data)

# Step 2: Load the Richtersveld boundary shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Step 3: Convert the CSV to an sf object (species data)
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Step 4: Reproject Richtersveld to UTM Zone 34S (EPSG:32734) for meter-based resolution
richtersveld <- st_transform(richtersveld, crs = 32734)
plot(richtersveld)

# Step 5: Reproject species data to the same CRS as the Richtersveld boundary
species_sf <- st_transform(species_sf, crs = 32734)

# Step 6: Clip the species data to the Richtersveld boundary (keeping points within Richtersveld)
species_clipped <- st_intersection(species_sf, richtersveld)

# Step 7: Create a raster grid for the heatmap
species_extent <- st_bbox(richtersveld)  # Get the extent of Richtersveld

# Specify the resolution in meters (e.g., 100 meters)
raster_grid <- raster(extent(richtersveld), res = 6324.56)  # Adjust resolution as needed

# Step 8: Set the CRS of the raster grid to match the Richtersveld boundary's CRS
crs(raster_grid) <- st_crs(richtersveld)$proj4string

# Step 9: Rasterize the species points, counting occurrences
species_raster <- rasterize(st_coordinates(species_clipped), raster_grid, field = 1, fun = "count")

library(ggplot2)
library(ggspatial)  # For adding scale bar

# Step 10: Convert raster to a data frame for ggplot
raster_df <- as.data.frame(species_raster, xy = TRUE, na.rm = TRUE)
raster_df

# Step 11: Plot the heatmap using ggplot2 without coordinates
ABC<-ggplot(data = raster_df) +
  geom_tile(aes(x = x, y = y, fill = layer)) +  # Plot the heatmap
  scale_fill_viridis_c(breaks = seq(0, max(raster_df$layer), by = 2), name = NULL) +  # Set legend breaks to increments of 2
  coord_fixed() +  # Keep aspect ratio 1:1
  ggtitle("Additive Best-case Scenario") +  # Add title
  geom_sf(data = richtersveld, fill = NA, color = "black", size = 1) +  # Add Richtersveld boundary
  theme_minimal() +  # Use a clean theme
  theme(legend.position = "bottom",  # Position legend at the bottom
        axis.text = element_blank(),  # Remove axis text (coordinates)
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank()) +  # Remove axis ticks
  annotation_scale(location = "bl", width_hint = 0.1)+
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"))  # Add north arrow

ABC


###WORST CASE
# Load the species occurrence data (assuming it has 'longitude', 'latitude' and 'species_name')
species_data <- read.csv("Final_occurrences_clean.csv")  # Adjust path to your species data CSV
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Load the Richtersveld shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Ensure the CRS of species data and Richtersveld match
species_sf <- st_transform(species_sf, crs = st_crs(richtersveld))

# Define the species you want to include (example: "Species_A", "Species_B", "Species_C") - change for bc and wc
selected_species <- c("Albuca etesiogaripensis U.Müll.-Doblies", "Aloe meyeri Van Jaarsv." , "Astridia citrina (L.Bolus) L.Bolus",  "Ceropegia herrei  (A.C. White & B. Sloane)  Bruyns", "Cheilanthes namaquensis (Baker) Schelpe & N.C.Anthony", "Codon schenckii Schinz", "Crassula subacaulis Schönland & Baker f.", "Crotalaria humilis Eckl. & Zeyh.", "Cucumis rigidus E.Mey. ex Sond.", "Diclis petiolaris Benth.", "Gymnosporia gariepensis M. Jordaan", "Hemarthria altissima (Poir.) Stapf & C.E.Hubb.", "Isolepis hemiuncialis (C.B.Clarke) J.Raynal", "Manulea robusta Pilg.", "Mesembryanthemum pellitum Friedrich",  "Ozoroa crassinervia (Engl.) R.Fern. & A.Fern.", "Pelargonium antidysentericum (Eckl. & Zeyh.) Kostel. ", "Pelargonium echinatum Curtis", "Ruschia glauca L.Bolus", "Schwantesia herrei L.Bolus", "Stipagrostis geminifolia Nees", "Trachyandra aridimontana J.C.Manning", "Tylecodon reticulatus (L.f.) Toelken") # Modify with the species of your choice

# Filter species data to include only the selected species
species_filtered <- species_sf[species_sf$X.Species %in% selected_species, ]

# Clip the filtered species data to the Richtersveld boundary
species_clipped <- st_intersection(species_filtered, richtersveld)

# Extract coordinates from species_clipped
coordinates <- st_coordinates(species_clipped)

# Convert the sf object to a data frame, including the coordinates
species_df <- as.data.frame(species_clipped)

# Add the coordinates (longitude and latitude) to the data frame
species_df$DDE <- coordinates[, 1]
species_df$DDS <- coordinates[, 2]

species_df

# Save the data frame as a CSV
write.csv(species_df, "AWCspecies_clipped.csv", row.names = FALSE)

library(sf)
library(raster)
library(ggplot2)

# Step 1: Read the CSV file containing species data with coordinates
species_data <- read.csv("AWCspecies_clipped.csv", row.names = NULL)
str(species_data)

# Step 2: Load the Richtersveld boundary shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Step 3: Convert the CSV to an sf object (species data)
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Step 4: Reproject Richtersveld to UTM Zone 34S (EPSG:32734) for meter-based resolution
richtersveld <- st_transform(richtersveld, crs = 32734)
plot(richtersveld)

# Step 5: Reproject species data to the same CRS as the Richtersveld boundary
species_sf <- st_transform(species_sf, crs = 32734)

# Step 6: Clip the species data to the Richtersveld boundary (keeping points within Richtersveld)
species_clipped <- st_intersection(species_sf, richtersveld)

# Step 7: Create a raster grid for the heatmap
species_extent <- st_bbox(richtersveld)  # Get the extent of Richtersveld

# Specify the resolution in meters (e.g., 100 meters)
raster_grid <- raster(extent(richtersveld), res = 6324.56)  # Adjust resolution as needed

# Step 8: Set the CRS of the raster grid to match the Richtersveld boundary's CRS
crs(raster_grid) <- st_crs(richtersveld)$proj4string

# Step 9: Rasterize the species points, counting occurrences
species_raster <- rasterize(st_coordinates(species_clipped), raster_grid, field = 1, fun = "count")

library(ggplot2)
library(ggspatial)  # For adding scale bar

# Step 10: Convert raster to a data frame for ggplot
raster_df <- as.data.frame(species_raster, xy = TRUE, na.rm = TRUE)
raster_df

# Step 11: Plot the heatmap using ggplot2 without coordinates
AWC<-ggplot(data = raster_df) +
  geom_tile(aes(x = x, y = y, fill = layer)) +  # Plot the heatmap
  scale_fill_viridis_c(breaks = seq(0, max(raster_df$layer), by = 2), name = NULL) +  # Set legend breaks to increments of 2
  coord_fixed() +  # Keep aspect ratio 1:1
  ggtitle("Additive Worst-case Scenario") +  # Add title
  geom_sf(data = richtersveld, fill = NA, color = "black", size = 1) +  # Add Richtersveld boundary
  theme_minimal() +  # Use a clean theme
  theme(legend.position = "bottom",  # Position legend at the bottom
        axis.text = element_blank(),  # Remove axis text (coordinates)
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank() ) +  # Remove axis ticks
  annotation_scale(location = "bl", width_hint = 0.1)+
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"))  # Add north arrow

AWC

library(patchwork)
ABC+AWC




##ORDINAL
###BEST CASE
# Load the species occurrence data (assuming it has 'longitude', 'latitude' and 'species_name')
species_data <- read.csv("Final_occurrences_clean.csv")  # Adjust path to your species data CSV
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Load the Richtersveld shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Ensure the CRS of species data and Richtersveld match
species_sf <- st_transform(species_sf, crs = st_crs(richtersveld))

# Define the species you want to include (example: "Species_A", "Species_B", "Species_C") - change for bc and wc
selected_species <- c("Adromischus marianiae (Marloth) A.Berger", "Albuca longipes Baker", "Anacampseros albissima Marloth", "Bromus pectinatus Thunb.","Ceropegia articulata (Aiton) Bruyns", "Ceropegia perlata (Dinter) Bruyns", "Crotalaria humilis Eckl. & Zeyh.", "Cucumis rigidus E.Mey. ex Sond.", "Euphorbia gummifera Boiss.", "Helichrysum leontonyx DC.", "Hermannia eenii Baker f.", "Lapeirousia littoralis Baker","Manulea robusta Pilg.","Mesembryanthemum pellitum Friedrich" , "Pelargonium antidysentericum (Eckl. & Zeyh.) Kostel.", "Ursinia nana DC.") # Modify with the species of your choice

# Filter species data to include only the selected species
species_filtered <- species_sf[species_sf$X.Species %in% selected_species, ]

# Clip the filtered species data to the Richtersveld boundary
species_clipped <- st_intersection(species_filtered, richtersveld)

# Extract coordinates from species_clipped
coordinates <- st_coordinates(species_clipped)

# Convert the sf object to a data frame, including the coordinates
species_df <- as.data.frame(species_clipped)

# Add the coordinates (longitude and latitude) to the data frame
species_df$DDE <- coordinates[, 1]
species_df$DDS <- coordinates[, 2]

species_df

# Save the data frame as a CSV
write.csv(species_df, "OBCspecies_clipped.csv", row.names = FALSE)

library(sf)
library(raster)
library(ggplot2)

# Step 1: Read the CSV file containing species data with coordinates
species_data <- read.csv("OBCspecies_clipped.csv", row.names = NULL)
str(species_data)

# Step 2: Load the Richtersveld boundary shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Step 3: Convert the CSV to an sf object (species data)
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Step 4: Reproject Richtersveld to UTM Zone 34S (EPSG:32734) for meter-based resolution
richtersveld <- st_transform(richtersveld, crs = 32734)
plot(richtersveld)

# Step 5: Reproject species data to the same CRS as the Richtersveld boundary
species_sf <- st_transform(species_sf, crs = 32734)

# Step 6: Clip the species data to the Richtersveld boundary (keeping points within Richtersveld)
species_clipped <- st_intersection(species_sf, richtersveld)

# Step 7: Create a raster grid for the heatmap
species_extent <- st_bbox(richtersveld)  # Get the extent of Richtersveld

# Specify the resolution in meters (e.g., 100 meters)
raster_grid <- raster(extent(richtersveld), res = 6324.56)  # Adjust resolution as needed

# Step 8: Set the CRS of the raster grid to match the Richtersveld boundary's CRS
crs(raster_grid) <- st_crs(richtersveld)$proj4string

# Step 9: Rasterize the species points, counting occurrences
species_raster <- rasterize(st_coordinates(species_clipped), raster_grid, field = 1, fun = "count")

library(ggplot2)
library(ggspatial)  # For adding scale bar

# Step 10: Convert raster to a data frame for ggplot
raster_df <- as.data.frame(species_raster, xy = TRUE, na.rm = TRUE)
raster_df

# Step 11: Plot the heatmap using ggplot2 without coordinates
OBC<-ggplot(data = raster_df) +
  geom_tile(aes(x = x, y = y, fill = layer)) +  # Plot the heatmap
  scale_fill_viridis_c(breaks = seq(0, max(raster_df$layer), by = 2), name = NULL) +  # Set legend breaks to increments of 2
  coord_fixed() +  # Keep aspect ratio 1:1
  ggtitle("Ordinal Best-case Scenario") +  # Add title
  geom_sf(data = richtersveld, fill = NA, color = "black", size = 1) +  # Add Richtersveld boundary
  theme_minimal() +  # Use a clean theme
  theme(legend.position = "bottom",  # Position legend at the bottom
        axis.text = element_blank(),  # Remove axis text (coordinates)
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank()) +  # Remove axis ticks
  annotation_scale(location = "bl", width_hint = 0.1)+
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"))  # Add north arrow

OBC



#WORST CASE
# Load the species occurrence data (assuming it has 'longitude', 'latitude' and 'species_name')
species_data <- read.csv("Final_occurrences_clean.csv")  # Adjust path to your species data CSV
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Load the Richtersveld shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Ensure the CRS of species data and Richtersveld match
species_sf <- st_transform(species_sf, crs = st_crs(richtersveld))

# Define the species you want to include (example: "Species_A", "Species_B", "Species_C") - change for bc and wc
selected_species <- c("Adromischus marianiae (Marloth) A.Berger", "Albuca longipes Baker", "Aloe meyeri Van Jaarsv.", "Anacampseros albissima Marloth", "Antizoma miersiana Harv.", "Astridia citrina (L.Bolus) L.Bolus", "Ceropegia articulata (Aiton) Bruyns ","Ceropegia perlata (Dinter) Bruyns", "Cheilanthes namaquensis (Baker) Schelpe & N.C.Anthony", "Codon schenckii Schinz", "Crassula subacaulis Schönland & Baker f.", "Cyrtanthus herrei (F.M.Leight.) R.A.Dyer", "Cucumis rigidus E.Mey. ex Sond.", "Gymnosporia gariepensis M. Jordaan", "Helichrysum leontonyx DC.", "Hemarthria altissima (Poir.) Stapf & C.E.Hubb.", "Hermannia eenii Baker f.", "Isolepis hemiuncialis (C.B.Clarke) J.Raynal", "Lapeirousia dolomitica Dinter", "Lapeirousia littoralis Baker", "Manulea robusta Pilg.", "Mesembryanthemum pellitum Friedrich", "Oxalis sonderiana (Kuntze) T.M.Salter", "Ozoroa crassinervia (Engl.) R.Fern. & A.Fern.", "Pelargonium echinatum Curtis", "Ruschia glauca L.Bolus", "Schwantesia herrei L.Bolus", "Trachyandra aridimontana J.C.Manning", "Tylecodon reticulatus (L.f.) Toelken", "Ursinia nana DC.", "Vahlia capensis (L.f.) Thunb.") # Modify with the species of your choice

# Filter species data to include only the selected species
species_filtered <- species_sf[species_sf$X.Species %in% selected_species, ]

# Clip the filtered species data to the Richtersveld boundary
species_clipped <- st_intersection(species_filtered, richtersveld)

# Extract coordinates from species_clipped
coordinates <- st_coordinates(species_clipped)

# Convert the sf object to a data frame, including the coordinates
species_df <- as.data.frame(species_clipped)

# Add the coordinates (longitude and latitude) to the data frame
species_df$DDE <- coordinates[, 1]
species_df$DDS <- coordinates[, 2]

species_df

# Save the data frame as a CSV
write.csv(species_df, "OWCspecies_clipped.csv", row.names = FALSE)

library(sf)
library(raster)
library(ggplot2)

# Step 1: Read the CSV file containing species data with coordinates
species_data <- read.csv("OWCspecies_clipped.csv", row.names = NULL)
str(species_data)

# Step 2: Load the Richtersveld boundary shapefile
richtersveld <- st_read("richtersveld.shp")  # Adjust path to the Richtersveld shapefile

# Step 3: Convert the CSV to an sf object (species data)
species_sf <- st_as_sf(species_data, coords = c("DDE", "DDS"), crs = 4326)

# Step 4: Reproject Richtersveld to UTM Zone 34S (EPSG:32734) for meter-based resolution
richtersveld <- st_transform(richtersveld, crs = 32734)
plot(richtersveld)

# Step 5: Reproject species data to the same CRS as the Richtersveld boundary
species_sf <- st_transform(species_sf, crs = 32734)

# Step 6: Clip the species data to the Richtersveld boundary (keeping points within Richtersveld)
species_clipped <- st_intersection(species_sf, richtersveld)

# Step 7: Create a raster grid for the heatmap
species_extent <- st_bbox(richtersveld)  # Get the extent of Richtersveld

# Specify the resolution in meters (e.g., 100 meters)
raster_grid <- raster(extent(richtersveld), res = 6324.56)  # Adjust resolution as needed

# Step 8: Set the CRS of the raster grid to match the Richtersveld boundary's CRS
crs(raster_grid) <- st_crs(richtersveld)$proj4string

# Step 9: Rasterize the species points, counting occurrences
species_raster <- rasterize(st_coordinates(species_clipped), raster_grid, field = 1, fun = "count")

library(ggplot2)
library(ggspatial)  # For adding scale bar

# Step 10: Convert raster to a data frame for ggplot
raster_df <- as.data.frame(species_raster, xy = TRUE, na.rm = TRUE)
raster_df

# Step 11: Plot the heatmap using ggplot2 without coordinates
OWC<-ggplot(data = raster_df) +
  geom_tile(aes(x = x, y = y, fill = layer)) +  # Plot the heatmap
  scale_fill_viridis_c(breaks = seq(0, max(raster_df$layer), by = 2), name = NULL) +  # Set legend breaks to increments of 2
  coord_fixed() +  # Keep aspect ratio 1:1
  ggtitle("Ordinal Worst-case Scenario") +  # Add title
  geom_sf(data = richtersveld, fill = NA, color = "black", size = 1) +  # Add Richtersveld boundary
  theme_minimal() +  # Use a clean theme
  theme(legend.position = "bottom",  # Position legend at the bottom
        axis.text = element_blank(),  # Remove axis text (coordinates)
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank() ) +  # Remove axis ticks
  annotation_scale(location = "bl", width_hint = 0.1)+
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(1, "cm"), width = unit(1, "cm"))  # Add north arrow

OWC

library(patchwork)
OBC+OWC

# 6. SPEARMAN’S CORRELATION
# Load the necessary library
library(stats)

# Read the CSV file (adjust the file path as needed)
data <- read.csv("Spearman.csv")

# View the first few rows of the data to ensure it is correct
head(data)

# Assuming your CSV has two columns: Rank_X and Rank_Y
# Replace "Rank_X" and "Rank_Y" with the actual column names in your dataset

# Extract the rank columns
rank_x <- data$IUCNRANK
rank_y <- data$BCRANK
rank_Z <- data$WCRANK

# Calculate Spearman's Rank Correlation Coefficient
BC <- cor.test(rank_x, rank_y, method = "spearman")
WC <- cor.test(rank_x, rank_Z, method = "spearman")


# Print the result
cat("BC Spearman's Rank Correlation Coefficient (rs): ", BC$estimate, "\n")
cat("BC p-value: ", BC$p.value, "\n")
cat("WC Spearman's Rank Correlation Coefficient (rs): ", WC$estimate, "\n")
cat("WC p-value: ", WC$p.value, "\n")



# 7. CIRCULAR DIAGRAMS:
library(tidyverse)
library(ggplot2)

#data path
data <- read.csv("albuca.csv")

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)




#IF UNKNOWNS other than G1 or G2
# Adjust transparency for selected individuals - those with unknown score
# Assign 0.5 transparency to specific individuals, 1 to others
# Assign transparency values manually

data$alpha <- case_when( data$individual %in% c("B5") #unknowns beside G1 or G2
~ 0.1, TRUE ~ 1) # Full transparency for others

# Define custom colors for each individual, if G1 or G2 was unknown, made colour burlywood to make manually lighter, usually orange
custom_colors <- c(
  "A1" = "brown", "A2" = "brown", "B1" = "brown", "B2" = "brown", "B3" = "brown", "B4" = "brown", "B5" = "brown", "B6" = "brown",  "B7" = "brown", "C1" = "brown", "D1" ="brown", "D2" = "brown", "D3" = "brown", "D4" = "brown", "E1" = "lightgreen", "F1" = "lightgreen",  "F2" = "lightgreen", "G1" = "burlywood", "G2" = "burlywood")

p <- ggplot(data, aes(x=as.factor(id), y=value, fill= individual)) +       
  geom_bar(stat="identity") +  # Removed alpha=0.5 here
  ylim(-4,6) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) + 
  geom_hline(yintercept = seq(0, 3, by = 1), colour = "white", linewidth = 0.5) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.5, label=individual, hjust=hjust), 
            color="black", fontface="bold", alpha=0.8, size= 5, 
            angle= label_data$angle, inherit.aes = FALSE )+ 
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  scale_alpha_continuous(range = c(0.5, 1))  # Define alpha range

p

ggsave("trachyandra.png", plot = p, width = 10, height = 10, dpi = 600)




#IF NO UNKNOWNS
# Define custom colors for each individual, if G1 or G2 was unknown, made colour burlywood to make manually lighter, usually orange
custom_colors <- c(
  "A1" = "brown", "A2" = "brown", "B1" = "brown", "B2" = "brown", "B3" = "brown", "B4" = "brown", "B5" = "brown", "B6" = "brown",  "B7" = "brown", "C1" = "brown", "D1" ="brown", "D2" = "brown", "D3" = "brown", "D4" = "brown", "E1" = "lightgreen", "F1" = "lightgreen",  "F2" = "lightgreen", "G1" = "burlywood", "G2"= "burlywood")

p <- ggplot(data, aes(x=as.factor(id), y=value, fill=individual)) +       
  geom_bar(stat="identity") +  # Removed alpha from aes()
  ylim(-4, 6) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) + 
  geom_hline(yintercept = seq(0, 3, by = 1), colour = "white", linewidth = 1) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.5, label=individual, hjust=hjust), 
            color="black", fontface="bold", size=5, angle=label_data$angle, inherit.aes = FALSE) + 
  scale_fill_manual(values = custom_colors)  # Apply custom colors

p

ggsave("albuca.png", plot = p, width = 10, height = 10, dpi = 600)

#HEATMAP OUTPUTS
#ORDINAL <img width="940" height="638" alt="image" src="https://github.com/user-attachments/assets/d6b6bd37-65d5-4ec8-af5d-528687460f34" />
#ADDITIVE <img width="940" height="635" alt="image" src="https://github.com/user-attachments/assets/b6c05842-f7d2-4912-82dc-4f1bb2a68d3a" />

#CIRCULAR VULNERABILITY DIAGRAMS
<img width="311" height="423" alt="image" src="https://github.com/user-attachments/assets/fdebcfb0-7a56-4fa9-9d77-8af3bce13269" />



