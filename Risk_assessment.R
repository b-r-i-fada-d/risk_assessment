#### --- Set Up --- ####

library("tidyverse")
library("terra")
library("sf")
library("raster")
library("rnaturalearth")
library("sf")
library("ggplot2")

setwd("C:/Users/boconnor/OneDrive - Marine Institute/Documents/Rwd/Risk_assessment")

gillnet <- raster("fishing2025_rasters/gillnet.tif")
otter <- raster("fishing2025_rasters/otter.tif")
beam <- raster("fishing2025_rasters/beam.tif")

countries <- ne_countries(scale = "medium", returnclass = "sf")

### --- Risk assessments using threshold for noise in SDM output, single species --- ####

# set threshold of minimum predicted presence probability from jSDM
threshold <- 0.005

# load species raster
jsdm_raster <- raster("species_rasters/Raja brachyura.tif")

# apply threshold species presence
species_mask <- jsdm_raster
species_mask[species_mask < threshold] <- NA
species_mask[species_mask >= threshold] <- 1

# resample fishing effort rasters to jSDM raster resolution
gillnet_res <- resample(gillnet, jsdm_raster, method = "bilinear")
otter_res <- resample(otter, jsdm_raster, method = "bilinear")
beam_res <- resample(beam, jsdm_raster, method = "bilinear")

# replace NAs in gear rasters with 0
gillnet_res[is.na(gillnet_res)] <- 0
otter_res[is.na(otter_res)] <- 0
beam_res[is.na(beam_res)] <- 0

# mask gear rasters by species presence
risk_g <- gillnet_res*species_mask
risk_o <- otter_res*species_mask
risk_b <- beam_res*species_mask

# normalize each gear to 0-1 for plotting
# I did this without normalising and the map was the same but the legend was different
# think about this
normalize <- function(x) {
  rng <- range(values(x), na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
  }

risk_g_norm <- normalize(risk_g)
risk_o_norm <- normalize(risk_o)
risk_b_norm <- normalize(risk_b)

# convert to dfs
df_g <- as.data.frame(rasterToPoints(risk_g_norm)); df_g$gear <- "Gillnet"
df_o <- as.data.frame(rasterToPoints(risk_o_norm)); df_o$gear <- "Otter"
df_b <- as.data.frame(rasterToPoints(risk_b_norm)); df_b$gear <- "Beam"

risk_df <- rbind(df_g, df_o, df_b)
colnames(risk_df) <- c("lon", "lat", "risk", "gear")
risk_df$gear <- factor(risk_df$gear, levels = c("Gillnet", "Otter", "Beam"))

# keep only non-zero risk cells for plotting
risk_df <- risk_df[risk_df$risk > 0, ]

# load countries
countries <- ne_countries(scale = "medium", returnclass = "sf")

# --- Plot --- #
ggplot(risk_df) +
  geom_tile(aes(x = lon, y = lat, fill = risk)) + 
  geom_sf(data = countries, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,1)) +
  facet_wrap(~gear) +
  coord_sf(xlim = c(-20,0), ylim = c(45,60)) +
  labs(title = "Bycatch risk by gear",
       fill = "Risk (0-1)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 10))

#### --- Risk assessment with normalised risk across gears --- ####

species_folder <- "species_rasters/"
species_files <- list.files(species_folder, pattern="\\.tif$", full.names = TRUE)

# set threshold of minimum predicted presence probability from jSDM
threshold <- 0.005

# min-max normalisation: changes each gear to 0-1 for plotting
# I did this without normalising and the map was the same but the legend was different
# think about this

normalize <- function(x){
  rng <- range(values(x), na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
  }

# Function to convert raster to dataframe (with gear risk normalised)
raster_to_df <- function(r, gear_name){  # r is raster object
# if anything produces all NAs, fill in with 0s
# this is needed because deepwater species won't interact with the beam trawl
# and it will break the loop  
  if (all(is.na(values(r)))){           
    return(data.frame(lon = numeric(0),
                      lat = numeric(0),
                      risk = numeric(0),
                      gear = factor(levels=c("Gillnet", "Otter", "Beam"))))
    } else {
      df <- as.data.frame(rasterToPoints(r))
      df$gear <- gear_name
      return(df)
    }
}


# empty data frame to store all risks
all_risk_df <- data.frame()

# loop over species
for(species_file in species_files){
  species_name <- tools::file_path_sans_ext(basename(species_file))
  jsdm_raster <- raster(species_file) # load species raster
  species_mask <- jsdm_raster # threshold
  species_mask[species_mask < threshold] <- NA
  species_mask[species_mask >= threshold] <- 1
# resample fishing effort rasters to jSDM raster resolution
  gillnet_res <- resample(gillnet, jsdm_raster, method = "bilinear")
  otter_res <- resample(otter, jsdm_raster, method = "bilinear")
  beam_res <- resample(beam, jsdm_raster, method = "bilinear")
# replace NAs in gear rasters with 0
  gillnet_res[is.na(gillnet_res)] <- 0
  otter_res[is.na(otter_res)] <- 0
  beam_res[is.na(beam_res)] <- 0
# mask gear rasters by species presence
  risk_g <- gillnet_res*species_mask
  risk_o <- otter_res*species_mask
  risk_b <- beam_res*species_mask
# normalize each gear to 0-1
  risk_g_norm <- normalize(risk_g)
  risk_o_norm <- normalize(risk_o)
  risk_b_norm <- normalize(risk_b)
# convert to data frames
  df_g <- raster_to_df(risk_g_norm, "Gillnet")
  df_o <- raster_to_df(risk_o_norm, "Otter")
  df_b <- raster_to_df(risk_b_norm, "Beam")
  species_risk_df <- rbind(df_g, df_o, df_b)
  colnames(species_risk_df) <- c("lon", "lat", "risk", "gear")
  species_risk_df$gear <- factor(species_risk_df$gear, 
                                 levels = c("Gillnet","Otter","Beam"))
  species_risk_df$species <- species_name
# keep only non-zero risk cells
  species_risk_df <- species_risk_df[species_risk_df$risk > 0, ]
  all_risk_df <- rbind(all_risk_df, species_risk_df)
}

#### --- Take 2: not normalised --- ####

species_folder <- "species_rasters/"
species_files <- list.files(species_folder, pattern="\\.tif$", full.names = TRUE)

# set threshold of minimum predicted presence probability from jSDM
threshold <- 0.005

# Function to convert raster to dataframe (keep absolute risk)
raster_to_df_nonnorm <- function(r, gear_name){
  if(all(is.na(values(r)))){
    df <- data.frame(
      lon = NA,
      lat = NA,
      risk = 0,
      gear = factor(gear_name, levels = c("Gillnet","Otter","Beam"))
    )
    return(df)
  } else {
    df <- as.data.frame(rasterToPoints(r))
    colnames(df)[1:3] <- c("lon","lat","risk")  # rename x, y, layer
    df$gear <- factor(gear_name, levels = c("Gillnet","Otter","Beam"))
    return(df)
  }
}

# empty data frame to store all risks
all_risk_df_nonnorm <- data.frame()

# loop over species
for(species_file in species_files){
  species_name <- tools::file_path_sans_ext(basename(species_file))
  jsdm_raster <- raster(species_file)
  # create presence mask (1 = presence, NA = absence)
  species_mask <- jsdm_raster
  species_mask[species_mask < threshold] <- NA
  species_mask[species_mask >= threshold] <- 1
  # resample fishing effort rasters to jSDM raster resolution
  gillnet_res <- resample(gillnet, jsdm_raster, method = "bilinear")
  otter_res  <- resample(otter, jsdm_raster, method = "bilinear")
  beam_res   <- resample(beam, jsdm_raster, method = "bilinear")
  # replace NAs with 0 (no effort)
  gillnet_res[is.na(gillnet_res)] <- 0
  otter_res[is.na(otter_res)] <- 0
  beam_res[is.na(beam_res)] <- 0
  # mask gear rasters by species presence
  risk_g <- gillnet_res * species_mask
  risk_o <- otter_res  * species_mask
  risk_b <- beam_res   * species_mask
  # convert to dataframes (keep absolute magnitude of gear-specific pressure)
  df_g <- raster_to_df_nonnorm(risk_g, "Gillnet")
  df_o <- raster_to_df_nonnorm(risk_o, "Otter")
  df_b <- raster_to_df_nonnorm(risk_b, "Beam")
  # combine gear dataframes
  species_risk_df <- rbind(df_g, df_o, df_b)
  species_risk_df$species <- species_name
  all_risk_df_nonnorm <- rbind(all_risk_df_nonnorm, species_risk_df)
}

#### --- Plot faceted by species (rows) and gear (columns) --- ####

# filter only valid lon/lat and positive risk
nonnorm_plot <- all_risk_df_nonnorm[!is.na(all_risk_df_nonnorm$lon) & all_risk_df_nonnorm$risk > 0, ]

p <- ggplot(all_risk_df) +
  geom_tile(aes(x = lon, y = lat, fill = risk)) +
  geom_sf(data = countries, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma") +
  facet_grid(species ~ gear) +
  coord_sf(xlim = c(-20, 0), ylim = c(45, 60)) +
  labs(title = "Bycatch Risk by Species and Gear", fill = "Risk") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 10))

ggsave("output/bycatch_risk_species_gear_norm.png", plot=p, # change here
       width = 10, height = 4 + length(unique(all_risk_df$species))*1.5,  # adjust height by #species
       dpi=300)

#### --- Cumulative risk --- ####

# --- Compute cumulative risk per species ---
cumulative_risk_df <- all_risk_df %>%
  group_by(species, lon, lat) %>%
  summarize(cum_risk = sum(risk, na.rm=TRUE), .groups="drop")  # sum across gears

# --- Normalize cumulative risk to 0-1 for plotting ---
cumulative_risk_df <- cumulative_risk_df %>%
  group_by(species) %>%
  mutate(cum_risk_norm = (cum_risk - min(cum_risk, na.rm=TRUE)) / 
           (max(cum_risk, na.rm=TRUE) - min(cum_risk, na.rm=TRUE))) %>%
  ungroup()

#### --- Plot cumulative risk per species --- ####

species_order <- c(
  "Deania calcea",
  "Centrophorus squamosus",
  "Centroscymnus coelolepis",
  "Dipturus oxyrinchus",
  "Centroselachus crepidater",
  "Scymnodon ringens",
  "Centroscyllium fabricii",  
  "Etmopterus princeps",
  "Etmopterus spinax",
  "Hexanchus griseus",
  "Galeus melastomus",
  "Leucoraja circularis",
  "Raja clavata",
  "Squalus acanthias",
  "Raja montagui",
  "Leucoraja naevus",
  "Apristurus aphyodes",
  "Apristurus manis",
  "Raja undulata",
  "Raja brachyura"
)

cumulative_risk_df$species <- factor(cumulative_risk_df$species,
                                     levels = species_order)

# 4x5 layout
p <- ggplot(cumulative_risk_df) +
  geom_tile(aes(x = lon, y = lat, fill = cum_risk_norm)) +
  geom_sf(data = countries, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,1)) +
  facet_wrap(~species, ncol = 4) +  # 4 columns, species follow factor order
  coord_sf(xlim = c(-20, 0), ylim = c(45, 60)) +
  labs(title = "Cumulative bycatch risk across all gears",
       fill = "Cumulative Risk (0-1)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 10))

ggsave("output/cumulative_bycatch_risk.png",
       plot = p,        
       width = 14,
       height = 16,
       dpi = 300)

#### --- Statistics --- ####

# --- Summary statistics per species
summary <- all_risk_df %>%
  group_by(species) %>%
  summarize(
    n_cells = n(),
    mean_risk = mean(risk),
    max_risk = max(risk),
    median_risk = median(risk),
    high_risk_frac = mean(risk > 0.5)
  )

# --- Gear-specific comparisons
gear_risks <- all_risk_df %>%
  group_by(species, gear) %>%
  summarize(total_risk = sum(risk)) %>%
  mutate(percent_contrib = total_risk / sum(total_risk) * 100)

#### --- Graphing --- #####

# Create a data frame mapping species to habitat
species_group <- data.frame(
  Species = c(deepwater_species, shallow_species),
  Habitat = c(rep("Deepwater", length(deepwater_species)),
              rep("Shallow", length(shallow_species)))
)

ggplot(gear_risks, 
       aes(x = factor(species), 
           y = percent_contrib, 
           fill = gear)) +
  geom_col() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 10)) + 
  scale_fill_viridis_d(option = "mako", direction = -1,
                       begin = 0.15, end = 0.85
                       ) +
  labs(y = "Percentage of contribution to overall risk",
       fill = "Gear") +
  geom_col(color = "grey20", linewidth = 0.2) +
  # scale_fill_brewer(palette = "Blues") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        plot.margin = margin(t = 5, r = 5, b = 25, l = 5))

# 
# ggsave("output/cumulative_bycatch_risk.png",
#        plot = p,        
#        width = 14,
#        height = 16,
#        dpi = 300)


# --- Spatial overlap / hotspots
# Identify cells where multiple species have high risk (hotspots for conservation action).
# Could calculate correlation of risk maps between species to see if high-risk areas coincide.
# # Example: correlation between two species
# cor(subset(all_risk_df, species == "SpeciesA")$risk,
#     subset(all_risk_df, species == "SpeciesB")$risk)

# --- Metrics like species richness in high-risk areas (number of species with risk > threshold in each cell).

# --- Threshold sensitivity
# Check how results change with different presence probability thresholds (like your 0.005).
# Could present a table of total area at risk vs. threshold.

# --- Cumulative risk across species
# Sum risk across all species per cell 
# Could quantify proportion of study area that has cumulative risk above certain levels.


# --- Spatial statistics: Moran’s I, hotspot analysis, or clustering of high-risk areas.
 
# --- Typical Tables/Figures
# Table: Risk summary per species + gear.
# Bar chart: Gear contribution per species.

# --- Heatmap: Species richness or high-risk hotspot overlap.
