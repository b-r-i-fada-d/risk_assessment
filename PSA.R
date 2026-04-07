#### --- Set Up --- ####

library("tidyverse")
library("classInt")
library("ggplot2")
library("terra")
library("purrr")
library("sf")

setwd("C:/Users/boconnor/OneDrive - Marine Institute/Documents/Rwd/Risk_assessment")

rmax <- read_csv("ch3_rmax.csv")

#### --- Productivity --- ####

# --- calculate 3 Jenks breaks
breaks <- classIntervals(rmax$rmax, n = 3, style = "jenks")

# findCols returns 1 for the lowest rmax and 3 for the highest
# Since low rmax = high risk, subtract from 4 to flip
rmax <- rmax %>% mutate(Productivity_Score = 4 - findCols(breaks))

#### --- Susceptibility --- ####

# file paths
species_folder <- "species_rasters_thresholded" 
species_files  <- list.files(species_folder, pattern = "\\.tif$", full.names = TRUE)

# load gears
gear_list <- list(
  "otter"    = rast("fishing2025_rasters/otter.tif"),
  "gillnet"  = rast("fishing2025_rasters/gillnet.tif"),
  "beam" = rast("fishing2025_rasters/beam.tif")
)

all_results <- list()

# --- Calculate species-gear overlap loop

for(gear_name in names(gear_list)){
  current_gear <- gear_list[[gear_name]]
  fishing_mask <- current_gear > 0 # create a binary mask of fishing presence 
  for(spec_path in species_files) {
    spec_rast <- rast(spec_path)
    spec_name <- tools::file_path_sans_ext(basename(spec_path))
    # ensure species and gear align perfectly (resampling if necessary)
    if(!compareGeom(spec_rast, fishing_mask, stopOnError = FALSE)) {
      spec_rast <- resample(spec_rast, fishing_mask, method = "near")
    }
    total_pixels <- global(spec_rast, "sum", na.rm = TRUE)$sum # counts all pixels where species is present
    overlap_rast <- spec_rast * fishing_mask # finds where fish & gear overlap
    overlap_pixels <- global(overlap_rast, "sum", na.rm = TRUE)$sum
    overlap_pct <- (overlap_pixels / total_pixels) * 100 # what % of species' dist overlaps with gear gear
    all_results[[paste(gear_name, spec_name, sep="_")]] <- data.frame(
      Species = spec_name,
      Gear = gear_name,
      Overlap_Percent = as.numeric(overlap_pct)
    )
  }
}

# save as table
overlap <- bind_rows(all_results) %>%
  pivot_wider(names_from = Gear, values_from = Overlap_Percent)

# --- Jenks

# turn matrix back into a single vector of all overlap values
all_overlaps <- c(overlap$beam, overlap$gillnet, overlap$otter)

# calculate the 3 Jenks breaks for Susceptibility
breaks_vector <- classIntervals(all_overlaps, n = 3, style = "jenks")$brks

# function to apply specific breaks to df
assign_s_score <- function(x, brks_to_use) {
  findCols(classIntervals(x, n = 3, style = "fixed", fixedBreaks = brks_to_use))
}

# create final df
susceptibility <- overlap %>%
  mutate(
    S_beam    = assign_s_score(beam, breaks_vector),
    S_gillnet = assign_s_score(gillnet, breaks_vector),
    S_otter   = assign_s_score(otter, breaks_vector)
  )

psa <- full_join(rmax, susceptibility)

#### --- PSA --- ####

# --- Calculate vulnerability score for each gear type
# subtract 1 from the scores so that the 'origin' is 0
psa <- psa %>%
  mutate(
    V_beam = sqrt((Productivity_Score - 1)^2 + (S_beam - 1)^2),
    V_gillnet = sqrt((Productivity_Score - 1)^2 + (S_gillnet - 1)^2),
    V_otter  = sqrt((Productivity_Score - 1)^2 + (S_otter - 1)^2)
  )

psa <- psa %>%
  # Replace NAs in the Vulnerability columns with 0
  # This assumes NA means "No overlap, therefore No Risk"
  mutate(across(starts_with("V_"), ~replace_na(., 0))) %>%
  
  # Now the math will work perfectly
  mutate(
    V_Overall_Max = pmax(V_beam, V_gillnet, V_otter, na.rm = TRUE),
    V_Overall_Avg = rowMeans(dplyr::select(., V_beam, V_gillnet, V_otter), na.rm = TRUE)
  )

# write_csv(psa, "PSA_scores.csv")

#### --- Graphing, in progress --- ####

library("ggrepel")

plot_data <- psa %>%
  dplyr::select(Species, Productivity_Score, S_beam, S_gillnet, S_otter) %>%
  pivot_longer(cols = starts_with("S_"), 
               names_to = "Gear", 
               values_to = "S_score") %>%
  mutate(Gear = gsub("S_", "", Gear)) 


grid <- expand.grid(P = seq(1, 3, length.out = 100), 
                    S = seq(1, 3, length.out = 100))
grid$V <- sqrt((grid$P - 1)^2 + (grid$S - 1)^2)

ggplot() +
  geom_contour_filled(data = grid, aes(x = P, y = S, z = V), alpha = 0.2, breaks = c(0, 1, 2, 3)) +
  scale_fill_brewer(palette = "YlOrRd", name = "Vulnerability Index") +
  geom_point(data = plot_data, aes(x = Productivity_Score, y = S_score), size = 3, color = "black") +
  geom_text_repel(data = plot_data, aes(x = Productivity_Score, y = S_score, label = Species), 
                  size = 3, max.overlaps = 10) +
  facet_wrap(~Gear) +
  scale_x_continuous(breaks = 1:3, limits = c(1, 3)) +
  scale_y_continuous(breaks = 1:3, limits = c(1, 3)) +
  labs(title = "Productivity-Susceptibility Analysis by Gear Type",
       x = "Productivity Score (1 = High Productivity, 3 = Low)",
       y = "Susceptibility Score (1 = Low Overlap, 3 = High)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))