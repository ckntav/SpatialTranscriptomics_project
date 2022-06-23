setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(tidyverse)

workdir <- "/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project"

library_id <- paste0("sample", c("A1", "B1", "C1", "D1"))
molecule_h5 <- file.path(workdir, "output", library_id, "outs/molecule_info.h5")
cloupe_file <- file.path(workdir, "output", library_id, "outs/cloupe.cloupe")
spatial_folder <- file.path(workdir, "output", library_id, "outs/spatial")
tissue <- c("Prostate", "LympNode", "Prostate", "Prostate")

df <- data.frame(library_id = library_id,
                 molecule_h5 = molecule_h5,
                 cloupe_file = cloupe_file,
                 spatial_folder = spatial_folder,
                 tissue = tissue) %>% as_tibble


aggrcsv_filepath <- file.path("input", "ABCD_aggrcsv.csv")
write_csv(df, file = aggrcsv_filepath)            
