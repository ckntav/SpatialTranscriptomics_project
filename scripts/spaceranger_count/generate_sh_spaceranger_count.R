# setwd("/Users/chris/Desktop/SpatialTranscriptomics_project")

library(tidyverse)

sSheet <- read_tsv("input/sampleSheet_spatialTranscriptomics_prostate.txt")
area_list <- sSheet %>% pull(area)

#
header_sh <- c("#!/bin/sh",
               "#SBATCH --time=12:00:00",
               "#SBATCH --cpus-per-task=32",
               "#SBATCH --mem=128G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")
workdir <- "/home/chris11/projects/def-stbil30/chris11/SpatialTranscriptomics_project"

# arguments
transcriptome_path <- file.path(workdir, "input", "refdata-gex-GRCh38-2020-A")
fastq_folder <- file.path(workdir, "raw", "fastq")
slidefile_path <- file.path(workdir, "input", "V10J27-102.gpr")
slide <- "V10J27-102"

for (areax in area_list) {
  area_line <- sSheet %>% dplyr::filter(area == areax)
  
  idx <- area_line %>% pull(id)
  samplex <- area_line %>% pull(sample)
  imagex <- area_line %>% pull(image)
  
  message('# ', areax, " | ", idx, " | ", samplex)
  
  call_spaceranger_count <- paste("spaceranger" , "count",
                                  "--reorient-images",
                                  paste0("--id=", idx),
                                  paste0("--transcriptome=", transcriptome_path),
                                  paste0("--fastqs=", fastq_folder),
                                  paste0("--sample=", samplex),
                                  paste0("--image=", file.path(workdir, "raw", "image", imagex)),
                                  paste0("--slidefile=", slidefile_path),
                                  paste0("--slide=", slide),
                                  paste0("--area=", areax),
                                  "--localcores=32",
                                  "--localmem=128")
  
  message(call_spaceranger_count)
  # system(call_bamCoverage)
  
  header_sh_specific <- c(paste0("#SBATCH --job-name=", paste0("spaceranger_count_", areax)),
                       paste0("#SBATCH --output=", file.path(workdir, "output", "log", paste0("spaceranger_count_", areax, "_log.txt"))))
  header_sh_x <- c(header_sh, header_sh_specific)
  
  
  
  
  file_sh <- file.path("scripts/spaceranger_count/batch_sh",
                       paste0("spaceranger_count_", areax, ".sh"))
  fileConn <- file(file_sh)
  writeLines(c(header_sh_x, "\n", call_spaceranger_count), fileConn)
  close(fileConn)
}