# setwd("/Users/chris/Desktop/gr_project")

library(lubridate)

#####
get_today <- function() {
  today() %>% as.character %>% gsub(pattern = "-", replacement = "")
}

#####
savePlotPDF <- function(plot, output_dir, output_file, width_val = 25, height_val = 22) {
  system(paste("mkdir", "-p", output_dir))
  output_filepath <- file.path(output_dir, paste0(output_file, ".pdf"))
  pdf(file = output_filepath, width = width_val, height = height_val)
  print(plot)
  dev.off()
  message(" > Plot (pdf) saved in ", output_filepath)
}

#####
savePlotPNG <- function(plot, output_dir, output_file, width_val = 25, height_val = 22, res_val = 1200) {
  system(paste("mkdir", "-p", output_dir))
  output_filepath <- file.path(output_dir, paste0(output_file, ".png"))
  png(file = output_filepath, width = width_val, height = height_val, units = "in", res = res_val)
  print(plot)
  dev.off()
  message(" > Plot (png) saved in ", output_filepath)
}

#####
savePlotEPS <- function(plot, output_dir, output_file, width_val = 25, height_val = 22) {
  system(paste("mkdir", "-p", output_dir))
  output_filepath <- file.path(output_dir, paste0(output_file, ".eps"))
  setEPS()
  postscript(file = output_filepath, width = width_val, height = height_val)
  print(plot)
  dev.off()
  message(" > Plot (eps) saved in ", output_filepath)
}