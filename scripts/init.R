
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(broom)

# GLOBAL SETTINGS ####

## Set working directory. Set here the path of the root of the project
setwd("C:/Users/guillert/Desktop/Projets/sourdough_evolution/")

## Growth parameter estimation (std, nls, bayes)
growth_est_mthd = "nls"

## Figure output configuration ####
plot_background = "white"
plot_resolution = 300
plot_ext = "png"

strain_types_cols <- c(Bioreal = "red",
                       Hirondelle = "green",
                       Instant = "blue",
                       sourdough = "transparent")

# CUSTOM FUNCTIONS ####

## Custom function to save ggplots with custom default parameters defined above
myggsave <- function(plot = last_plot(), filename, width = 6, height = 5, 
                     bg = plot_background, 
                     dpi = plot_resolution, 
                     device = plot_ext, ...) {
  
  file = paste0(filename,".",device)
  ggsave(plot = plot, filename = file, width = width, height = height, bg = bg, dpi = 300, device = device, ...)
  
}

## Add units to the parameters name
add_units <- function(column) {
  
  tmp <- case_when(column == "co2max" ~ paste(column, "(g)"),
                   column == "vmax" ~ paste(column, "(g/h)"),
                   column == "tvmax" ~ paste(column, "(h)"),
                   column == "lag" ~ paste(column, "(h)"),
                   TRUE ~ column)
  
  return(tmp)
}


