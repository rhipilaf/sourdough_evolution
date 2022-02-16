
library(magrittr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(broom)
library(forcats)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpmisc)
library(FactoMineR)
library(factoextra)
library(sf)

# GLOBAL SETTINGS ####

## Set working directory. Set here the path of the root of the project
setwd("C:/Users/guillert/Desktop/Projets/sourdough_evolution/")

## Growth parameter estimation method (std, nls, bayes)
growth_est_mthd = "nls" # parameter not used yet


# FIGURE OUTPUT SETTINGS ####

plot_background = "white"
plot_resolution = 300
plot_ext = "png"

## Colors ####
strain_types_cols <- c(Bioreal = alpha("red", 0.7),
                       Hirondelle = alpha("green", 0.7),
                       Instant = alpha("blue", 0.7),
                       sourdough = ("transparent"))


# CUSTOM FUNCTIONS ####

## Abbreviates estimates names

abb_data <- readxl::read_xlsx("data/abbrev.xlsx", col_types = "text")
abbrev <- function(x, data) {
  stringi::stri_replace_all_regex(x,
                         pattern=data$long,
                         replacement=data$short,
                         vectorize=FALSE)
}

## To retrieve errors and warnings (from https://stackoverflow.com/a/4952908/2161065)
factory <- function(fun) {
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(fun(...), error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(res, warn=warn, err=err)
  }
}
  

## Custom function to save ggplots with custom default parameters defined above ####
myggsave <- function(plot = last_plot(), filename, width = 6, height = 5, 
                     bg = plot_background, 
                     dpi = plot_resolution, 
                     device = plot_ext, ...) {
  
  file = paste0(filename,".",device)
  ggsave(plot = plot, filename = file, width = width, height = height, bg = bg, dpi = 300, device = device, ...)
  
}

kb_style <- function(kable) {
  kable %>%
    kableExtra::kable_styling() %>%
    kableExtra::scroll_box(width = "100%", height = "500px")
}

## Add units to the parameters name ####
add_units <- function(column) {
  
  tmp <- case_when(column == "co2max" ~ paste(column, "(g)"),
                   column == "vmax" ~ paste(column, "(g/h)"),
                   column == "tvmax" ~ paste(column, "(h)"),
                   column == "lag" ~ paste(column, "(h)"),
                   TRUE ~ column)
  
  return(tmp)
}


signif <- function(p.values) {
  
  tmp <- case_when(p.values < 1 & p.values > 0.1 ~ "",
            p.values <= 0.1 & p.values > 0.05 ~ ".",
            p.values <= 0.05 & p.values > 0.01 ~ "*",
            p.values <= 0.01 & p.values > 0.001 ~ "**",
            p.values <= 0.001 ~ "")
  return(tmp)
  
}


mod_vers <- function(full.model, force = NULL) {

  # In this version of the function, interaction terms must be explicitely specified. For example, the b * a should be written b + a + b:a.
  # The 'force' argument must be a vector giving the position of forced elements in the formula.

  
  model_elements <- strsplit(as.character(full.model)[3], split = " + ", 
                             fixed = T)[[1]]
  model_elements_forced <- model_elements[force]
  model_elements_totest <- setdiff(model_elements, model_elements_forced)
  
  model_versions <- sapply(1:length(model_elements_totest), 
                            FUN = function(x) combn(model_elements_totest, x, 
                                                    simplify = F)) %>% 
    flatten() %>% 
    map(~ paste(.x, collapse = " + ")) %>% 
    map(~ paste0("value ~ ", .x)) %>% 
    map(~ gsub("~ \\(", "~ 1 + \\(", .x)) %>% #Adds a '1 + ' for the intercept when only random effects are specified
    flatten_chr()
  
  model_versions_final <- paste0(c(model_versions, "value ~ 1"), paste0(c("", model_elements_forced), collapse = " + "))
  
  return(model_versions_final)
  
}


