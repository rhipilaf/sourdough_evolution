
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(purrr)
library(broom)
library(forcats)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggpmisc)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(sf)
library(rlist)
library(data.table)
library(lubridate)
library(kableExtra)
library(latex2exp)
library(ggnewscale)
library(scales)


# GLOBAL SETTINGS ####

## Set working directory. Set here the path of the root of the project
setwd("C:/Users/guillert/Desktop/Projets/sourdough_evolution/")


# FIGURE OUTPUT SETTINGS ####

plot_background = "white"
plot_resolution = 300
plot_ext = "png"

## Colors ####

### Palette from https://sashamaps.net/docs/resources/20-colors/
my_cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                  '#ffffff', '#000000')

### For commercial strains
strain_types_cols <- c(Bioreal = alpha("red", 0.7),
                       Hirondelle = alpha("green", 0.7),
                       Instant = alpha("blue", 0.7),
                       sourdough = ("transparent"))


# CUSTOM FUNCTIONS ####

## Abbreviates estimates names

eqt_abbrev <- readxl::read_xlsx("data/eqtable_abbrev.xlsx", col_types = "text")
abbrev <- function(x, data) {
  stringi::stri_replace_all_regex(x,
                                  pattern=data$long,
                                  replacement=data$short,
                                  vectorize=FALSE)
}


convnames <- function(x, data, from, to="std") {
  
  if(!all(c(from,to) %in% names(data)))
    "from and to arguments must be names of data."
  
  unname(unlist(data[,to])[match(x, unlist(data[,from]))])
  
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

## Annotate significance of p.values
p_signif <- function(p.values) {
  
  tmp <- case_when(p.values < 1 & p.values > 0.1 ~ "",
                   p.values <= 0.1 & p.values > 0.05 ~ ".",
                   p.values <= 0.05 & p.values > 0.01 ~ "*",
                   p.values <= 0.01 & p.values > 0.001 ~ "**",
                   p.values <= 0.001 ~ "")
  
  return(tmp)
  
}

## All versions of a model
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



## To get modes of a statistical distribution

get_modes <- function(x, adjust, signifi, from, to) {
  
  den <- density(x, kernel=c("gaussian"),adjust=adjust,from=from,to=to)
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  den.sign <- sign(s.1$y)
  a<-c(1,1+which(diff(den.sign)!=0))
  b<-rle(den.sign)$values
  df<-data.frame(a,b)
  df = df[which(df$b %in% -1),]
  modes<-s.1$x[df$a]
  density<-s.0$y[df$a]
  
  df2<-data.frame(modes,density)
  df2$sig<-signif(df2$density,signifi)
  df2<-df2[with(df2, order(-sig)), ] 
  
  df<-as.data.frame(df2 %>% 
                      mutate(m = min_rank(desc(sig)) ) %>% #, count = sum(n)) %>% 
                      group_by(m) %>% 
                      summarize(a = paste(format(round(modes,2),nsmall=2), collapse = ',')) %>%
                      spread(m, a, sep = ''))
  colnames(df)<-paste0("m",1:length(colnames(df)))
  
  return(df)
  
}

## Function from the package spMisc

make.filenames <- function (s, replacement = "_", allow.space = TRUE) {
  replacement <- as.character(replacement)
  if (grepl(replacement, "\\\\/\\:\\*\\?\\\"\\<\\>\\|") == 
      TRUE | nchar(replacement) != 1) {
    stop(sprintf("Replacement symbol '%s' is not allowed.", replacement))
  }
  s <- gsub("[\\\\/\\:\\*\\?\\\"\\<\\>\\|]", replacement, s)
  if (allow.space == FALSE) 
    s <- gsub("[[:space:]]", replacement, s)
  return(s)
}

