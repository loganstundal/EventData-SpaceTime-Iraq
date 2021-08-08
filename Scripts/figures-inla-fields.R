#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 06, 2021
# Purpose:       INLA space-time model figures
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#
#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# ADMINISTRATIVE                                                          ----
#-----------------------------------------------------------------------------#
#---------------------------#
# Clear working environment
#---------------------------#
rm(list = ls())
#---------------------------#

#---------------------------#
# Load required packages
#---------------------------#
library(dplyr)
library(tidyr)
library(sf)
library(INLA)
library(ggplot2)
library(scales)
library(stars)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("Results/Estimates/results-spde.Rdata")
#---------------------------#

#---------------------------#
# Functions
#---------------------------#
source("Scripts/fn_inla-field.R")

field_map <- function(field, sf, lims = NULL){
  # ----------------------------------- #
  # Description
  # ----------------------------------- #
  # A local function to produce consistent field maps for all models. Takes
  # a field data frame which is the unmodified output of inla_fields() which
  # is sourced above.
  # ----------------------------------- #
  ggplot(data = field) +
    geom_raster(aes(x = x, y = y, fill = val)) +
    geom_sf(data = sf, fill = "transparent", color = "black", size = 0.1) +
    scale_fill_gradient2(low      = muted("green"),
                         high     = muted("red"),
                         mid      = "white",
                         midpoint = 0,
                         limits   = lims) +
    facet_wrap(~ time, ncol = 4, nrow = 2) +
    theme(
      panel.background = element_rect(fill  = "transparent",
                                      color = "black",
                                      size  = 0.1),
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      strip.background = element_rect(fill  = "gray80",
                                      color = "black",
                                      size  = 0.1))
}
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# FIELD DATA                                                              ----
#-----------------------------------------------------------------------------#

# Periods IDs for naming:
periods <- dat$time_id %>% unique %>% as.character

# ----------------------------------- #
# Collect field data
# ----------------------------------- #
# Conflict models
conflict_fields <- sapply(conflict, function(model){
  inla_fields(model        = model,
              mesh         = mesh,
              index        = iset,
              time_periods = periods,
              boundary     = irq0)
}, simplify = FALSE)

# Bias models
bias_fields <- sapply(bias, function(model){
  inla_fields(model        = model,
              mesh         = mesh,
              index        = iset,
              time_periods = periods,
              boundary     = irq0)
}, simplify = FALSE)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# FIELD MAPS                                                              ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Conflict fields
# ----------------------------------- #
# For comparison to SIGACTS, assign scales identical limits
sapply(conflict_fields, function(x){summary(x$mean$val)})

# SIGACTS - Field Mean
field_plt_sig <- field_map(field = conflict_fields$sigacts$mean,
                           sf    = irq0,
                           lims  = c(-3.05, 3.0)) +
  labs(title = "Conflict models: GMRF Means", subtitle = "SIGACTS")
field_plt_sig

# ICEWS - Field Mean
field_plt_ice <- field_map(field = conflict_fields$icews$mean,
                           sf    = irq0,
                           lims  = c(-3.05, 3.0)) +
  labs(title = "Conflict models: GMRF Means", subtitle = "ICEWS")
field_plt_ice

# GED - Field Mean
field_plt_ged <- field_map(field = conflict_fields$ged$mean,
                           sf    = irq0,
                           lims  = c(-3.05, 3.0)) +
  labs(title = "Conflict models: GMRF Means", subtitle = "GED")
field_plt_ged
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Save conflict field maps
# ----------------------------------- #
ggsave(plot     = field_plt_sig,
       filename = "Results/Figures/figure-field-sigacts.png",
       width    = 9.0,
       height   = 6.5,
       dpi      = 350,
       units    = "in")

ggsave(plot     = field_plt_ice,
       filename = "Results/Figures/figure-field-icews.png",
       width    = 9.0,
       height   = 6.5,
       dpi      = 350,
       units    = "in")

ggsave(plot     = field_plt_ged,
       filename = "Results/Figures/figure-field-ged.png",
       width    = 9.0,
       height   = 6.5,
       dpi      = 350,
       units    = "in")
# ----------------------------------- #


#rm(list = ls())
#-----------------------------------------------------------------------------#
