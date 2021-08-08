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
load("Results/Estimates/estimate-table-params.Rdata")
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
# COEF PLOT                                                               ----
#-----------------------------------------------------------------------------#
model_colors        <- c("#7CAE00", "#F8766D", "#00BFC4")
names(model_colors) <- c("SIGACTS","ICEWS","GED")


my_theme <- theme(axis.title         = element_blank(),
                  panel.background   = element_rect(fill   = "transparent",
                                                    colour = "black",
                                                    size   = 0.3),
                  panel.border       = element_rect(colour = "black",
                                                    fill   = NA,
                                                    size   = 0.3),
                  panel.grid         = element_blank(),
                  panel.grid.major.x = element_line(linetype = "dotted",
                                                    color    = "gray70",
                                                    size     = 0.2),
                  strip.background   = element_rect(fill     = NA,
                                                    linetype = "solid",
                                                    size     = 0.3),
                  strip.text         = element_text(size = 10),
                  legend.position    = "bottom",
                  legend.direction   = "horizontal",
                  legend.title       = element_blank(),
                  legend.key         = element_rect(fill = "transparent"))



inla_cplot <- function(data, type,
                       colors,
                       drop_vars    = NULL){
  pltdat <- data %>%
    bind_rows(., .id = "model") %>%
    # mutate(model = stringr::str_remove(model, "_bias_bias")) %>%
    mutate(model = stringr::str_to_upper(model)) %>%
    filter(type == !!type,
           !variable %in% drop_vars) %>%
    group_by(model, type) %>%
    mutate(y = 1:n()) %>%
    ungroup %>%
    mutate(z = case_when(model == "ICEWS" ~  0.1,
                         model == "GED"   ~ -0.1,
                         TRUE ~ 0)) %>%
    mutate(y = y + z) %>%
    # filter(!str_detect(model, "bias")) %>%
    mutate(model = factor(model, levels = c("ICEWS","SIGACTS","GED")))

  vars    <- pltdat %>% filter(model == "SIGACTS") %>% pull(variable)


  plt <- ggplot(data = pltdat, aes(y = y, color = model)) +
    geom_point(aes(x = median)) +
    geom_errorbar(aes(xmin = lb, xmax = ub)) +
    scale_y_continuous(breaks = 1:length(vars),
                       labels = vars) +
    scale_color_manual(values = {{colors}}) +
    my_theme

  return(plt)
}

fixed <- inla_cplot(data      = params[4:5],
                    type      = "fixed",
                    color     = model_colors,
                    drop_vars = c("Ruzicka Spending (PC)","Intercept")) +
  theme(legend.position = "none") +
  labs(title    = "Conflict Models - Coefficient Estimates and 95% HPD",
       subtitle = "Fixed Effects")
fixed

hyper <- inla_cplot(data      = params[1:3],
                    type      = "hyper",
                    color     = model_colors,
                    drop_vars = c("Range")) +
  theme(legend.position = "none") +
  labs(subtitle = "GMRF Hyperparameters")
hyper

range <- inla_cplot(data      = params[1:3],
                    type      = "hyper",
                    color     = model_colors,
                    drop_vars = c("Rho","Sigma","Kappa")) +
  xlim(0,1000)
range

# LEGEND (extract from Range and update range to no legend)
plt_legend <- get_legend(range)
range      <- range + theme(legend.position = "none")

# Event - full plot
final_cplot <- plot_grid(fixed,
                         hyper,
                         range,
                         plt_legend,
                         rel_heights = c(3,1.5,.75,.3),
                         nrow  = 4,
                         align = "v",
                         axis  = "l")
final_cplot
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


# ----------------------------------- #
# Save conflict coef plot
# ----------------------------------- #
ggsave(plot     = final_cplot,
       filename = "Results/Figures/figure-cplot-conflict.png",
       width    = 6.5,
       height   = 8.5,
       dpi      = 350,
       units    = "in")
# ----------------------------------- #

#rm(list = ls())
#-----------------------------------------------------------------------------#
