#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          September 12, 2021
# Purpose:       Correct RW Field for APSA
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
library(INLA)
library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(scales)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("Results/Estimates/results-apsa.Rdata")
load("data/Oxford_HB_2021_APSA-Data.Rdata")
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

x      <- mod_list$m4_1_spde
spde   <- spde
dims   <- c( 200, 200)

proj   <- inla.mesh.projector(spde$mesh,dims = dims)
fields <- summary(x)$random.names[summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy"]
idx    <- which(summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy")
n      <- length(fields)
t      <- length(unique(dat$time_id))
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
means <- lapply(1:t, function(j){
  r <- inla.mesh.project(proj,
                         field = x$summary.random$i$mean[1:spde$n.spde + (j-1)*spde$n.spde]);
  return(r)})
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Create a raster object from projgrid:
rast    <- expand.grid(proj$x,
                       proj$y)
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
irq1 <- irq %>% st_transform(., crs = "+proj=longlat")

irq0 <- irq1 %>%
  st_union() %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  st_set_geometry(.,"geometry")

boundary = irq0
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
res <- lapply(1:length(means), function(z){
  time_id        <- paste0("mean_", unique(dat$time_id)[z])

  tmp <- rast %>%
    mutate(val = as.vector(means[[z]])) %>%
    rename(x = 1, y = 2)

  tmp <- st_as_stars(tmp) %>%
    st_set_crs(., st_crs(boundary))
  tmp <- suppressMessages({st_crop(tmp, boundary)})

  return(tmp)
})

names(res) <- paste0("mean_", unique(dat$time_id))
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
res2 <- lapply(names(res), function(id){

  d <- as.data.frame(res[[id]], xy = TRUE) %>% drop_na(val) %>%
    mutate(time_id = id)

  return(d)
}) %>%
  bind_rows %>%
  mutate(time_id = str_remove(time_id, "mean_")) %>%
  mutate(time_id = str_replace(time_id, "h", " - Half "))
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
summary(res2$val)

rw_1_mean <- ggplot(data = res2) +
  geom_raster(aes(x = x, y = y, fill = val)) +
  geom_sf(data = irq0, color = "black",  size = 0.1, fill = NA) +
  geom_sf(data = irq1, color = "gray50", size = 0.01, fill = NA) +
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
    legend.key.width = unit(0.5, "in"),
    legend.key.height= unit(0.1, "in"),
    # legend.text      = element_text(size  = (12/.pt)),
    strip.background = element_rect(fill  = "gray80",
                                    color = "black",
                                    size  = 0.1),
    plot.margin      = unit(c(1,0,1,0), "mm"),
    # strip.text       = element_text(size  = (14/.pt))
    ) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       breaks = seq(-1.5, 1.5, 0.5),
                       guide  = guide_colorbar(frame.colour = "black",
                                               ticks.colour = "white")) +
  facet_wrap(facets = ~time_id, ncol = 4, nrow = 2)
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
ggsave(plot     = rw_1_mean,
       filename = "Results/Figures/figure-gmrf_mean.png",
       width    = 6.5,
       height   = 4.5,
       units    = "in")

file.copy(
  from = "Results/Figures/figure-gmrf_mean.png",
  to   = "Drafts/Drafts/GAST20210912/figure-gmrf_mean.png",
  overwrite = TRUE
)

#-----------------------------------------------------------------------------#
