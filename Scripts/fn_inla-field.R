#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 06, 2021
# Purpose:       functions to extract field datt from estimated INLA models
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#      9/12/2021 -- needs updating. Bad code. Aggregates fields when plotted.
#      See: figures-inla-apsa_rw_corrected.R for time corrected version.
#
# NOTE NOTE NOTE --- THIS CODE IS BAD -- INLA.FIELDS() GRABBED INCORRECT VALS
# USE CORRECTED VERSION IN "figures-inla-fields_apsa_rw_corrected.R
#
# They corrected script and function produce the same rows - so the ggplot
# is wrong due to a missing facet. BUT -- the summary() result of the field
# value across all years differs!!! Why???
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# inla_fields() ~ INLA spacetime fields
#-----------------------------------------------------------------------------#

require(stars)

inla_fields <- function(model,
                        mesh,
                        index,
                        time_periods,
                        boundary,
                        dims = 200){
  # ----------------------------------- #
  # Description
  # ----------------------------------- #
  # This function takes an inla model object, mesh, and time index, and returns
  # a dataframe of field means and standard deviations at raster grid locations
  #
  # Note - increasing dims can greatly increase object size. Increase
  #        with caution.
  # ----------------------------------- #

  # Construct projector grid
  projgrid <- inla.mesh.projector(mesh       = mesh,
                                  dims       = c(dims, dims),
                                  projection = "longlat")

  # Create a raster object from projgrid:
  rast    <- expand.grid(projgrid$x,
                         projgrid$y)

  # Extract field data
  field <- lapply(1:length(time_periods), function(x){
    tp_mean <- paste0("mean_", as.character(time_periods[x]))
    tp_sd   <- paste0("sd_",   as.character(time_periods[x]))

    xm <- inla.mesh.project(
      projgrid, model$summary.random$i$mean[index$i.group == x])
    xs <- inla.mesh.project(
      projgrid, model$summary.random$i$sd[index$i.group == x])

    tmp       <- rast
    tmp[[tp_mean]] <- as.vector(xm)
    tmp[[tp_sd]]   <- as.vector(xs)

    tmp <- st_as_stars(tmp) %>%
      st_set_crs(., st_crs(boundary))
    tmp <- suppressMessages({st_crop(tmp, boundary)})
    tmp <- as.data.frame(tmp, xy = TRUE) %>% drop_na(!!tp_mean, !!tp_sd)
    return(tmp)
  })
  field <- suppressMessages({bind_cols(field)}) %>%
    rename(x = 1, y = 2) %>%
    dplyr::select(-contains("..")) %>%
    pivot_longer(.,
                 cols      = !c(x,y),
                 names_to  = "id",
                 values_to = "val") %>%
    separate(., id, into = c("type", "time"), sep = "_")

  field_mean <- field %>% filter(type == "mean")
  field_sd   <- field %>% filter(type == "sd")

  field <- list("mean" = field_mean,
                "sd"   = field_sd)

  return(field)
}
#-----------------------------------------------------------------------------#

