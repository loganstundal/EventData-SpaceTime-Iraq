#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 08, 2021
# Purpose:       GMRF Counterfactual analysis
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#      Calculate a prediction of the expected conflict events per capita on a
#      dense grid of Iraq
#
#      This script is a work in progress. Likely better if predictions are done
#      in the inla() estimation. Although this greatly increases computation
#      time, it will save any headaches associated with data alignment.
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
library(sf)
library(stars)
library(raster)
library(ggplot2)
library(INLA)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("Results/Estimates/results-spde.Rdata")
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# CREAD PREDICTION GRID                                                   ----
#-----------------------------------------------------------------------------#


irqr <- st_rasterize(sf = irq0)
irqr <- as.data.frame(irqr, xy = T)
irqr <- rasterFromXYZ(irqr, crs = "+proj=longlat")
in_irq <- which(!is.na(getValues(irqr)))

reference.coordinates <- coordinates(irqr)[in_irq,]


#make these into points and extract covariates for prediction grid
pred.points <- SpatialPoints(reference.coordinates, proj4string = crs(irqr))
covs <- list.files('covariates/', full.names = T) %>% stack()
pred.covs <- raster::extract(covs, pred.points, df=T)


#remake the A matrix for prediction
Aprediction <- inla.spde.make.A(mesh = mesh, loc = reference.coordinates);
dim(Aprediction)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
## using results from Model1
model = conflict$sigacts
## recall:: formula<-y ~ -1 +Intercept + f(spatial.field, model=spde) + Access + Elevation + EVI + LST_day

# Covariates for prediction points
Access<- pred.covs$Access
Elevation <-  pred.covs$Elevation
EVI <- pred.covs$EVI
LST_day <-  pred.covs$LST_day

#create the spatial structure
sfield_nodes <- model$summary.random$i['mean']
# field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1])
field <- (A %*% as.data.frame(sfield_nodes)[, 1])

#make empty matrix to fill predictions
pred <- matrix(NA, nrow = dim(A)[1], ncol = 1)

## Calculate Predicted values using regression formula
pred <- model$summary.fixed['Intercept', 'mean'] +
  model$summary.fixed['Access', 'mean'] * Access +
  model$summary.fixed['Elevation', 'mean'] * Elevation +
  model$summary.fixed['EVI', 'mean'] * EVI +
  model$summary.fixed['LST_day', 'mean'] * LST_day +
  field




# write results as a raster
x <- as.matrix(reference.coordinates)
z <- as.matrix(results)
pr.mdg.out <- rasterFromXYZ(cbind(x, z))

plot(pr.mdg.out, main = 'Prediction outside INLA')


#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# ANALYSIS                                                                ----
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
#save.image()
#rm(list = ls())
#-----------------------------------------------------------------------------#
