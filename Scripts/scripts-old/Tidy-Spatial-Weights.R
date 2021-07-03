#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          January 31, 2021
# Purpose:       Tidy - Spatial weights construction for (initial) discrete space-time analysis
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#       - Naming script: Tidy-Spatial-Weights.R for convention as will require
#         long-lat in the future, Spatial-Coords.R
#
#-----------------------------------------------------------------------------#



# ADMINISTRATIVE --------------------------------------------------------------

#---------------------------#
# Clear working environment
#---------------------------#
rm(list = ls())

#---------------------------#
# Load required packages
#---------------------------#
library(tidyverse)
library(sf)
library(spdep)
library(magic)

#---------------------------#
# Set working directory
#---------------------------#
setwd("c:/users/logan/googledrive/umn/research/ra_john/2021 - Time Series - RINLA")

#---------------------------#
# Load data
#---------------------------#
load("Data/Data-All-Tidy-20210131.Rdata")



#-----------------------------------------------------------------------------#
# SPATIAL WEIGHTS -------------------------------------------------------------


# ----------------------------------- #
# W - Construct cross-section matrix
# ----------------------------------- #
# Create a grid for one month. Units don't change over time (thank goodness!)
# and so Will be identical across time.
g <- d_cs %>%  select(District)
g <- as_Spatial(g)

W1 <- spdep::poly2nb(g)
W1 <- spdep::nb2mat(W1, style = "W")
# ----------------------------------- #


# ----------------------------------- #
# W - Construct space-time matrix.
# ----------------------------------- #
# Slide W along diagonal n times:
nrow(d_pn_mo) / nrow(W1)

W <- do.call(adiag, replicate((nrow(d_pn_mo) / nrow(W1)), W1, simplify=FALSE))
W <- as(W, "sparseMatrix")

colnames(W) <- with(d_pn_mo, paste(Year, Month, District, sep = "-"))
rownames(W) <- colnames(W)
# ----------------------------------- #


# ----------------------------------- #
# W - Construct space-time nbs list
# ----------------------------------- #
W_list <- mat2listw(W)
# ----------------------------------- #



#-----------------------------------------------------------------------------#
# W eigenvalues ---------------------------------------------------------------

e       <- eigen(W1, only.values = TRUE)$values
W_eigen <- rep(e, each = (nrow(d_pn_mo) / nrow(W1)))



#-----------------------------------------------------------------------------#
# Clean up --------------------------------------------------------------------
rm(list = setdiff(ls(), c("W", "W_list", "W_eigen")))


# SAVE ------------------------------------------------------------------------
save.image("Data/Data-Spatial-Weights-20210131.Rdata")
rm(list = ls())
