#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          January 31, 2021
# Purpose:       Analysis - Initial discrete space-time analysis
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



# ADMINISTRATIVE --------------------------------------------------------------

#---------------------------#
# Clear working environment
#---------------------------#
rm(list = ls())

#---------------------------#
# Load required packages
#---------------------------#
library(tidyverse)
library(spatialreg)
library(sf)
library(texreg)
library(sandwich)

#---------------------------#
# Set working directory
#---------------------------#
setwd("c:/users/logan/googledrive/umn/research/ra_john/2021 - Time Series - RINLA")

#---------------------------#
# Load data
#---------------------------#
load("Data/Data-All-Tidy-20210131.Rdata")
load("Data/Data-Spatial-Weights-20210131.Rdata")

rm(d_cs, d_ts, d_pn_yr)

#-----------------------------------------------------------------------------#
# MODEL FORMATTING ------------------------------------------------------------

# ----------------------------------- #
# Static formulas
# ----------------------------------- #
fs <- list("SIGACT" = formula(SIGACT_pc ~ SIGACT_pc_lag +
                                Pop_den +
                                # Pop_urban_pct + Unemp_rate +
                                CD_Coalition + CD_Insurgent + CD_Sectarian +
                                # Troop_Strgth + PRT + CMOC +
                                p_spentruzicka + p_spentcptotal),
           "ICEWS"  = formula(ICEWS_pc ~ ICEWS_pc_lag +
                                Pop_den +
                                # Pop_urban_pct + Unemp_rate +
                                CD_Coalition + CD_Insurgent + CD_Sectarian +
                                # Troop_Strgth + PRT + CMOC +
                                p_spentruzicka + p_spentcptotal),
           "GED"    = formula(GED_pc ~ GED_pc_lag +
                                Pop_den +
                                # Pop_urban_pct + Unemp_rate +
                                CD_Coalition + CD_Insurgent + CD_Sectarian +
                                # Troop_Strgth + PRT + CMOC +
                                p_spentruzicka + p_spentcptotal))
# ----------------------------------- #


# ----------------------------------- #
# Dynamic formulas - TO DO
# ----------------------------------- #

# To do
# ----------------------------------- #


# ----------------------------------- #
# Model Name formatting
# ----------------------------------- #
m_vnames = list("SIGACT_pc_lag"     = "DV - Time Lag",
                "ICEWS_pc_lag"      = "DV - Time Lag",
                "GED_pc_lag"        = "DV - Time Lag",
                "Pop_den"           = "Pop. Density",
                "Pop_urban_pct"     = "Pop. Urban",
                "Unemp_rate"        = "Unemployment",
                "CD_Coalition"      = "CD Coalition",
                "CD_Insurgent"      = "CD Insurgent",
                "CD_Sectarian"      = "CD Sectarian",
                "Troop_Strgth"      = "Troop Strength",
                "PRT"               = "PRT",
                "CMOC"              = "CMOC",
                "p_spentruzicka"    = "PC Ruzicka",
                "p_spentcptotal"    = "PC CERP",
                "(Intercept)"       = "Intercept")



#-----------------------------------------------------------------------------#
# MODELS: NON-SPATIAL ---------------------------------------------------------


# ----------------------------------- #
# Linear models (no dynamics, no FEs)
# ----------------------------------- #
m1 <- lapply(fs, function(x){lm(formula = x,
                                data    = d_pn_mo)})

m1_ses <- lapply(m1,           function(x){sqrt(diag(vcovCL(x, cluster = ~ District, type = "HC1")))})
m1_pvs <- lapply(1:length(m1), function(x){
  2 * (1 - pnorm(abs(coef(m1[[x]]))/m1_ses[[x]]))
})

# ----------------------------------- #
# Linear models (no dynamics, time, and fixed effects)
# ----------------------------------- #
m2 <- lapply(fs, function(x){lm(formula = update(x, . ~ . + Time_FE + District),
                                data    = d_pn_mo)})

m2_ses <- lapply(m2,           function(x){sqrt(diag(vcovCL(x, cluster = ~ District, type = "HC1")))})
m2_pvs <- lapply(1:length(m2), function(x){
  2 * (1 - pnorm(abs(coef(m2[[x]]))/m2_ses[[x]]))
})


screenreg(l = c(m1,m2),
          custom.model.names = rep(c("SIGACT","ICEWS","GED"),2),
          custom.header      = list("No FEs" = 1:3, "FEs" = 4:6),
          custom.coef.map    = m_vnames,
          override.se        = c(m1_ses, m2_ses),
          override.pvalues   = c(m1_pvs, m2_pvs))



#-----------------------------------------------------------------------------#
# SPATIAL DIAGNOSTICS ---------------------------------------------------------


# ----------------------------------- #
# Spatial diagnostic tests
# ----------------------------------- #
m2_diag <- lapply(m2, function(x){
  spdep::lm.LMtests(model = x,
                    listw = W_list,
                    test  = c('LMerr','LMlag','RLMerr','RLMlag'))
})

# The warning about row-standardization is a damn lie.

m2_diag$SIGACT  # Spatial lag - robust spatial lag highly significant
m2_diag$ICEWS   # Spatial lag - robust spatial lag significant, error more significant
m2_diag$GED     # Does not support spatial lag or spatial error



#-----------------------------------------------------------------------------#
# SPATIAL LAG MODELS ----------------------------------------------------------


# ----------------------------------- #
# Spatial linear model (no FEs)
# ----------------------------------- #
# use eigenvalue method (Ord, 1975) for computational boost
m3 <- lagsarlm(formula = fs$SIGACT,
               data    = d_pn_mo,
               listw   = W_list,
               method  = "eigen",
               control = list(pre_eig  = W_eigen,
                              OrdVsign = -1))

# This ^ works! Convergence is slow, but produces parameter estimates with no
# errors or warning messages. Therefore, this is a good baseline to compare
# the custom NLS model against.



#-----------------------------------------------------------------------------#
# Marginal Effects: STAR-Lag Model --------------------------------------------


# ----------------------------------- #
# NB - I have modified this to account for the fact that I forgot the
#      temporal lag the first time I executed m3
# ----------------------------------- #

n <- nrow(d_pn_mo)

rho  <- coef(m3)['rho']      # Extract spatial lag parameter from model
# phi  <- coef(m3)['ben95t1']  # Extract temporal lag parameter from model
beta <- coef(m3)["p_spentcptotal"]   # Variable of interest
Id   <- diag(1,n,n)                        # nxn Identity matrix (NOTE - HERE n SET one cross-section!

# Estimate Long-Run-Steady-State Multiplier Matrix
# M  <- solve(Id - rho*W - phi*Id)     # LRSS_multiplier
M <- solve(Id - rho*W)                                      # Ah yes, this solve() is computationally intense. Doable. But likely not in an monte carlo.
bM <- beta * M

# Like with our SAR model, we could looks at the LRSS for one state:
# Note - this is the long-run steady state after the system reachs an
# equilibrium from a one-unit change in X.
# Here is an example using California again and our SAR map function
# sar_eff_map(unit     = "CA",
#             map      = us,
#             bM       = bM,
#             join_var = 'STUSPS',
#             breaks   = c(0.001,0.01,0.05,0.1,0.6))
# round(sort(bM[,'CA'], decreasing = T), 3)

# Thus, relative to the regression coefficient on the variable of interest
# you can see that the response in y to a unit change in x is much larger
# than the coefficient alone would imply. This response accounts for the
# temporal and spatial dynamics implied by our model.

# Effects:
dir.eff <- mean(diag(bM))
ind.eff <- (sum(bM) - sum(diag(bM))) / n
total   <- sum(bM) / n
cbind('Direct' = dir.eff, 'Indirect' = ind.eff, 'Total' = total)






#-----------------------------------------------------------------------------#
# SAVE ------------------------------------------------------------------------
# save.image()
# rm(list = ls())



