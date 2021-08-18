#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 07, 2021
# Purpose:       INLA space-time model tables
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
library(INLA)
library(kableExtra)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("Results/Estimates/results-spde.Rdata")
rm(list = setdiff(ls(), c("conflict", "bias", "spde")))
#---------------------------#

#---------------------------#
# Functions
#---------------------------#
source("scripts/fn_inla-table.R")
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# TIDY PARAMETERS                                                         ----
#-----------------------------------------------------------------------------#
conflict_params <- inla_params(model_list = conflict, spde = spde)
bias_params     <- inla_params(model_list = bias,     spde = spde)
names(bias_params) <- paste0(names(bias_params),"_bias")

params <- c(conflict_params, bias_params)

var_names <- c(
  "sigacts_l"                 = "Phi",
  "icews_l"                   = "Phi",
  "ged_l"                     = "Phi",
  "sigacts2"                  = "SIGACTS",
  "p_spentcptotal_d"          = "Condolence Spending (PC)",
  "p_spentruzicka_d"          = "Ruzicka Spending (PC)",
  "coalitioncc_d"             = "Coalition Collateral Damage",
  "insurgentcc_d"             = "Insurgent Collateral Damage",
  "p_spentcerpsmall_noncp_d"  = "Other Small CERP Spending",
  "p_spentusaid_nonruzicka_d" = "Other USAID Spending",
  "a_of_batt_d"               = "Coalition Troop Strenght",
  "cmoc"                      = "CMOC Presence",
  "dis_usprt"                 = "PRT Presence"
)


params <- lapply(params, function(x){
  x <- x %>%
    mutate(variable = recode(variable, !!!var_names))
})


save(params, file = "Results/Estimates/estimate-table-params.Rdata")
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# COEFFICIENT PLOTS                                                       ----
#-----------------------------------------------------------------------------#

conflict_tables <- inla_table(params = conflict_params, format = "html")
print(conflict_tables)
#-----------------------------------------------------------------------------#


