#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          June 29, 2021
# Purpose:       STAR Extension of Silverman (2021) + ICEWS, GED comparatives
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
# Load required packages
#---------------------------#
library(tidyverse)
library(sf)
library(spatialreg)
library(texreg)

#---------------------------#
# Load data
#---------------------------#
load("data/Oxford_HB_2021_APSA-Data.Rdata")

#---------------------------#
# Local functions
#---------------------------#
var_names <- function(models){

  # This function takes a models list and formats the coefficient names
  # locating "Phi" (temporal lag) and "Rho" (spatial lag - if applicable)
  # at the top of the list and the intercept at the bottom.

  # For use as an input to screenreg() or texreg() in the custom.coef.map
  # parameter.

  nms <- unique(sapply(models, function(x){
    names(coef(x))}) %>% c %>% as.list)

  # Set names
  names(nms) <- nms

  # Replace lags with consistent phi label and correct intercept label
  nms[str_detect(nms, "_lag")]        <- "Phi"
  if("rho" %in% nms){nms[str_detect(nms, "rho")] <- "Rho"}
  nms[str_detect(nms, "(Intercept)")] <- "Intercept"

  # Move phi and rho to top and intercept to bottom
  nms <- c(nms[nms == "Rho"], nms[nms != "Rho"])
  nms <- c(nms[nms == "Phi"], nms[nms != "Phi"])
  nms <- c(nms[nms != "Intercept"], nms[nms == "Intercept"])

  # Return
  return(nms)
}

lm_test <- function(models, weights){

  # This function returns a tidied data frame of Lagrange multiplier
  # tests of spatial dependence on an input list of linear models.

  res <- lapply(models, function(x){
    tmp <- spdep::lm.LMtests(model = x,
                             listw = weights,
                             test  = "all")

    tmp <- summary(tmp)$results %>% rownames_to_column(var = "Test")
    tmp <- sapply(tmp, unlist) %>% as_tibble %>%
      mutate(Stat = as.numeric(statistic),
             Pval = as.numeric(p.value)) %>%
      select(-parameter, -statistic, -p.value)
      # pivot_longer(!Test, names_to = "Stat", values_to = "Value") %>%
      # mutate(Value = as.numeric(Value))

    return(tmp)
  })

  res <- bind_rows(res, .id = "DV")

  return(res)
}
#---------------------------#
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# MODEL SETUP                                                             ----
#-----------------------------------------------------------------------------#

# Orig DV: p_S1_d
# My DVL : p_S1_d_lag

# Original:
f <- formula(. ~ p_spentcptotal_d + p_spentruzicka_d +
               coalitioncc_d + insurgentcc_d +
               p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
               a_of_batt_d + cmoc + dis_usprt)

# f <- formula(. ~ p_spentcptotal + p_spentruzicka +
#                coalitioncc + insurgentcc +
#                p_spentcerpsmall_noncp + p_spentusaid_nonruzicka +
#                a_of_batt + cmoc + dis_usprt)
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# LINEAR MODELS                                                           ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Half-year models
# ----------------------------------- #
mods_linear_hy <- list(
  "SIGACT" = lm(update(f, p_s1_d    ~ p_s1_d_lag    + .), data = irq_halfyr),
  "ICEWS"  = lm(update(f, p_icews_d ~ p_icews_d_lag + .), data = irq_halfyr),
  "GED"    = lm(update(f, p_ged_d   ~ p_ged_d_lag   + .), data = irq_halfyr)
)

screenreg(l = mods_linear_hy,
          custom.coef.map = var_names(mods_linear_hy))
# ----------------------------------- #


# ----------------------------------- #
# Monthly models
# ----------------------------------- #
mods_linear_mo <- list(
  "SIGACT" = lm(update(f, p_s1_d    ~ p_s1_d_lag    + .), data = irq_month),
  "ICEWS"  = lm(update(f, p_icews_d ~ p_icews_d_lag + .), data = irq_month),
  "GED"    = lm(update(f, p_ged_d   ~ p_ged_d_lag   + .), data = irq_month)
)

screenreg(l = mods_linear_mo,
          custom.coef.map = var_names(mods_linear_mo))
# ----------------------------------- #


# ----------------------------------- #
# Lagrange multiplier tests
# ----------------------------------- #
lagrange <- list(
  "halfyr" = lm_test(mods_linear_hy, weights = sp_wts$hy_lw),
  "month"  = lm_test(mods_linear_mo, weights = sp_wts$mo_lw)
)

lagrange2 <- lagrange %>%
  map(function(x){x %>% filter(str_detect(Test, "RLM"))} ) %>%
  bind_rows(.id = "Agg") %>%
  pivot_longer(cols = Stat:Pval,
               names_to = "stat",
               values_to = "val") %>%
  mutate(val  = format(val, digits = 0, nsmall = 3,
                       scientific = F, trim = T)) %>%
  mutate(val  = case_when(stat == "Pval" ~ sprintf("(%s)",val), TRUE ~ val),
         Agg  = case_when(Agg == "halfyr" ~ "Half-year", TRUE ~ "Month"),
         Test = case_when(Test == "RLMerr" ~ "Error", TRUE ~ "Lag")) %>%
  pivot_wider(id_cols = c(Agg, Test, stat),
              names_from = DV,
              values_from = val) %>%
  select(-stat)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SPATIAL MODELS - MLE                                                    ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Half-year models
# ----------------------------------- #
mods_spatial_hy <- list(
  "SIGACT" = lagsarlm(formula = update(f, p_s1_d ~ p_s1_d_lag + .),
                      data    = irq_halfyr,
                      listw   = sp_wts$hy_lw,
                      method  = "eigen",
                      control = list(pre_eig = sp_wts$hy_ev)),
  "ICEWS"  = lagsarlm(formula = update(f, p_icews_d ~ p_icews_d_lag + .),
                      data    = irq_halfyr,
                      listw   = sp_wts$hy_lw,
                      method  = "eigen",
                      control = list(pre_eig = sp_wts$hy_ev)),
  "GED"    = lagsarlm(formula = update(f, p_ged_d ~ p_ged_d_lag + .),
                      data    = irq_halfyr,
                      listw   = sp_wts$hy_lw,
                      method  = "eigen",
                      control = list(pre_eig = sp_wts$hy_ev))
)

screenreg(l = mods_spatial_hy,
          custom.coef.map = var_names(mods_spatial_hy))
# ----------------------------------- #


# ----------------------------------- #
# Monthly models
# ----------------------------------- #
# Code difference due to long computation time for monthly models.
# This allows for individual model estimation.

SIGACT  <- lagsarlm(formula = update(f, p_s1_d ~ p_s1_d_lag + .),
                    data    = irq_month,
                    listw   = sp_wts$mo_lw,
                    method  = "eigen",
                    control = list(pre_eig = sp_wts$mo_ev))

ICEWS   <- lagsarlm(formula = update(f, p_icews_d ~ p_icews_d_lag + .),
                    data    = irq_month,
                    listw   = sp_wts$mo_lw,
                    method  = "eigen",
                    control = list(pre_eig = sp_wts$mo_ev))

GED     <- lagsarlm(formula = update(f, p_ged_d ~ p_ged_d_lag + .),
                    data    = irq_month,
                    listw   = sp_wts$mo_lw,
                    method  = "eigen",
                    control = list(pre_eig = sp_wts$mo_ev))

mods_spatial_mo <- list(
  "SIGACT" = SIGACT,
  "ICEWS"  = ICEWS,
  "GED"    = GED
)
rm(SIGACT, ICEWS, GED)

screenreg(l = mods_spatial_mo,
          custom.coef.map = var_names(mods_spatial_mo))
# ----------------------------------- #


# ----------------------------------- #
# Create tidy spatial models for table
# ----------------------------------- #
mods_spatial_hy_tidy <- lapply(mods_spatial_hy, function(x){
  tmp <- broom::tidy(lmtest::coeftest(x))
  tmp <- createTexreg(coef.names = tmp$term,
                      coef       = tmp$estimate,
                      se         = tmp$std.error,
                      pvalues    = tmp$p.value,
                      gof.names  = c("Num. obs."),
                      gof        = c(832))
})

mods_spatial_mo_tidy <- lapply(mods_spatial_mo, function(x){
  tmp <- broom::tidy(lmtest::coeftest(x))
  tmp <- createTexreg(coef.names = tmp$term,
                      coef       = tmp$estimate,
                      se         = tmp$std.error,
                      pvalues    = tmp$p.value,
                      gof.names  = c("Num. obs."),
                      gof        = c(6032))
})
# ----------------------------------- #
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SPATIAL MODELS - BAYES                                                  ----
#-----------------------------------------------------------------------------#
# library(brms)
#
# f <- formula(. ~ p_spentcptotal + p_spentruzicka +
#                # coalitioncc + insurgentcc +
#                # p_spentcerpsmall_noncp + p_spentusaid_nonruzicka +
#                # a_of_batt + cmoc + dis_usprt +
#                # sar(nbs_hy_w, type = "error") +
#                sar(hy_lw, type = "lag"))
#
# m <- brm(formula = update(f, p_s1_d ~ p_s1_d_lag + .),
#          data    = irq_halfyr,
#          data2   = list(hy_lw = sp_wts$hy_lw),
#          chains  = 2,
#          iter    = 200,
#          warmup  = 100,
#          cores   = 2)
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Variable names for models
# ----------------------------------- #
names_linear <- var_names(mods_linear_hy)
names_linear[4:12] <- c("Condolence spending",
                        "Ruzicka spending",
                        "Coalition collateral damage",
                        "Insurgent collateral damage",
                        "Other small CERP spending",
                        "Other USAID spending",
                        "Coalition troop strength",
                        "CMOC presence",
                        "PRT presence")

names_spatial <- var_names(mods_spatial_hy)
names_spatial[5:13] <- c("Condolence spending",
                        "Ruzicka spending",
                        "Coalition collateral damage",
                        "Insurgent collateral damage",
                        "Other small CERP spending",
                        "Other USAID spending",
                        "Coalition troop strength",
                        "CMOC presence",
                        "PRT presence")
# ----------------------------------- #

# save(list = ls()[str_detect(ls(), "mods_|lagrange|var_|names_")],
#      file = "data/Oxford_HB_2021_APSA-Models-Discrete.Rdata")
rm(list = ls())
#-----------------------------------------------------------------------------#
