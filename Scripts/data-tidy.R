#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          June 29, 2021
# Purpose:       Tidy Silverman (2021) replication data in preparation for
#                both STAR and GMRF extensions.
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
library(skimr)
library(janitor)

#---------------------------#
# Load data
#---------------------------#
irq_month  <- read_csv("Data/IRQmonthly.csv")
irq_halfyr <- read_csv("Data/IRQhalfYear.csv")
irq        <- read_sf("data/esoc_v3/gis boundary files/iraq_district_boundaries_utm.shp")
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# TIDY DATA                                                               ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Monthly
# ----------------------------------- #
irq_month <- irq_month %>%
  select(-starts_with(c("halfyr", "su_vh"))) %>%

  # Generate appropriate first differences and drop NAs on DV
  group_by(district_name) %>%
  mutate(across(starts_with(c("p_S1", "p_icews", "p_ged")),
                .fns = list(d   = ~ . - dplyr::lag(., 1)))) %>%

  # Create a temporal lag of DVs levels for future use
  mutate(across(starts_with(c("p_S1", "p_icews", "p_ged")),
                .fns = ~ dplyr::lag(.,1), .names="{.col}_lag")) %>%

  # Create first differences of non-endogenous vars:
  mutate(across(c(p_spentcptotal, p_spentruzicka, coalitioncc, insurgentcc,
                  p_spentcerpsmall_noncp, p_spentusaid_nonruzicka, a_of_batt),
                .fns = list(d = ~ . - dplyr::lag(., 1)))) %>%

  # Drop missing values on key variables
  drop_na(p_S1_d_lag, a_of_batt, dis_usprt) %>%

  # create time-id and convert ID vars to factors
  mutate(time_id = paste(year, month, sep = ".")) %>%
  # convert to factors
  mutate(district_name = as_factor(district_name),
         month = as_factor(month),
         half  = as_factor(half),
         year  = as_factor(year),
         time_id = as_factor(time_id)) %>%

  # Select variables and clean-up:
  select(district_name, month, half, year, time_id,
         contains(c("_S1", "SIG", "_icews", "_ged")),
         everything()) %>%
  clean_names() %>%    # janitor
  ungroup()
# ----------------------------------- #


# ----------------------------------- #
# Half-year
# ----------------------------------- #
irq_halfyr <- irq_halfyr %>%
  select(-num_range("halfyr", 1:10),
         -num_range("su_vh", 1:10)) %>%

  # Generate appropriate first differences and drop NAs on DV
  group_by(district_name) %>%
  mutate(across(starts_with(c("p_S1", "p_icews", "p_ged")),
                .fns = list(d = ~ . - dplyr::lag(., 1)))) %>%

  # Create a temporal lag of DVs levels for future use
  mutate(across(starts_with(c("p_S1", "p_icews", "p_ged")),
                ~ dplyr::lag(.,1), .names="{.col}_lag")) %>%

  # Create first differences of non-endogenous vars:
  mutate(across(c(p_spentcptotal, p_spentruzicka, coalitioncc, insurgentcc,
                  p_spentcerpsmall_noncp, p_spentusaid_nonruzicka, a_of_batt),
                .fns = list(d = ~ . - dplyr::lag(., 1)))) %>%

  # Drop missing values on key variables
  drop_na(p_S1_d_lag) %>%

  # create time-id and convert ID vars to factors
  rename(time_id = halfyr) %>%
  # convert to factors
  mutate(district_name = as_factor(district_name),
         time_id       = as_factor(time_id)) %>%

  # Select variables and clean-up:
  select(district_name, time_id,
         contains(c("_S1", "SIG", "_icews", "_ged")),
         everything(),
         -c("_merge","halfyrid")) %>%
  clean_names() %>%    # janitor
  ungroup()
# ----------------------------------- #


# ----------------------------------- #
# Tidy spatial data
# ----------------------------------- #
irq <- irq %>%
  rename(district_name = ADM3NAME) %>%
  select(district_name, AREA_KM2, PERIM_KM) %>%
  clean_names()
# ----------------------------------- #


# ----------------------------------- #
# Review tidy data
# ----------------------------------- #
# skim(irq_month)
# skim(irq_halfyr)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# RECONCILE NAME DIFFERENCES                                              ----
#-----------------------------------------------------------------------------#

# Adjust spatial data names and reorder
irq <- irq %>%
    mutate(district_name = case_when(
      district_name == "Adhamiya"      ~ "Al-Adhamiya",
      district_name == "Akre"          ~ "Aqra",
      district_name == "Al-Ba'aj"      ~ "Al-Baaj",
      district_name == "Al-Ka'im"      ~ "Al-Kaim",
      district_name == "Al-Mahawil"    ~ "Al-Mahaweel",
      district_name == "Al-Musayab"    ~ "Al-Mussyab",
      district_name == "Al-Na'maniya"  ~ "Al-Namaniya",
      district_name == "Al Resafa"     ~ "Al-Risafa",
      district_name == "Al Sadr"       ~ "Al-Thawra",
      district_name == "Ba'quba"       ~ "Baquba",
      district_name == "Baiji"         ~ "Beygee",
      district_name == "Baladrooz"     ~ "Baladruz",
      district_name == "Basrah"        ~ "Al-Basrah",
      district_name == "Dahuk"         ~ "Duhok",
      district_name == "Diwaniya"      ~ "Al-Diwaniya",
      district_name == "Falluja"       ~ "Al-Falluja",
      district_name == "Halabja"       ~ "Halabcha",
      district_name == "Hashimiya"     ~ "Al-Hashimiya",
      district_name == "Hatra"         ~ "Al-Hatra",
      district_name == "Hilla"         ~ "Al-Hilla",
      district_name == "Karkh"         ~ "Al-Karkh",
      district_name == "Kerbala"       ~ "Kerbela",
      district_name == "Khadamiya"     ~ "Al-Kadhmiyah",
      district_name == "Koisnjaq"      ~ "Koysinjaq",
      district_name == "Kufa"          ~ "Al-Kufa",
      district_name == "Kut"           ~ "Al-Kut",
      district_name == "Mada'in"       ~ "Al-Mada'in",
      district_name == "Mahmoudiya"    ~ "Al-Mahmoudiya",
      district_name == "Makhmur"       ~ "Makhmour",
      district_name == "Mergasur"      ~ "Al-Zibar",
      district_name == "Mosul"         ~ "Al-Mosul",
      district_name == "Najaf"         ~ "Al-Najaf",
      district_name == "Nassriya"      ~ "Al-Nasiriya",
      district_name == "Penjwin"       ~ "Panjwin",
      district_name == "Ramadi"        ~ "Al-Ramadi",
      district_name == "Shatt Al-Arab" ~ "Shat Al-Arab",
      district_name == "Soran"         ~ "Rawanduz",
      district_name == "Sulaymaniya"   ~ "Al-Sulaymaniyah",
      district_name == "Sumel"         ~ "Sumail",
      district_name == "Tilkaif"       ~ "Tilkaef",
      district_name == "Tooz"          ~ "Tooz Khurmato",
      TRUE ~ district_name)) %>%
  arrange(district_name)

# Reorder half-year and monthly data by district_name
irq_halfyr <- irq_halfyr %>% arrange(time_id, district_name)
irq_month  <- irq_month  %>% arrange(time_id, district_name)
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SPATIAL WEIGHTS MATRICES                                               ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Create cross-section W
# ----------------------------------- #
# Create spatialpolygonsdataframe
sp_dat <- as_Spatial(irq)

# Create neighbors list
nbs <- spdep::poly2nb(pl        = sp_dat,
                      row.names = sp_dat$district_name,
                      queen     = TRUE)

# Create matrix of cross-section weights
nbs_cs_w <- spdep::nb2mat(neighbours = nbs,
                          style      = "W")
colnames(nbs_cs_w) <- rownames(nbs_cs_w)

# Create a list w for cross-section (may not need)
nbs_cs_lw <- spdep::mat2listw(nbs_cs_w, style = "W")
# ----------------------------------- #


# ----------------------------------- #
# Create district-half-year weights
# ----------------------------------- #
# Create a space-time weights matrix
nbs_hy_w <- do.call(what = magic::adiag,
                    args = replicate(n = length(unique(irq_halfyr$time_id)),
                                     expr     = nbs_cs_w,
                                     simplify = FALSE))

ids <- expand.grid(rownames(nbs_cs_w), unique(irq_halfyr$time_id))
rownames(nbs_hy_w) <- apply(ids, 1, paste, collapse=".")
colnames(nbs_hy_w) <- rownames(nbs_hy_w)

# Create a list w object
nbs_hy_lw <- spdep::mat2listw(nbs_hy_w, style = "W")

# Convert to sparse matrix for space saving
nbs_hy_w  <- as(nbs_hy_w, "sparseMatrix")
# ----------------------------------- #


# ----------------------------------- #
# Create district-month weights
# ----------------------------------- #
nbs_mo_w <- do.call(what = magic::adiag,
                    args = replicate(n = length(unique(irq_month$time_id)),
                                     expr      = nbs_cs_w,
                                     simplify = FALSE))

ids <- expand.grid(rownames(nbs_cs_w), unique(irq_month$time_id))
rownames(nbs_mo_w) <- apply(ids, 1, paste, collapse=".")
colnames(nbs_mo_w) <- rownames(nbs_mo_w)

# Create a list w object
nbs_mo_lw <- spdep::mat2listw(nbs_mo_w, style = "W")

# Convert to sparse matrix for space saving
nbs_mo_w  <- as(nbs_mo_w, "sparseMatrix")
# ----------------------------------- #


# ----------------------------------- #
# W Eigenvalues
# ----------------------------------- #
# Construct eigenvalues of spatial weights for speedier estimation:

e_cs <- eigen(nbs_cs_w, only.values = TRUE)$values
e_hy <- rep(e_cs, each = length(unique(irq_halfyr$time_id)))
e_mo <- rep(e_cs, each = length(unique(irq_month$time_id)))
# ----------------------------------- #


# ----------------------------------- #
# Clean up
# ----------------------------------- #

# Organize spatial elements into named list:
sp_wts <- list("nbs"   = nbs,
               "cs_lw" = nbs_cs_lw,
               "cs_w"  = nbs_cs_w,
               "cs_ev" = e_cs,
               "hy_lw" = nbs_hy_lw,
               "hy_w"  = nbs_hy_w,
               "hy_ev" = e_hy,
               "mo_lw" = nbs_mo_lw,
               "mo_w"  = nbs_mo_w,
               "mo_ev" = e_mo)

rm(ids, sp_dat, e_cs, e_hy, e_mo,
   nbs, nbs_cs_lw, nbs_cs_w, nbs_hy_lw, nbs_hy_w, nbs_mo_lw, nbs_mo_w)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SPATIAL COORDINATES (INLA)                                              ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Extract latitude and longitude values (district centroids)
# ----------------------------------- #
latlongs <- irq %>%
  st_centroid() %>%
  st_transform(crs = "+proj=longlat") %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(lat = X,
         long = Y) %>%
  bind_cols(., irq) %>%
  select(lat, long, district_name)
# ----------------------------------- #


# ----------------------------------- #
# Merge to data
# ----------------------------------- #
irq_halfyr <- left_join(x  = irq_halfyr,
                        y  = latlongs,
                        by = "district_name")

irq_month  <- left_join(x  = irq_month,
                        y  = latlongs,
                        by = "district_name")

rm(latlongs)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



# ______________________________----



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
# save.image(file = "data/Oxford_HB_2021_APSA-Data.Rdata")
rm(list = ls())
#-----------------------------------------------------------------------------#
