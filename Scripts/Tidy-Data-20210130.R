#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          January 30, 2021
# Purpose:       Tidy - Initial data tidy of Silverman (IO 2020) replication
#                       data (+ ICEWS and GED from Ben) into spatial panels.
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#         # Modify code to only make use of SIGACTS (non-criminal) data. Per Ben
#           ICEWS & GED constructions exclude criminal event codes. (SIG_1)
#
#         # In line with above am dropping SIGACTS which includes criminal events
#           from output. Therefore, ONLY retaining SIGACTS NON-CRIMINAL.
#
#-----------------------------------------------------------------------------#



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
library(lubridate)
library(sf)

#---------------------------#
# Set working directory
#---------------------------#
setwd("c:/users/logan/googledrive/umn/research/ra_john/2021 - Time Series - RINLA")

#---------------------------#
# Load data
#---------------------------#
monthly <- read_csv("data/IRQmonthly.csv")

iraq <- st_read("data/iraq_humdata_shp/irq_admbnda_adm2_cso_20190603.shp") %>%
  select(ADM1_EN, ADM2_EN) %>%
  mutate_if(is.factor, as.character)


#-----------------------------------------------------------------------------#
# TIDY DATA (Monthly) ---------------------------------------------------------

# ----------------------------------- #
# Drop January 2004 - per Ben email of Dec. 13, this month only includes final week of January.
# ----------------------------------- #
monthly <- monthly %>%
  filter(month != 1 | year != 2004)
# ----------------------------------- #

# ----------------------------------- #
# Consolidate district_name for merge with Iraq map data, summarize new duplicates based on
# convention in "descriptions" below
# ----------------------------------- #

# Note - Abu Ghraib, Choman, and Tarmia are being aggregated into their
# respective districts based on google maps coordinates these are:
# "Al-Kadhmiyah" for Abu Ghraib and Tarmia and "Rawanduz" for Choman.

monthly <- monthly %>%
  mutate(district_name = case_when(district_name == "Abu Ghraib"       ~ "Al-Kadhmiyah",
                                   district_name == "Al-Mejar Al-Kabi" ~ "Al-Mejar Al-Kabir",
                                   district_name == "Al-Rifa'i"        ~ "Al-Rifai",
                                   district_name == "Amara"            ~ "Al-Amara",
                                   district_name == "Amedi"            ~ "Al-Amadiya",
                                   district_name == "Choman"           ~ "Rawanduz",
                                   district_name == "Darbandihkan"     ~ "Derbendikhan",
                                   district_name == "Fao"              ~ "Al-Faw",
                                   district_name == "Hamza"            ~ "Al-Hamza",
                                   district_name == "Qal'at Saleh"     ~ "Qalat Saleh",
                                   district_name == "Tarmia"           ~ "Al-Kadhmiyah",
                                   TRUE ~ district_name)) %>%
  group_by(year, month, district_name) %>%
  summarize(half = half[1],
            across(c(POP, contains("cc"), a_of_batt, SIGACT, SIG_1, icews_reb, ged_reb), sum),
            across(c(dis_usprt, cmoc, starts_with("halfy"), starts_with("su_v")), max),
            across(c(starts_with("p_"), urate, pop_den, pct_urban), mean))

vdesc = c("FORMAT"                  = "Aggregation Function -- Variable Type -- Variable short description -- Long escription",
          "half"                    = "NA   -- Identifier -- half-year ID",
          "year"                    = "NA   -- Identifier -- Identifier - year",
          "month"                   = "NA   -- Identifier -- Identifier - month",
          "district_name"           = "NA   -- Identifier -- Identifier - district name",

          "POP"                     = "SUM  -- Count      -- Population, total",
          "coalitioncc"             = "SUM  -- Count      -- Collateral damage incidents (count) attributable to coalition forces. ('Iraqi civilian casualties' (Silverman, p. 857)) ",
          "insurgentcc"             = "SUM  -- Count      -- Collateral damage incidents (count) attributable to insurgent forces.",
          "sectariancc"             = "SUM  -- Count      -- Collateral damage incidents (count) attributable to sectarian forces.",
          "unknowncc"               = "SUM  -- Count      -- Collateral damage incidents (count) attributable to unknown forces.",
          "a_of_batt"               = "SUM  -- Count      -- Coalition troop strength - (count, continuous?) - battalion count? range [0-15]",
          "SIGACT"                  = "SUM  -- Count      -- Sum of all SIGACTS (all attacks) events in district-month",
          "SIG_1"                   = "SUM  -- Count      -- Sum of all SIGACTS (excluding criminal) in district-month. Used to construct p_S1, Silverman's main DV.",
          "icews_reb"               = "SUM  -- Count      -- Sum of all ICEWS events in district-month ~ 'material violence incidents, Rebel/Insurgent > (civilians + gov)' (Ben, email of Dec. 13 2020)",
          "ged_reb"                 = "SUM  -- Count      -- Sum of all GED   events in district-month ~ 'violent incidents, Rebel/Insurgent > (civilians + gov)' (Ben, email of Dec. 13 2020)",

          "dis_usprt"               = "MAX  -- Binary     -- PRT presence - binary [0,1], 'Provincial Reconstruction Team, a measure of local Coalition expertise' (Silverman, p. 861)",
          "cmoc"                    = "MAX  -- Binary     -- CMOC presence - binary [0,1], Civil-military operations center - 'a center usually established by a military force for coordinating civil-military operations in an area of operatoins.' (Wikipedia) (Data collected by Silverman)",
          "halfyr[1:11]"            = "MAX  -- Binary     -- Half-year dummy indicator",
          "su_vh[1:11]"             = "MAX  -- Binary     -- ''a set of interactions between governorate-level vote shares for Sunni Arab parties in the 2005 parliamentary elections and half years meant to pick up broad sectarian shifts such as the Sunni Awakening.' (Silverman, p. 864)",

          "p_S1"                    = "MEAN -- Rate       -- SIGACT (noncriminal) events per 1000 population. Silverman's main DV. -- NOTE: based on my calc. likely per 1e5 pop for WEEKLY, per 1e3 only for half-year",
          "p_spentruzicka"          = "MEAN -- Rate       -- USAID funds per capita (via Marla Ruzicka Iraq War Victims Fund program) - to 'provide civilians in-kind vocational training and livelihood assistance rather than cash asfter suffering harm by Coalition forces (Silverman, p.856)'",
          "p_spentcptotal"          = "MEAN -- Rate       -- CERP (brigate-level officer distributed) funds (via DOD). Condolence - 'compensation for injury, death, or property damage by Coalition forces.' (Silverman, p. 856)",
          "p_spentcerpsmall_noncp"  = "MEAN -- Rate       -- 'noncondolence small CERP spending (Silverman, p. 861)' per capita",
          "p_spentusaid_nonruzicka" = "MEAN -- Rate       -- 'non-Ruzicka USAID spending (Silverman, p. 861)' per capita",
          "p_spentrecon_noncerp"    = "MEAN -- Rate       -- Non-CERP reconstruction spending per capita",
          "p_spentcpall"            = "MEAN -- Rate       -- Combined spending per capita",
          "p_icews"                 = "MEAN -- Rate       -- ICEWS events per 100,000 population",
          "p_ged"                   = "MEAN -- Rate       -- GED   events per 100,000 population",
          "urate"                   = "MEAN -- Rate       -- Unemployment rate - continuous [0-1]",
          "pop_den"                 = "MEAN -- Rate       -- Population density - contniuous [0-14] ~ per square (km?)",
          "pct_urban"               = "MEAN -- Rate       -- Percent urban population - continuous [0 - 1.1]")


# ----------------------------------- #
# Rename variables for consistent convention
# ----------------------------------- #
monthly <- monthly %>%
  rename(Half_Yr_ID    = half,
         Year          = year,
         Month         = month,
         District      = district_name,
         Pop_ct        = POP,
         Pop_urban_pct = pct_urban,
         Pop_den       = pop_den,
         Unemp_rate    = urate,

         CD_Coalition  = coalitioncc,
         CD_Insurgent  = insurgentcc,
         CD_Sectarian  = sectariancc,
         CD_Unknown    = unknowncc,

         Troop_Strgth  = a_of_batt,

         # SIGACT_ct_all = SIGACT,
         SIGACT_ct     = SIG_1,
         SIGACT_pc     = p_S1,
         ICEWS_ct      = icews_reb,
         ICEWS_pc      = p_icews,
         GED_ct        = ged_reb,
         GED_pc        = p_ged,

         PRT           = dis_usprt,
         CMOC          = cmoc       ) %>%

  # Dropping SIGACTS with criminal events here
  select(everything(), -SIGACT)


# ----------------------------------- #
# Monthly lags (levels and first diffs.), drop lag NAs
# ----------------------------------- #
monthly <- monthly %>%
  arrange(District, Year, Month) %>%
  group_by(District) %>%
  mutate(across(starts_with(c("SIG","ICE","GED")), .fns = list(lag      = ~ dplyr::lag(., 1),
                                                               diff     = ~ . - dplyr::lag(., 1),
                                                               diff_lag = ~ dplyr::lag(. - dplyr::lag(., 1), 1)))) %>%
  ungroup() %>%
  tidyr::drop_na(contains("_lag")) %>%
  select(Year, Month, District, Half_Yr_ID,
         starts_with(c("SIG","ICE","GED")),
         starts_with("CD_"),
         Troop_Strgth,
         PRT,
         CMOC,
         starts_with("Pop"),
         Unemp_rate,
         starts_with("p_"),
         everything())
# ----------------------------------- #



#-----------------------------------------------------------------------------#
# DATA FRAME AGGREGATES -------------------------------------------------------


# ----------------------------------- #
# Time Series (monthly) [plotting only]
# ----------------------------------- #
d_ts <- monthly %>%
  group_by(Year, Month) %>%
  summarize(across(ends_with("_ct"), sum)) %>%
  ungroup() %>%
  mutate(Date = ymd(sprintf("%s-%s-%s", Year, Month, 15))) %>%
  select(Date, ends_with("_ct"), -Pop_ct) %>%
  pivot_longer(cols      = ends_with("_ct"),
               names_to  = "Var",
               values_to = "Value")


# ----------------------------------- #
# Cross section [plotting only]
# ----------------------------------- #
d_cs <- monthly %>%
  group_by(District) %>%
  summarize(across(ends_with("_ct"), sum),
            across(ends_with("_pc"), mean)) %>%
  ungroup() %>%
  select(District,
         ends_with(c("_ct","pc")),
         -Pop_ct)

# ----------------------------------- #
# Panel (district-year)
# ----------------------------------- #
d_pn_yr <- monthly %>%
  group_by(Year, District) %>%
  summarize(across(ends_with("_ct"), sum),
            across(ends_with("_pc"), mean)) %>%
  ungroup() %>%
  select(Year, District,
         ends_with(c("_ct","pc")),
         -Pop_ct)

# ----------------------------------- #
# Panel (district-month) + time fixed-effect
# ----------------------------------- #
d_pn_mo <- monthly %>%
  mutate(Time_FE = paste(Year, Month, sep = "-"))

rm(monthly)



#-----------------------------------------------------------------------------#
# SPATIAL DATA MERGING --------------------------------------------------------


# ----------------------------------- #
# Merge map geometry to cross-section and panel dataframes
# ----------------------------------- #
iraq <- iraq %>% rename(Governorate = ADM1_EN,
                        District    = ADM2_EN)

d_cs <- d_cs %>%
  left_join(., iraq, by = "District") %>%
  select(Governorate, everything()) %>%
  as.data.frame() %>%
  st_set_geometry(., "geometry")

d_pn_mo <- d_pn_mo %>%
  left_join(., iraq, by = "District") %>%
  select(Year, Month, Governorate, everything()) %>%
  as.data.frame() %>%
  st_set_geometry(., "geometry")

d_pn_yr <- d_pn_yr %>%
  left_join(., iraq, by = "District") %>%
  select(Year, Governorate, everything()) %>%
  as.data.frame() %>%
  st_set_geometry(., "geometry")


rm(iraq)


# SAVE ------------------------------------------------------------------------
save.image("Data/Data-All-Tidy-20210131.Rdata")
rm(list = ls())

