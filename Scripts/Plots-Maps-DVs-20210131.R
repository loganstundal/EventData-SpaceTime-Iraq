#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          January 31, 2021
# Purpose:       Producing maps of SIGACTS, ICEWS, and GED event counts
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#     Jan. 31 - Major revisions to original data tidy. Code below will need to
#               be updated to reflect streamlined variable naming convention.
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
# library(tidyverse)
# library(lubridate)
# library(sf)
# library(ggpmisc)

#---------------------------#
# Set working directory
#---------------------------#
# setwd()

#---------------------------#
# Load data
#---------------------------#
# load()

#---------------------------#
# Load functions
#---------------------------#


#-----------------------------------------------------------------------------#

tbl_ts <- as.data.frame(t(sapply(unique(dat_ts$Var),
                                 function(x){
                                   Value = dat_ts %>% filter(Var == x) %>% pull(Value)
                                   c(summary(Value),"SD" = sd(Value))
                                 }))) %>%
  rownames_to_column() %>%
  mutate(Source   = str_extract(rowname, "[^_]+"),
         Variable = rowname) %>%
  select(`Variable`,`Mean`,`Min.`,`Max.`, `SD`, `Source`)

kbl(tbl_ts %>% select(-Source),
    digits = 3,
    format = "html",
    table.attr = "style = \"color: black;\"",
    caption    = "Summary statistics, Aggregation: Month") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  pack_rows(index = c("SIGACT" = 3, "ICEWS" = 2, "GED" = 2))
rm(tbl_ts)



tbl_cs <- as.data.frame(t(sapply(dat_cs %>% select(-governorate,-district_name, -contains("brk")) %>% st_drop_geometry(),
                                 function(x){
                                   c(summary(x), "SD" = sd(x), "n>0" = sum(x>0))
                                 }))) %>%
  rownames_to_column() %>%
  mutate(Source   = str_extract(rowname, "[^_]+"),
         Variable = rowname) %>%
  select(`Variable`,`n>0`,`Mean`,`Min.`,`Max.`, `SD`, `Source`)

kbl(tbl_cs %>% select(-Source),
    digits = 3,
    format = "html",
    table.attr = "style = \"color: black;\"",
    caption    = "Summary statistics, Aggregation: District") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  pack_rows(index = c("SIGACT" = 3, "ICEWS" = 2, "GED" = 2))
rm(tbl_cs)



tbl_pn_yr <- as.data.frame(t(sapply(dat_pn_dist_yr %>% select(-year,-governorate,-district_name,-contains("brk"))
                                    %>% st_drop_geometry(),
                                    function(x){
                                      c(summary(x), "SD" = sd(x), "n>0" = sum(x>0))
                                    }))) %>%
  rownames_to_column() %>%
  mutate(Source   = str_extract(rowname, "[^_]+"),
         Variable = rowname) %>%
  select(`Variable`,`n>0`,`Mean`,`Min.`,`Max.`, `SD`, `Source`)

kbl(tbl_pn_yr %>% select(-Source),
    digits = 3,
    format = "html",
    table.attr = "style = \"color: black;\"",
    caption    = "Summary statistics, Aggregation: District-Year") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  pack_rows(index = c("SIGACT" = 3, "ICEWS" = 2, "GED" = 2))
rm(tbl_pn_yr)



tbl_pn_mo <- as.data.frame(t(sapply(dat_pn_dist_mo %>% select(-year,-month,-governorate,-district_name,-contains("brk"))
                                    %>% st_drop_geometry(),
                                    function(x){
                                      c(summary(x), "SD" = sd(x), "n>0" = sum(x>0))
                                    }))) %>%
  rownames_to_column() %>%
  mutate(Source   = str_extract(rowname, "[^_]+"),
         Variable = rowname) %>%
  select(`Variable`,`n>0`,`Mean`,`Min.`,`Max.`, `SD`, `Source`)

kbl(tbl_pn_mo %>% select(-Source),
    digits = 3,
    format = "html",
    table.attr = "style = \"color: black;\"",
    caption    = "Summary statistics, Aggregation: District-Month") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  pack_rows(index = c("SIGACT" = 3, "ICEWS" = 2, "GED" = 2))
rm(tbl_pn_mo)




ggplot(data = dat_ts %>%
         filter(str_detect(Var, pattern = "monthly_sum")),
       aes(x = date)) +
  geom_line(aes(y = Value)) +
  facet_wrap(~ Var, scales = "free") +
  labs(title = "Reported Monthly Events", y = "", x = "") +
  theme_ts






#-----------------------------------------------------------------------------#
# MAP SPECIFIC VARS -----------------------------------------------------------

# ----------------------------------- #
# Bins for choropleth maps
# ----------------------------------- #
pltbreaks = c(0, 10, 50, 100, 250, 1000, Inf)

dat_cs <- dat_cs %>%
  mutate(SIGACT_monthly_sum_brk = cut(SIGACT_monthly_sum, breaks = pltbreaks, include.lowest = TRUE),
         ICEWS_monthly_sum_brk  = cut(ICEWS_monthly_sum,  breaks = pltbreaks, include.lowest = TRUE),
         GED_monthly_sum_brk    = cut(GED_monthly_sum,    breaks = pltbreaks, include.lowest = TRUE)) %>%
  mutate(SIGACT_monthly_sum_brk = plyr::revalue(SIGACT_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         SIGACT_monthly_sum_brk = plyr::revalue(SIGACT_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")),

         ICEWS_monthly_sum_brk = plyr::revalue(ICEWS_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         ICEWS_monthly_sum_brk = plyr::revalue(ICEWS_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")),

         GED_monthly_sum_brk = plyr::revalue(GED_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         GED_monthly_sum_brk = plyr::revalue(GED_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")) )

dat_pn_dist_mo <- dat_pn_dist_mo %>%
  mutate(SIGACT_brk = cut(SIGACT, breaks = pltbreaks, include.lowest = TRUE),
         ICEWS_brk  = cut(ICEWS,  breaks = pltbreaks, include.lowest = TRUE),
         GED_brk    = cut(GED,    breaks = pltbreaks, include.lowest = TRUE)) %>%
  mutate(SIGACT_brk = plyr::revalue(SIGACT_brk, c("(250,1e+03]" = "(250,1000]")),
         SIGACT_brk = plyr::revalue(SIGACT_brk, c("(1e+03,Inf]" = "(1000+]")),

         ICEWS_brk = plyr::revalue(ICEWS_brk, c("(250,1e+03]" = "(250,1000]")),
         ICEWS_brk = plyr::revalue(ICEWS_brk, c("(1e+03,Inf]" = "(1000+]")),

         GED_brk = plyr::revalue(GED_brk, c("(250,1e+03]" = "(250,1000]")),
         GED_brk = plyr::revalue(GED_brk, c("(1e+03,Inf]" = "(1000+]")) )

dat_pn_dist_yr <- dat_pn_dist_yr %>%
  mutate(SIGACT_monthly_sum_brk = cut(SIGACT_monthly_sum, breaks = pltbreaks, include.lowest = TRUE),
         ICEWS_monthly_sum_brk  = cut(ICEWS_monthly_sum,  breaks = pltbreaks, include.lowest = TRUE),
         GED_monthly_sum_brk    = cut(GED_monthly_sum,    breaks = pltbreaks, include.lowest = TRUE)) %>%
  mutate(SIGACT_monthly_sum_brk = plyr::revalue(SIGACT_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         SIGACT_monthly_sum_brk = plyr::revalue(SIGACT_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")),

         ICEWS_monthly_sum_brk = plyr::revalue(ICEWS_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         ICEWS_monthly_sum_brk = plyr::revalue(ICEWS_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")),

         GED_monthly_sum_brk = plyr::revalue(GED_monthly_sum_brk, c("(250,1e+03]" = "(250,1000]")),
         GED_monthly_sum_brk = plyr::revalue(GED_monthly_sum_brk, c("(1e+03,Inf]" = "(1000+]")) )

rm(pltbreaks)







# Map formatting options
# Map colors - fixed
map_cols = viridis::viridis(n = 6)
map_lvls = levels(dat_pn_dist_yr$SIGACT_monthly_sum_brk)

# Baghdad - for inset box
bdad = st_bbox(dat_pn_dist_yr %>% filter(governorate == "Baghdad"))



# SIGACTS Faceted Maps
map_main <-
  ggplot(data = dat_pn_dist_yr) +
  geom_sf(aes(fill = SIGACT_monthly_sum_brk),
          lwd    = 0,
          colour = "white") +
  scale_fill_manual(NULL,
                    values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  labs(title = "SIGACT Reported Events") +
  guides(fill=guide_legend(nrow  = 2,
                           byrow = TRUE)) +
  geom_rect(aes(xmin = bdad[1],
                xmax = bdad[3],
                ymin = bdad[2],
                ymax = bdad[4]), fill = "NA", color = "black", size = 1) +
  theme_void() +
  theme(strip.text.x = element_text( margin = margin( b = 10, t = 0) ),
        strip.text   = element_text(face = "bold", size = 12),
        plot.margin = unit(rep(1, 4), unit = "cm"),
        legend.justification = "center",
        legend.position = "bottom") +
  facet_wrap( ~ year)

map_inset <-
  ggplot(data = dat_pn_dist_yr %>% filter(governorate == "Baghdad")) +
  geom_sf(aes(fill = SIGACT_monthly_sum_brk),
          lwd = 0,
          colour = "white") +
  scale_fill_manual(values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.position = "none") +
  facet_wrap( ~ year)


# Inset maps to list
pp <- map(unique(dat_pn_dist_yr$year), function(x) {
  map_inset$data <- map_inset$data %>% filter(year == x)
  map_inset +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none",
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
})

# Inset kist to dataframe with relevant ploting value for insets
inset <- tibble(x = c(rep(1.00, 6)),
                y = c(rep(1.00, 6)),
                plot = pp,
                year = unique(dat_pn_dist_yr$year))


# Print map
map_main +
  ggpmisc::geom_plot_npc(data = inset, aes(npcx = x, npcy = y, label = plot, vp.width = 0.25, vp.height = 0.25))







# ICEWS Faceted Maps
map_main <-
  ggplot(data = dat_pn_dist_yr) +
  geom_sf(aes(fill = ICEWS_monthly_sum_brk),
          lwd    = 0,
          colour = "white") +
  scale_fill_manual(NULL,
                    values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  labs(title = "ICEWS Reported Events") +
  guides(fill=guide_legend(nrow  = 2,
                           byrow = TRUE)) +
  geom_rect(aes(xmin = bdad[1],
                xmax = bdad[3],
                ymin = bdad[2],
                ymax = bdad[4]), fill = "NA", color = "black", size = 1) +
  theme_void() +
  theme(strip.text.x = element_text( margin = margin( b = 10, t = 0) ),
        strip.text   = element_text(face = "bold", size = 12),
        plot.margin = unit(rep(1, 4), unit = "cm"),
        legend.justification = "center",
        legend.position = "bottom") +
  facet_wrap( ~ year)

map_inset <-
  ggplot(data = dat_pn_dist_yr %>% filter(governorate == "Baghdad")) +
  geom_sf(aes(fill = ICEWS_monthly_sum_brk),
          lwd = 0,
          colour = "white") +
  scale_fill_manual(values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.position = "none") +
  facet_wrap( ~ year)


# Inset maps to list
pp <- map(unique(dat_pn_dist_yr$year), function(x) {
  map_inset$data <- map_inset$data %>% filter(year == x)
  map_inset +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none",
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
})

# Inset kist to dataframe with relevant ploting value for insets
inset <- tibble(x = c(rep(1.00, 6)),
                y = c(rep(1.00, 6)),
                plot = pp,
                year = unique(dat_pn_dist_yr$year))


# Print map
map_main +
  ggpmisc::geom_plot_npc(data = inset, aes(npcx = x, npcy = y, label = plot, vp.width = 0.25, vp.height = 0.25))





# GED Faceted Maps
map_main <-
  ggplot(data = dat_pn_dist_yr) +
  geom_sf(aes(fill = GED_monthly_sum_brk),
          lwd    = 0,
          colour = "white") +
  scale_fill_manual(NULL,
                    values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  labs(title = "GED Reported Events") +
  guides(fill=guide_legend(nrow  = 2,
                           byrow = TRUE)) +
  geom_rect(aes(xmin = bdad[1],
                xmax = bdad[3],
                ymin = bdad[2],
                ymax = bdad[4]), fill = "NA", color = "black", size = 1) +
  theme_void() +
  theme(strip.text.x = element_text( margin = margin( b = 10, t = 0) ),
        strip.text   = element_text(face = "bold", size = 12),
        plot.margin = unit(rep(1, 4), unit = "cm"),
        legend.justification = "center",
        legend.position = "bottom") +
  facet_wrap( ~ year)

map_inset <-
  ggplot(data = dat_pn_dist_yr %>% filter(governorate == "Baghdad")) +
  geom_sf(aes(fill = GED_monthly_sum_brk),
          lwd = 0,
          colour = "white") +
  scale_fill_manual(values = setNames(map_cols, map_lvls)) +
  coord_sf(expand = TRUE) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.position = "none") +
  facet_wrap( ~ year)


# Inset maps to list
pp <- map(unique(dat_pn_dist_yr$year), function(x) {
  map_inset$data <- map_inset$data %>% filter(year == x)
  map_inset +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none",
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
})

# Inset kist to dataframe with relevant ploting value for insets
inset <- tibble(x = c(rep(1.00, 6)),
                y = c(rep(1.00, 6)),
                plot = pp,
                year = unique(dat_pn_dist_yr$year))


# Print map
map_main +
  ggpmisc::geom_plot_npc(data = inset, aes(npcx = x, npcy = y, label = plot, vp.width = 0.25, vp.height = 0.25))





