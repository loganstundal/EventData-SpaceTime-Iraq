#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          September 05, 2021
# Purpose:       descriptives - APSA paper figures and descriptive statistics
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
library(tidyverse)
library(sf)
library(kableExtra)
library(skimr)
library(biscale)
library(cowplot)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("data/Oxford_HB_2021_APSA-Data.Rdata")
#---------------------------#

#---------------------------#
# Locals
#---------------------------#
map_theme <- theme(panel.background = element_rect(fill = "transparent", color = "black",
                                                   size = 0.2),
                   panel.grid       = element_blank(),
                   axis.text        = element_blank(),
                   axis.ticks       = element_blank(),
                   legend.title     = element_blank(),
                   legend.position  = "bottom",
                   legend.direction = "horizontal",
                   strip.background = element_rect(fill = "gray90", color = "black",
                                                   size = 0.2))
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# DESCRIPTIVE STATISTICS                                                  ----
#-----------------------------------------------------------------------------#
d <- irq_halfyr %>% drop_na(p_s1_d_lag) %>%
  select(p_s1_d, p_spentcptotal_d, p_spentruzicka_d, coalitioncc_d,
         insurgentcc_d, p_spentcerpsmall_noncp_d, p_spentusaid_nonruzicka_d,
         a_of_batt_d, cmoc, dis_usprt)

# skimr::skim(d)

x <- list("p_s1_d"                    = "Insurgent Events",
          "p_spentcptotal_d"          = "Condolence Spending (PC)",
          "p_spentruzicka_d"          = "Ruzicka Spending (PC)",
          "coalitioncc_d"             = "Coalition Collateral Damage",
          "insurgentcc_d"             = "Insurgent Collateral Damage",
          "p_spentcerpsmall_noncp_d"  = "Other Small CERP Spending",
          "p_spentusaid_nonruzicka_d" = "Other USAID Spending",
          "a_of_batt_d"               = "Coalition Troop Strength",
          "cmoc"                      = "CMOC Presence",
          "dis_usprt"                 = "PRT Presence")

descr <- d %>% summarize(
  across(.vars = everything(),
         .fns  = ~c("Mean" = mean(.x),
                    "SD"   = sd(.x),
                    "Min." = min(.x),
                    "Max." = max(.x)))) %>% t %>% as.data.frame %>%
  rename("Mean" = 1,
         "SD"   = 2,
         "Min." = 3,
         "Max." = 4) %>%
  mutate(across(everything(), ~round(.x, 3))) %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = recode(Variable, !!!x))

# Generate latex code for descriptive statistics
kbl(x = descr, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = T) %>%
  column_spec(., 1, width = "1.5in")
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# AVERAGE DISTRICT BORDER DISTANCE                                        ----
#-----------------------------------------------------------------------------#

# NOTE - use this estimate to motivate GMRF model and range which is greater
# than average border distance implying non-localized conflict processes.

# st_centroid(irq0)
# irq1_prj <- st_transform(irq1,
#                          "+proj=aea +lat_1=29 +lat_2=37 +lon_0=33.05 +datum=NAD83 +units=m +no_defs")

x      <- irq$district_name
res_km <- lapply(x, function(id){
  tmp <- irq %>% filter(district_name == id) %>% select(district_name)
  cnt <- suppressWarnings({st_centroid(tmp)})

  tmp_bdr <- tmp %>% st_union() %>% st_cast(., "POINT")
  res <- st_distance(tmp_bdr, cnt) %>% mean %>% as.numeric
  res <- res / 1e3
  return(res)
}) %>% unlist
summary(res_km)
quantile(res_km, probs = c(0.025, 0.5, 0.975))

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# MAPS                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Tidy half-year id for plot:
# ----------------------------------- #
irq_halfyr <- irq_halfyr %>%
  mutate(time_id = str_replace(time_id, "h", " - Half "))
# ----------------------------------- #


# ----------------------------------- #
# Map data
# ----------------------------------- #
mp_dat <- irq_halfyr %>%
  select(district_name, time_id, p_s1_d, p_spentcptotal_d, p_spentruzicka_d) %>%
  left_join(., irq[,c("district_name","geometry")], by = "district_name") %>%
  st_set_geometry(., "geometry")
# ----------------------------------- #


# ----------------------------------- #
# Insurgent Events
# ----------------------------------- #
insurgent_attacks <- ggplot(data = mp_dat) +
  geom_sf(aes(fill = p_s1_d), size = 0.2) +
  scale_fill_gradient2(low      = scales::muted("green"),
                       high     = scales::muted("red"),
                       mid      = "white",
                       midpoint = 0,
                       limits   = c(-15, 10)) +
  facet_wrap(~time_id) +
  map_theme +
  labs(
    # title    = "Insurgent Attacks",
    # subtitle = "SIGACTS-III Data",
    caption  = "Map produced with data from Silverman (2021) and SIGACTS-III database.")

# Save map
ggsave(plot     = insurgent_attacks,
       filename = "Results/Figures/figure-insurgent_attacks.png",
       width    = 6.5,
       height   = 7.0,
       units    = "in",
       dpi      = 350)

# file.copy(
#   from = "Results/Figures/figure-insurgent_attacks.png",
#   to   = "Drafts/Drafts/GAST20210906/figure-insurgent_attacks.png",
#   overwrite = TRUE)
# ----------------------------------- #


# ----------------------------------- #
# Spending
# ----------------------------------- #
      # d <- irq_halfyr %>%
      #   select(district_name, p_spentcptotal_d, p_spentruzicka_d) %>%
      #   group_by(district_name) %>%
      #   summarize(across(everything(), mean)) %>%
      #   ungroup %>%
      #   # mutate(p_spentcptotal_d = p_spentcptotal_d - mean(p_spentcptotal_d),
      #   #        p_spentruzicka_d = p_spentruzicka_d - mean(p_spentruzicka_d)) %>%
      #   left_join(., irq[,c("district_name","geometry")], by = "district_name") %>%
      #   st_set_geometry(., "geometry")
      #
      # x <- 0.01
      # table(abs(irq_halfyr$p_spentcptotal_d) > x) %>% prop.table() %>% round(.,2)
      # table(abs(irq_halfyr$p_spentruzicka_d) > x) %>% prop.table() %>% round(.,2)
      #
      #
      # summary(d$p_spentcptotal_d) %>% round(.,3)

# Using cross-section average:
mp_dat2 <- mp_dat %>%
  st_drop_geometry %>%
  select(-time_id) %>%
  group_by(district_name) %>%
  summarize(across(everything(), median)) %>%
  ungroup %>%
  left_join(., irq[,c("district_name", "geometry")], by = "district_name") %>%
  st_set_geometry(., "geometry")

spend_cp <- ggplot(data = mp_dat2) +
  geom_sf(aes(fill = p_spentcptotal_d), size = 0.2) +
  scale_fill_gradient2(low      = scales::muted("red"),
                       high     = scales::muted("green"),
                       mid      = "white",
                       midpoint = 0,
                       limits   = c(-0.1, 0.06)) +
  map_theme +
  theme(title = element_text(size = 12/.pt),
        text  = element_text(size = 10/.pt),
        legend.key.size = unit(5, "mm")) +
  labs(title    = "Condolence Spending",
       subtitle = "Median per capita spending between 2004.H2 - 2008.H2")

ggsave(plot     = spend_cp,
       filename = "Results/Figures/figure-spending_cp.png",
       width    = 4.0,
       height   = 4.0,
       units    = "in",
       dpi      = 350)


table(irq_halfyr$p_spentruzicka_d == 0.00) %>% prop.table() %>% round(. , 3)
table(abs(irq_halfyr$p_spentruzicka_d) > 0.10) %>% prop.table() %>% round(. , 3)

lapply(unique(irq_halfyr$time_id),
       function(x){
         irq_halfyr %>%
           filter(time_id == x) %>%
           pull(p_spentruzicka_d) %>% table(. == 0)})

# ggplot(data = mp_dat2) +
#   geom_sf(aes(fill = p_spentruzicka_d), size = 0.2) +
#   scale_fill_gradient2(low = scales::muted("red"),
#                        high = scales::muted("green"),
#                        mid = "white",
#                        midpoint = 0) +
#   map_theme
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# BIVARIATE MAP                                                           ----
#-----------------------------------------------------------------------------#
dat <- irq_halfyr

dat %>% select(time_id, p_s1_d) %>%
  group_by(time_id) %>%
  summarize(across(.fns = c("mean" = mean, "sd" = sd, "min" = min, "max" = max, "iqr" = IQR),
                   .names = "{.fn}")) %>%
  mutate("TEST" = iqr < max)

dat %>% select(time_id, p_spentcptotal_d) %>%
  group_by(time_id) %>%
  summarize(across(.fns = c("mean" = mean, "sd" = sd, "min" = min, "max" = max, "iqr" = IQR),
                   .names = "{.fn}")) %>%
  mutate("TEST" = iqr < max)

vals <- function(data, variable){
  data <- data %>%
    select(sym(variable)) %>%
    summarize(across(.fns = c("mean" = mean, "sd" = sd, "min" = min, "max" = max, "iqr" = IQR),
                    .names = "{.fn}"))
  return(data)
}

# vals(data = dat %>% filter(time_id == "2004h2"), "p_s1_d")
# vals(data = dat %>% filter(time_id == "2004h2"), "p_spentcptotal_d")

ids <- as.character(unique(dat$time_id))
chk <- lapply(ids, function(id){
  tmp <- dat %>% filter(time_id == id) %>% select(district_name, p_s1_d, p_spentcptotal_d) %>%
    rename(dv = 2, iv = 3)

  vals_dv <- vals(data = tmp, variable = "dv")
  vals_iv <- vals(data = tmp, variable = "iv")

  tmp <- tmp %>%
    mutate(cut_dv = case_when(dv <= 0 ~ "1",
                              dv >  0 & dv < vals_dv$iqr ~ "2",
                              dv >= vals_dv$iqr ~ "3"),
           cut_iv = case_when(iv <= 0 ~ "1",
                              iv >  0 & iv < vals_iv$iqr ~ "2",
                              iv >= vals_iv$iqr ~ "3")) %>%
    mutate(bi_class = paste(cut_iv, cut_dv, sep = "-"),
           time_id  = id)

  return(tmp)
})

d <- bind_rows(chk)
# head(d)
# table(d$bi_class)

d <- left_join(d, irq[,c("district_name","geometry")], by = "district_name") %>% st_set_geometry(., "geometry")

pal <- bi_pal_manual(val_1_3 = "#cc0024", val_2_3 = "#8a274a", val_3_3 = "#4b264d",
                     val_1_2 = "#dd7c8a", val_2_2 = "#8d6c8f", val_3_2 = "#4a4779",
                     val_1_1 = "#dddddd", val_2_1 = "#7bb3d1", val_3_1 = "#016eae",
                     # preview = TRUE
                     )

leg <- bi_legend(pal = pal,
                 dim = 3,
                 xlab = "Condolence Spending ",
                 ylab = "Insurgent Violence ",
                 size = 8)

mp <- ggplot(data = d) +
  geom_sf(aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = pal, dim = 3) +
  facet_wrap(~time_id) +
  map_theme +
  theme(plot.margin=unit(c(0, 0.1, 0, 0),"in"))

leg2   <- plot_grid(NULL, leg,  NULL, rel_heights = c(6.5/3, 1, 0.05), nrow = 3)
final  <- plot_grid(mp,   leg2, rel_widths = c(6.5/3, 1), ncol = 2, align = "v", axis = "b")

windows(width = 6.5, height = 5.5);final
dev.off(2)

ggsave(plot     = final,
       filename = "Results/Figures/figure-bimap.png",
       width    = 6.5,
       height   = 5.5,
       units    = "in")

# ggsave(plot     = leg,
#        filename = "Results/Figures/figure-bimap-leg.png",
#        width    = 6.0,
#        height   = 6.0,
#        units    = "in")

file.copy(
  from = "Results/Figures/figure-bimap.png",
  to   = "Drafts/Drafts/GAST20210913/figure-bimap.png",
  overwrite = TRUE
)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
#save.image()
#rm(list = ls())
#-----------------------------------------------------------------------------#
