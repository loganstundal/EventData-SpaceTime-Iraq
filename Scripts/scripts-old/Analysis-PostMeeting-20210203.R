#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          February 03, 2021
# Purpose:       Post-meeting analysis actions
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
        "To do list:
          1. GENERATE NEW VARIABLES
              - BINARY ~ EVENT / NO EVENT

          2. NEW DESCRIPTIVES
              - # DISTRICTS WITH NO TEMPORAL VARIATION
              - DO COUNTS RESOLVE THIS
              - DOES ADDITIONAL AGGREGATION (i.e., HALF-YEAR?)

          3. GENERATE CONFUSION MATRIX (SUCH AS IN COLOMBIA PAPER)
              - INCLUDE KAPPA MEASURES, SIGACTS [BINARY EVENT/NO EVENT]
                AS TRUTH

          4. PROBIT MODELS W/O POPULATION CORRECTION FOR ALL 3 DVS

          5. Estimate STAR w/ w/o fixed effects in LEVELS AND DIFFS
              - PERFORM IN HALF-YEARS AND MONTHLY AGGREGATIONS

          FOR MORE INFO., SEE Meeting-20210203.docx IN PROJECT DIRECTORY."


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
library(sf)
library(lubridate)

#---------------------------#
# Set working directory
#---------------------------#
# setwd()

#---------------------------#
# Load data
#---------------------------#
load("../Data/Data-All-Tidy-20210131.Rdata")
load("../Data/Data-Spatial-Weights-20210131.Rdata")

# rm(W, W_list, W_eigen, d_ts, d_cs)


# ----------------------------------- #
# Whatever this function I wrote does...
# ----------------------------------- #
kappa_cm <- function(observed,
                     pred   = NULL,
                     p_bin  = NULL,
                     cutoff = 0.5,
                     sig    = 0.95,
                     silent = TRUE){
  if(is.null(pred) & is.null(p_bin)){
    stop('Either "pred" or "p_bin" are required.')
  }

  if(is.null(p_bin)){
    p_bin = ifelse(pred >= cutoff, 1, 0)
  }

  tab = table(observed, p_bin)

  p_agree  = (tab[1] + tab[4]) / sum(tab)
  p_random = (((tab[1] + tab[3]) / sum(tab)) * ((tab[1] + tab[2]) / sum(tab))) +
    (((tab[2] + tab[4]) / sum(tab)) * ((tab[3] + tab[4]) / sum(tab)))

  k    = (p_agree - p_random) / (1 - p_random)
  se_k = sqrt(( p_agree * (1 - p_agree) ) / ( length(observed) * (1 - p_random)^2  ))

  crit = qnorm(1 - (1-sig)/2)

  k_lci = k - crit*se_k
  k_uci = k + crit*se_k

  return(list('Kappa' = k,
              'Kappa_lci' = k_lci,
              'Kappa_uci' = k_uci))

  do.call(sprintf, c(list('Kappa: %s\n
                Kappa CI: [%s, %s]')),
          paste(round(c(k, k_lci, k_uci),3)))

  if(silent == FALSE){
    cat(sprintf('Kappa: %s\nKappa CI: [%s, %s]', round(k,3),round(k_lci,3),round(k_uci,3)))
  }
}
# ----------------------------------- #



#-----------------------------------------------------------------------------#
# NEW VARIABLES | REDUCE MAP GEOMETRY -----------------------------------------

"NOTE: integrate this component into 'Scrips/Tidy-Data-YYYYMMDD.R' on next full-data compile."
m <- d_pn_mo %>%
  mutate(# BINARY EVENTS
         SIGACT_BINARY  = case_when(SIGACT_ct > 0 ~ 1, TRUE ~ 0),
         ICEWS_BINARY   = case_when(ICEWS_ct  > 0 ~ 1, TRUE ~ 0),
         GED_BINARY     = case_when(GED_ct    > 0 ~ 1, TRUE ~ 0),

         # UNDER-REPORTING
         ICEWS_UNDER = case_when(ICEWS_BINARY == 0 & SIGACT_BINARY == 1 ~ 1, TRUE ~ 0),
         GED_UNDER   = case_when(GED_BINARY   == 0 & SIGACT_BINARY == 1 ~ 1, TRUE ~ 0),

         YRMON = ymd(paste(Year, Month, 01, sep = "-")))

# with(m, table(SIGACT_BINARY, ICEWS_BINARY, GED_BINARY))


# MIGRATE THIS TO TIDY SCRIPT AT NEXT DATA UPDATE
print(object.size(m), units = "Mb")
m <- rmapshaper::ms_simplify(m, keep = 0.01)
print(object.size(m), units = "Mb")


#-----------------------------------------------------------------------------#
# UNDERREPORTING MAPS ---------------------------------------------------------
under_dvs <- list('ICEWS' = 'ICEWS_UNDER','GED' = 'GED_UNDER')


mplt <- m %>%
  filter(Year %in% c(2006)) %>%
  select(Year, Month, ICEWS_UNDER, GED_UNDER) %>%
  pivot_longer(cols = ends_with("_UNDER"),
               names_to = "VAR",
               values_to = "Bias_Under") %>%
  mutate(Bias_Under = as_factor(Bias_Under)) %>%
  st_set_geometry(., "geometry")


ggplot(data = mplt) +
  geom_sf(aes(fill = Bias_Under)) +
  facet_wrap(vars(Year,Month, VAR), ncol = 6) +
  theme_minimal()






# latitude of Baghdad, Iraq is 33.312805, and the longitude is 44.361488
{baghdad <- data.frame('x' = 44.361488,
                       'y' = 33.312805)
  baghdad <- st_as_sf(baghdad, coords = c('x','y'), crs = '+proj=longlat +datum=WGS84')
  baghdad <- st_transform(baghdad, crs = st_crs(m))}

cols = c("gray90",viridis::magma(n = 8)[2])

under_maps <- lapply(1:length(under_dvs), function(x){
  var = under_dvs[[x]]

  ggplot(data = m %>% mutate(tmp = as_factor(m[[var]]))) +
    geom_sf(aes(fill = tmp), size = 0.01, alpha = 0.9) +
    geom_sf(data = m %>% filter(Governorate == "Baghdad"), aes(color = "Baghdad"), size = 0.1, fill = NA) +
    # geom_sf(data = baghdad, aes(shape = 'Baghdad'), size = 3, fill = cols[3],
    #         show.legend = 'point') +
    scale_fill_manual(name   = '',
                      values = cols,
                      guide  = guide_legend(override.aes = list(linetype = "blank", shape = NA)),
                      labels = c('Agreement','Underreporting')) +
    scale_color_manual(name   = "",
                         values = "red") +
    # guides(color = guide_legend(title = NULL)) +
    theme_void() +
    facet_wrap(~ YRMON) +
    labs(title = sprintf('%s: Underreporting', names(under_dvs)[x])
         # caption = sprintf('Plot date: %s.', format(Sys.Date(), "%B %d, %Y"))
    )
})

under_maps[[1]]
under_maps[[2]]

cowplot::plot_grid(under_maps[[1]], under_maps[[2]], ncol = 2)

# ----------------------------------- #











#-----------------------------------------------------------------------------#
# CONFUSION MATRICES ----------------------------------------------------------
# Kappa statistics
d <- sapply(c('ICEWS_BINARY','GED_BINARY'), function(x){
  kappa_cm(observed = m$SIGACT_BINARY,
           p_bin    = m[[x]])
}) %>% t %>% apply(., 2, as.numeric) %>% as.data.frame


kappa_plot <- {ggplot(data = d, aes(x = c(0,1))) +
    geom_point(aes(y = Kappa, colour = c('ICEWS','GED'))) +
    geom_errorbar(aes(ymin = Kappa_lci,
                      ymax = Kappa_uci,
                      colour = c('ICEWS','GED')),
                  width = 0.10) +
    geom_hline(aes(yintercept = 0), linetype = 'solid', size = 1, color = 'gray70') +
    scale_x_continuous(name   = '',
                       limits = c(-0.25, 0.25)) +
    scale_y_continuous(name   = '',
                       limits = c(-0.01, 0.5)) +
    # scale_color_manual(values = c('ICEWS' = cols[1],
    #                               'GED'   = cols[2])) +
    # coord_equal(expand = FALSE) +
    theme_minimal() +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title       = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.text.x        = element_blank(),
          axis.title.x       = element_blank(),
          legend.position    = 'bottom',
          legend.direction   = 'horizontal') +
    labs(title    = 'Kappa values and 95% confidence intervals',
         subtitle = 'ICEWS and GED reliabilty compared to CINEP',
         caption  = sprintf('Plot date: %s.', format(Sys.Date(), "%B %d, %Y")))
}
# ggsave(filename = 'paper drafts/memo_20200325/plots_20200325/RawData_Kappa.png',
#        plot     = kappa_plot,
#        width    = 5,
#        height   = 7,
#        dpi      = 320)

# ----------------------------------- #

# Confusion Matrices
lapply(c('ICEWS_BINARY','GED_BINARY'), function(x){
  caret::confusionMatrix(as_factor(m$SIGACT_BINARY), as_factor(m[[x]]))
})



#-----------------------------------------------------------------------------#
# DATA SPARSITY ---------------------------------------------------------------
m <- st_drop_geometry(m)

m2 <- m %>% group_by(District) %>%
  mutate(sig_tst = case_when(SIGACT_BINARY != lag(SIGACT_BINARY, 1) ~ 1, TRUE ~ 0),
         ice_tst = case_when(ICEWS_BINARY  != lag(ICEWS_BINARY, 1)  ~ 1, TRUE ~ 0),
         ged_tst = case_when(GED_BINARY    != lag(GED_BINARY, 1)    ~ 1, TRUE ~ 0)) %>%
  ungroup() %>%
  select(sig_tst,
         ice_tst,
         ged_tst) %>%
  lapply(X = ., FUN = function(x){round(prop.table(table(x)),2)});m2
#-----------------------------------------------------------------------------#











#-----------------------------------------------------------------------------#
# SAVE ------------------------------------------------------------------------
# save.image()
# rm(list = ls())
