#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          June 17, 2021
# Purpose:       This script is a first pass extension of Silverman (2020) to
#                account for spatio-temproal dynamics in sigacts data
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
library(spatialreg)
library(sf)
library(magic)
library(texreg)

#---------------------------#
# Set working directory
#---------------------------#
# setwd()

#---------------------------#
# Load data
#---------------------------#
load("scripts/silverman-replication/silverman_halfyr.Rdata")
irq <- read_sf("data/esoc_v3/gis boundary files/iraq_district_boundaries_utm.shp")
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# TIDY DATA                                                               ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Create spatial data frame
# ----------------------------------- #
# Tidy spatial data
irq %<>%
  rename(district = ADM3NAME) %>%
  select(district, AREA_KM2, PERIM_KM)

d <- left_join(silver, irq, by = 'district') %>%
  ungroup() %>%
  drop_na(p_S1_lag) %>%
  arrange(halfyr, district) %>%
  st_set_geometry("geometry")
# ----------------------------------- #


# ----------------------------------- #
# Create spatio-temporal weights
# ----------------------------------- #
cs   <- as_Spatial(d %>% filter(halfyr == 90))
nbs  <- spdep::poly2nb(cs, queen = TRUE, row.names = cs$district)

W_cs <- spdep::nb2mat(neighbours = nbs, style = "W")
colnames(W_cs) <- rownames(W_cs)

W_ts <- do.call(what = magic::adiag,
                args = replicate(n = length(unique(d$halfyr)),
                                 expr = W_cs, simplify = FALSE))

ids <- expand.grid(rownames(W_cs), min(d$halfyr):max(d$halfyr))
rownames(W_ts) <- apply(ids, 1, paste, collapse=".")
colnames(W_ts) <- rownames(W_ts)

nb_cs <- spdep::mat2listw(W_cs)
nb_ts <- spdep::mat2listw(W_ts)
# ----------------------------------- #


# ----------------------------------- #
# Cleanup
# ----------------------------------- #
rm(cs, ids, irq, silver)
# ----------------------------------- #

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# VISUALIZE EVENTS                                                        ----
#-----------------------------------------------------------------------------#
summary(d$p_S1_d)

d <- d %>%
  mutate(p_S1_d_breaks = cut(p_S1_d,
                             breaks = c(-Inf, -5, -1, 0, 1, 5, Inf),
                             include.lowest = T),
         time_label = case_when(halfyr == 90 ~ "2005-1",
                                halfyr == 91 ~ "2005-2",
                                halfyr == 92 ~ "2006-1",
                                halfyr == 93 ~ "2006-2",
                                halfyr == 94 ~ "2007-1",
                                halfyr == 95 ~ "2007-2",
                                halfyr == 96 ~ "2008-1",
                                halfyr == 97 ~ "2008-2")) %>%
  mutate(time_label = factor(time_label,
                             levels = c("2005-1","2005-2","2006-1","2006-2",
                                        "2007-1","2007-2","2008-1","2008-2")),

         p_S1_d_breaks = as.character(p_S1_d_breaks)) %>%
  mutate(p_S1_d_breaks = case_when(p_S1_d        == 0 ~ "No change",
                                   p_S1_d_breaks == "(0,1]" ~ "[0,1]",
                                   p_S1_d_breaks == "[-Inf,-5]" ~ "-5+",
                                   p_S1_d_breaks == "(5, Inf]" ~ "5+",
                                   TRUE ~ p_S1_d_breaks)) %>%
  mutate(p_S1_d_breaks = factor(p_S1_d_breaks,
                                levels = c("-5+","(-5,-1]","(-1,0]",
                                           "No change",
                                           "[0,1]","(1,5]","5+")))

ggplot(data = d) +
  geom_sf(aes(fill = p_S1_d_breaks), size = 0.01) +
  scale_fill_brewer(palette   = "RdBu",
                    drop      = FALSE,
                    direction = -1) +
  facet_wrap(~ time_label, nrow = 2, ncol = 4) +
  theme_minimal() +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title     = element_blank(),
        legend.key.size  = unit(c(2,2), "mm"),
        panel.grid       = element_blank(),
        panel.background = element_rect(fill  = "transparent",
                                        color = "black"),
        axis.text        = element_blank(),
        strip.background = element_rect(fill  = "gray90",
                                        color = "black")) +
  labs(title = "Change in SIGACT events from previous half-year") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(filename = "scripts/silverman-replication/silverman-writeup/sigact-map.png",
       width    = 08,
       height   = 05,
       units    = "in",
       dpi      = 500,
       scale    = 1.0)

# ggplot(data = d) +
#   geom_histogram(aes(x = p_S1_d, fill = time_label))
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# ANALYSIS                                                                ----
#-----------------------------------------------------------------------------#
f1 <- formula(p_S1_d ~ -1 + p_S1_d_lag + p_spentcptotal_d + p_spentruzicka_d +
                     coalitioncc_d + insurgentcc_d +
                     p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
                     a_of_batt_d + cmoc + dis_usprt +
                su_vh2 + su_vh3 + su_vh4 + su_vh5 + su_vh6 + su_vh7 +
                su_vh8 + su_vh9 + su_vh10 + district + as.factor(halfyr))

# f1 <- formula(p_S1 ~ p_S1_lag + p_spentcptotal + p_spentruzicka +
#                 coalitioncc + insurgentcc +
#                 p_spentcerpsmall_noncp + p_spentusaid_nonruzicka +
#                 a_of_batt + cmoc + dis_usprt)

m1ex <- lm(formula = f1, data = d)
# summary(m1es)
# lmtest::coeftest(m1es, vcov = sandwich::vcovCL(m1ex, cluster = d$district))

spdiag <- spdep::lm.LMtests(model = m1ex,
                            listw = nb_ts,
                            test  = c('LMerr','LMlag','RLMerr','RLMlag'))
summary(spdiag)

m1sp <- lagsarlm(formula = f1,
                 data    = d,
                 listw   = nb_ts)
# summary(m1sp)
# confint(m1sp)[1:4,]

# Because fuck me that's why
m1sp_tidy <- broom::tidy(lmtest::coeftest(m1sp))
m1sp_tidy <- createTexreg(coef.names = m1sp_tidy$term,
                          coef       = m1sp_tidy$estimate,
                          se         = m1sp_tidy$std.error,
                          pvalues    = m1sp_tidy$p.value,
                          gof.names  = c("Num. obs."),
                          gof        = c(round(nrow(d),0)))

screenreg(list(m1ex,m1sp), custom.coef.map = c(vars, rho = "Rho", p_S1_d_lag = "Phi"),
          include.aic = F,
          include.rsquared = F,
          include.adjrs = F,
          include.lr = F,
          include.log = F,
          include.nobs = F,
          custom.model.names = c("Model 4(d) + Temporal Lag", "Spatio-Temporal"))

confint(m1ex)[1:3,]
coef(m1ex)[1:3]
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# EFFECTS - LRSS                                                          ----
#-----------------------------------------------------------------------------#

# SpatialNLS::impacts

# ----------------------------------- #
# Setup - Point Estimates
# ----------------------------------- #
rho        <- coef(m1sp)['rho']
phi        <- coef(m1sp)['p_S1_d_lag']
names(phi) <- "phi"
variable   <- "p_spentcptotal_d"
# variable <- names(coef(m1sp))
# variable <- variable[!str_detect(str_to_lower(variable),
#                                  pattern = c("intercept|rho|lag|district|factor"))]
bs <- coef(m1sp)[variable]
n  <- nrow(W_cs)
Id <- diag(1, n, n)
M  <- solve(Id - rho * W_cs - phi * Id)

im <- sapply(variable, function(x) {
  b <- coef(m1sp)[x]
  bM <- b * M
  dir <- mean(diag(bM))
  ind <- (sum(bM) - sum(diag(bM)))/n
  tot <- sum(bM)/n
  return(cbind(Direct = dir, Indirect = ind, Total = tot))
}, simplify = FALSE)
im <- do.call(rbind, im) %>% as.data.frame()
rownames(im) <- variable
# ----------------------------------- #


# ----------------------------------- #
# LRSS - Simulation for inference
# ----------------------------------- #
m1sp_vcv <- vcov(m1sp)[c('rho','p_S1_d_lag',variable),
                       c('rho','p_S1_d_lag',variable)]

sims <- 1000

# Create array to store n*n matrices over 1000 sims
sim.effs <- array(data     = 0,
                  dim      = c(nrow(im),n,n,sims),
                  dimnames = list(rownames(im),
                                  d$district[1:n],
                                  d$district[1:n],
                                  sprintf("sim%s",1:sims)))
print(object.size(sim.effs), units = "Mb")

for(i in 1:sims){
  draws <- MASS::mvrnorm(n         = 1,
                         mu        = c(rho, phi, bs),
                         Sigma     = m1sp_vcv)

  sim.rho   <- draws['rho']
  # sim.beta = draws['wage95']
  sim.phi   <- draws['phi']
  sim.betas <- draws[!str_detect(str_to_lower(names(draws)),
                                pattern = c("rho|phi"))]
  sim.m     <- solve(Id - sim.rho * W_cs - sim.phi * Id)

  for(name in names(sim.betas)){
    sim.b  <- sim.betas[[name]]
    sim.bM <- sim.b * sim.m
    sim.effs[name,,,i] <- sim.bM
  }
}
rm(draws, sim.rho, sim.phi, sim.betas, sim.m, sim.b, sim.bM)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# LRSS - TIDY                                                             ----
#-----------------------------------------------------------------------------#

# EVALUATE SIMULATION RESULTS

# ----------------------------------- #
# VALUE OF SPATIO-TEMPORAL EXAMPLES
# ----------------------------------- #
# # Median indirect effect of unit-change in 'p_spentcptotal_d' in Najaf on 'p_S1_d' in Al-Salman
# names(sim.effs[,1,1,1])
# median(sim.effs["p_spentcptotal_d","Najaf","Al-Salman",])
# median(sim.effs["p_spentcptotal_d","Najaf","Zakho",])
# quantile(sim.effs["p_spentcptotal_d","Najaf","Al-Salman",], probs = c(0.025, 0.975))
# ----------------------------------- #


# ----------------------------------- #
# AVERAGE DIRECT EFFECT (ADE)
# ----------------------------------- #
ade.mean <- sapply(variable, function(x){
  mean(apply(sim.effs[x,,,], 3, function(y) median(diag(y))))
}) %>%
  as.data.frame() %>%
  rename(Mean = 1)

# 95% credibility interval
ade.ci   <- sapply(variable, function(x){
  quantile(apply(sim.effs[x,,,], 3, function(y) median(diag(y))),
           probs = c(0.025, 0.975))
}) %>%
  t %>% as.data.frame %>%
  rename(LB = 1,
         UB = 2)

ade <- bind_cols(ade.mean, ade.ci)
rownames(ade) <- "Direct"
# ----------------------------------- #


# ----------------------------------- #
# AVERAGE INDIRECT EFFECT (AIE)
# ----------------------------------- #
aie.mean <- sapply(variable, function(x){
  mean(apply(sim.effs[x,,,], 3, function(y) (sum(y) - sum(diag(y))) / n ))
}) %>%
  as.data.frame() %>%
  rename(Mean = 1)

# 95% credibility interval
aie.ci   <- sapply(variable, function(x){
  quantile(apply(sim.effs[x,,,], 3, function(y) (sum(y) - sum(diag(y))) / n  ),
           probs = c(0.025, 0.975))
}) %>%
  t %>% as.data.frame %>%
  rename(LB = 1,
         UB = 2)

aie <- bind_cols(aie.mean, aie.ci)
rownames(aie) <- "Inirect"
# ----------------------------------- #


# ----------------------------------- #
# AVERAGE TOTAL EFFECT (ATE)
# ----------------------------------- #
ate.mean <- sapply(variable, function(x){
  mean(apply(sim.effs[x,,,], 3, function(y) sum(y)  / n ))
}) %>%
  as.data.frame() %>%
  rename(Mean = 1)

# 95% credibility interval
ate.ci   <- sapply(variable, function(x){
  quantile(apply(sim.effs[x,,,], 3, function(y) sum(y) / n  ),
           probs = c(0.025, 0.975))
}) %>%
  t %>% as.data.frame %>%
  rename(LB = 1,
         UB = 2)

ate <- bind_cols(ate.mean, ate.ci)
rownames(ate) <- "Total"
# ----------------------------------- #


# ----------------------------------- #
# LRSS
# ----------------------------------- #
lrss <- bind_rows(ade, aie, ate) %>% select(LB, Mean, UB)
lrss <- lrss %>%
  rownames_to_column(var = "Type") %>%
  mutate(to_table = sprintf("%s [%s, %s]",
                            format(round(Mean, 3), nsmall = 3),
                            format(round(LB, 3),   nsmall = 3),
                            format(round(UB, 3),   nsmall = 3)))
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# EFFECTS - MARGINAL RESPONSE PATHS                                       ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# NOTE ----- NEED TO IMPROCE EFFICIENCY HERE.
# Presently 692 million calculations (832^2 * 1e3) / 1e6
# ----------------------------------- #

t    <- length(unique(d$halfyr))                # Number of time units
nobs <- n*t               # Now the number of observation is N x T, not N
Id   <- diag(1,nobs,nobs)

# Create the Time-Lag Generating Matrix, "L" with dimension (nt * nt):
L <- diag(nobs-n)                 # nt-n x nt-n matrix with diagonal = 1, 0 elsewhere
L <- rbind(matrix(0,n,nobs-n),L)  # Add n rows to top of 0s to this [first year lag]
L <- cbind(L,matrix(0,nobs,n))    # Add n columns to right of 0s.

# Estimate LRSS for Spatio-Temporal Data: -- NOTICE: "L" rather than "Id"
M  <- solve(Id - rho*W_ts - phi*L)

# Estimate response (marginal responses)
bM <- bs * M

bM[1:4, 1:4]
bM[476:480, 476:480]

# Tarmia experienced the greatest single half year incease in significant events
d[which.max(d$p_S1_d),]$district

report.neighbors <- function(w, unit){
  neighbors = rownames(subset(w,w[,unit]>0))
  return(neighbors)
}

report.neighbors(w = W_cs, unit = "Tarmia")

# Estimated (point) effect of a $100 increase in per capita condolence spending
# in Falluja in the first-half of 2005 on reported violence in Tarmia:
eff <- bM[,"Falluja.90"][str_detect(names(bM[,"Falluja.90"]), "Tarmia")]

eff <- data.frame(y    = cumsum(eff),
                  x    = 1:length(eff),
                  xlab = unique(d$time_label))

eff_plt <- ggplot(data = eff, aes(x = x, y = y)) +
  geom_line() +
  scale_x_continuous(breaks = 1:8,
                     labels = eff$xlab) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black", size = 0.2),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray80", size = 0.2)) +
  labs(title = "Response in SIGACTS per capita in Tarmia",
       subtitle = "to a $1 per capita increase in condolence spending in Falluja in 2005")


eff_mp <- d %>%
  filter(halfyr == 90) %>%
  select(district) %>%
  mutate(lab = case_when(district %in% c("Tarmia", "Falluja") ~ district))

eff_mp <- ggplot(data = eff_mp) + geom_sf(aes(fill = lab)) +
  scale_fill_discrete(na.translate = F) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        # panel.background = element_rect(fill = "white", color = "black", size = 0.2),
        plot.background = element_rect(fill = "white", color = "black", size = 0.2),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

library(cowplot)
eff <- ggdraw() +
  draw_plot(eff_plt) +
  draw_plot(eff_mp,
            x = 0.495, y = 0.24,
            hjust = 0, vjust = 0,
            valign = 0, halign = 0,
            width = 0.65, height = 0.65)

ggsave(plot = eff,
       filename = "scripts/silverman-replication/silverman-writeup/marginal-response-ex.png",
       width    = 08,
       height   = 05,
       units    = "in",
       dpi      = 500,
       scale    = 1.0)
# ----------------------------------- #
# Marginal Response - simulation for inference
# ----------------------------------- #
# sims  <- 10                  # Number of simulations to run
#
# sim.effs <- array(data     = 0,
#                   dim      = c(nobs,nobs,sims),
#                   dimnames = list(rownames(bM),
#                                   rownames(bM)))
#
# for(i in 1:sims){
#   # if(i %in% seq(100, sims, by = 100)){
#   #   print(sprintf('Working on iteration: %s', i))
#   # }
#
#   draws <- MASS::mvrnorm(n         = 1,
#                          mu        = c(rho, phi, bs),
#                          Sigma     = m1sp_vcv)
#
#   sim.rho  <- draws['rho']
#   sim.phi  <- draws['phi']
#   sim.beta <- draws[!str_detect(str_to_lower(names(draws)),
#                                 pattern = c("rho|phi"))]
#
#   sim.m  = solve(Id - sim.rho * W_ts - sim.phi * L)
#   sim.bm = sim.beta * sim.m
#
#   rownames(sim.bm) = colnames(sim.bm) =  rownames(bM)
#
#   sim.effs[,,i] = sim.bm
# }
# # Clean-up, remove unecessary objects left-over from simulation
# rm(draws, sim.rho, sim.beta, sim.m, sim.bm)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# MARGINAL RESPONSE - TIDY                                                ----
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
save(list = c("m1ex","spdiag","m1sp","m1sp_tidy","lrss"),
     file = "scripts/silverman-replication/silverman-writeup/spatial-extension.Rdata")
rm(list = ls())
#-----------------------------------------------------------------------------#
