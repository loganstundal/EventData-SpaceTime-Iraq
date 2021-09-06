#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 14, 2021
# Purpose:       All APSA Models
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
library(INLA)
library(dplyr)
library(forcats)
library(tidyr)
library(sf)
library(ggplot2)
library(sandwich)
library(texreg)
library(spatialreg)
library(scales)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("data/Oxford_HB_2021_APSA-Data.Rdata")
source("scripts/fn_inla-table.R")
source("Scripts/fn_inla-field.R")
rm(inla_table)
#---------------------------#

#---------------------------#
# Local functions
#---------------------------#
l_aic <- function(llik, k){-2 * llik + 2 * k}
l_bic <- function(llik, k, n){-2 * llik + log(n) * k}

tidy_mod <- function(model, prec = 2){
  # Combine parameters for inla mods
  if(class(model) == "inla"){
    inla_res <- bind_rows(model$summary.fixed, model$summary.hyperpar)
  }

  # Get point estimates:
  if(class(model) %in% c("lm", "sarlm")){
    m_cfs <- coef(model)
    m_cfs <- m_cfs[!is.na(m_cfs)]
  } else if(class(model) == "inla"){
    m_cfs <- inla_res[,"0.5quant"]
    names(m_cfs) <- rownames(inla_res)
  }


  # Get HPDs:
  if(class(model) == "lm"){
    m_vcv <- vcovCL(x = model, cluster = ~ district_name, type = "HC1")
    m_ci  <- lmtest::coefci(model, vcov. = m_vcv) %>%
      as.data.frame %>%
      rename(lb = 1, ub = 2)
  } else if(class(model) == "sarlm"){
    m_ci <- confint(model) %>%
      as.data.frame %>%
      rename(lb = 1, ub = 2)
  } else if (class(model) == "inla"){
    m_ci <- inla_res[,c("0.025quant","0.975quant")] %>%
      rename(lb = 1, ub = 2)
  }


  # Get ancillary parameters:
  if(class(model) %in% c("lm","sarlm")){
    n   <- ifelse(class(model) == "lm",
                  nrow(model.matrix(model)),
                  nrow(model[["X"]]))
    llk <- logLik(model) %>% as.numeric %>% round(., prec)
  } else if(class(model) == "inla"){
    n   <- model$.args$data[[1]] %>% length
    llk <- model$mlik[2]
  }

  k   <- length(m_cfs)
  aic <- l_aic(llik = llk, k = k) %>% round(., prec)
  bic <- l_bic(llik = llk, k = k, n = n) %>% round(., prec)



  # Return object
  res <- createTexreg(
    coef.names = names(m_cfs),
    coef       = m_cfs,
    ci.low     = m_ci$lb,
    ci.up      = m_ci$ub,
    gof.names  = c("n", "LogLik", "AIC", "BIC"),
    gof        = c(n, llk, aic, bic),
    gof.decimal= c(F, T, T, T))

  return(res)
}

map_theme <-   theme(
  panel.background = element_rect(fill  = "transparent",
                                  color = "black",
                                  size  = 0.1),
  panel.grid       = element_blank(),
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.title       = element_blank(),
  legend.position  = "bottom",
  legend.direction = "horizontal",
  legend.title     = element_blank(),
  strip.background = element_rect(fill  = "gray80",
                                  color = "black",
                                  size  = 0.1))
#---------------------------#
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SETUP                                                                   ----
#-----------------------------------------------------------------------------#
# Storage lists for estimated models
mod_list <- list()

# Rename lag to phi so matches in all models:
irq_halfyr <- irq_halfyr %>% rename(phi = p_s1_d_lag)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SILVERMAN REPLICATION                                                   ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# FE Model setup
# ----------------------------------- #
f <- as.formula(paste("p_s1_d ~ -1 + p_spentcptotal_d + p_spentruzicka_d +
                       coalitioncc_d + insurgentcc_d +
                       p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
                       a_of_batt_d + cmoc + dis_usprt +",
                      paste(sprintf("su_vh%s", 2:10), collapse = " + "), " + ",
                      "district_name + time_id"))
# ----------------------------------- #


# ----------------------------------- #
# Silverman - exact
# ----------------------------------- #
# Silverman drops Bagdhad
m1_dat <- irq_halfyr %>% filter(district_name != "Al-Karkh")

mod_list[["m1_silver"]] <- lm(formula = f,
                              data    = m1_dat,
                              weights = m1_dat$pop)
# ----------------------------------- #


# ----------------------------------- #
# Silverman - drop dv lag obs
# ----------------------------------- #
# Including lag drops observations
dat <- m1_dat %>% drop_na(phi)

mod_list[["m2_silver_sub"]] <- lm(formula = f,
                                  data    = dat,
                                  weights = dat$pop)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# STAR                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Create subset w & lw
# ----------------------------------- #
w      <- sp_wts$hy_w
dat$id <- paste(dat$district_name, dat$time_id, sep = ".")
id     <- colnames(w) %in% dat$id
w      <- w[id, id]
w      <- as.matrix(w)
w[w>0] <- 1
lw     <- spdep::mat2listw(w, row.names = row.names(w), style = "W")
w      <- spdep::listw2mat(lw)
# e     <- eigen(w, only.values = T)$values
rm(id)
# ----------------------------------- #


# ----------------------------------- #
# LM Test
# ----------------------------------- #
star_lmtst <- spdep::lm.LMtests(model = mod_list$m2_silver_sub,
                                listw = lw,
                                test  = c("LMerr", "LMlag",
                                          "RLMerr", "RLMlag"))
star_lmtst <- summary(star_lmtst)$results

star_lmtst <- createTexreg(
  coef.names = rownames(star_lmtst),
  coef       = star_lmtst$statistic %>% unlist,
  pvalues    = star_lmtst$p.value
)
# ----------------------------------- #


# ----------------------------------- #
# STAR model
# ----------------------------------- #
f <- p_s1_d ~ phi + p_spentcptotal_d + p_spentruzicka_d +
  coalitioncc_d + insurgentcc_d +
  p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
  a_of_batt_d + cmoc + dis_usprt

mod_list[["m3_star"]] <- lagsarlm(formula = f,
                                  data    = dat,
                                  listw   = lw)

# LRSS
sig <- vcov(mod_list$m3_star)
bs  <- mod_list$m3_star %>% coef
bs_names <- bs[!names(bs) %in% c("rho","(Intercept)","phi")] %>% names
n   <- nrow(dat %>% filter(time_id == "2005h1"))
k   <- length(bs_names)
I   <- diag(n)
wcs <- w[1:nrow(I), 1:ncol(I)]

sims <- 1000   # Set the number of iterations to simulate

# Create array to store n*n matrices over 1000 sims
sim.effs <- matrix(data     = 0,
                   nrow     = sims,
                   ncol     = k,
                   dimnames = list(NULL, bs_names))

# LRSS - Simulate
for(i in 1:sims){
  draws <- MASS::mvrnorm(n     = 1,
                         mu    = bs,
                         Sigma = sig)

  sim.rho   <- draws["rho"]
  sim.betas <- draws[!names(draws) %in% c("rho", "(Intercept)", "phi")]
  sim.phi   <- draws["phi"]

  sim.m  <- solve(I - sim.rho * wcs - sim.phi * I)
  # sim.bm <- sim.betas * sim.m
  sim.bm <- lapply(sim.betas, function(x){
    res <- sum(x * sim.m) / n
    return(res)
  })

  sim.effs[i,] <- unlist(sim.bm)
}

# LRSS - Tidy simulation results
sim_med <- apply(sim.effs, 2, median)
sim_res <- apply(sim.effs, 2, quantile, probs = c(0.025, 0.975)) %>%
  t %>%
  as.data.frame %>%
  rename(Lb = 1, Ub = 2) %>%
  mutate(Med = sim_med) %>%
  dplyr::select(Lb, Med, Ub)

sim_res <- createTexreg(
  coef.names = rownames(sim_res),
  coef       = sim_res$Med,
  ci.low     = sim_res$Lb,
  ci.up      = sim_res$Ub
)

mod_tex <- lapply(mod_list, tidy_mod)
mod_tex[["m3_star_lrss"]] <- sim_res

rm(i, sig, bs, bs_names, n, k, I, wcs, sims, sim_med, sim.effs, sim_res,
   draws, sim.rho, sim.betas, sim.phi, sim.m, sim.bm)
# ----------------------------------- #


# ----------------------------------- #
# TESTING inla star
# ----------------------------------- #
# inla.doc("slm")
# # re.idx = which(abs(Im(e)) < 1e-6)
# # rho.max = 1/max(Re(e[re.idx]))
# # rho.min = 1/min(Re(e[re.idx]))
# # rho = mean(c(rho.min, rho.max))
# f <- p_s1_d ~ p_s1_d_lag + p_spentcptotal_d + p_spentruzicka_d +
#   coalitioncc_d + insurgentcc_d +
#   p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
#   a_of_batt_d + cmoc + dis_usprt
# X <- model.matrix(f, dat)
#
# hyper = list(
#   prec = list(
#     prior = "loggamma",
#     param = c(0.01, 0.01)),
#   rho = list(
#     initial=0,
#     prior = "logitbeta",
#     param = c(1,1)))
#
# Q.beta = Diagonal(n=ncol(X), 0.0001)
#
#
# # Zero-variance for error term
# zero.variance <- list(prec=list(initial = 25, fixed=TRUE))
#
# dat <- dat %>% group_by(time_id) %>% mutate(id2 = 1:n()) %>% ungroup
#
# # Model formula
# fm <- p_s1_d ~ p_s1_d_lag + p_spentcptotal_d + p_spentruzicka_d +
#   coalitioncc_d + insurgentcc_d +
#   p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
#   a_of_batt_d + cmoc + dis_usprt +
#   f(id2,
#     model="slm",
#     args.slm=list(rho.min = -0.95,
#                   rho.max = 0.95,
#                   W=w,
#                   X=X,
#                   Q.beta=Q.beta),
#     hyper=hyper)
#
# slmm1 <- inla(formula = fm,
#               data    = dat,
#               family  = "gaussian",
#               control.fixed = list(correlation.matrix = TRUE))
#               # control.family = list(hyper=zero.variance),
#               # control.compute=list(dic=TRUE, cpo=TRUE))
#
# summary(slmm1)
#
# bs <- bind_rows(slmm1$summary.fixed,
#                 slmm1$summary.hyperpar["Rho for id2",]) %>%
#   tibble::rownames_to_column(var = "coef") %>%
#   mutate(coef = case_when(coef == "Rho for id2" ~ "rho", TRUE ~ coef)) %>%
#   select(coef, mode, contains("quant"))
#
# mod_tex[["m4_star_inla"]] <- createTexreg(
#   coef.names = bs$coef,
#   coef       = bs$mode,
#   ci.low     = bs$`0.025quant`,
#   ci.up      = bs$`0.975quant`
# )
# screenreg(mod_tex, omit.coef = "district_|time_|su_vh")
#
# sims <- 1e3
# library(MASS)
# sig <- slmm1$misc$lincomb.derived.covariance.matrix
#
# # Rho is not in sig - how to get vcv of rho wrt betas?
# for(i in 1:sims){
#   simbs <- mvrnorm(n = 1,
#                    mu = bs$`0.5quant`,
#                    Sigma = sig,
#                    empirical = T)
# };rm(i)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# INLA                                                                    ----
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# MESH                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Boundary Matrix and district coordinates
# ----------------------------------- #
irq1 <- irq %>% st_transform(., crs = "+proj=longlat")

irq0 <- irq1 %>%
  st_union() %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  st_set_geometry(.,"geometry")

irq_bound     <- as_Spatial(irq0)
irq_bound     <- inla.sp2segment(irq_bound)
irq_bound$loc <- inla.mesh.map(irq_bound$loc)
cords         <- cbind(dat$lat, dat$long)

rm(irq)
# ----------------------------------- #


# ----------------------------------- #
# Construct mesh
# ----------------------------------- #
# NOTES:
# max.edge ~ This sets maximum allowed triangle length in decimal degrees.
# Smaller numbers lead to more precision with a higher computational cost.
# More precision does not mean more accuracy!
# This is a tuple attribute taking 2 options:
# c(max inside boundary, max outside boundary)

# offset ~ This defines the distance (decimal degrees) beyond the inner
# domain (i.e., secondary boundary box).

# cutoff ~ to minimize too many triangles around tightly clustered (small)
# districts

# Nb. a mesh can be constructed with only the coordinates of spatial locations
# and no boundary resulting in a mesh based on a convex hull.

# At Baghdad Latitude (33.3152), 1 degree of longitude is approximately:
# (111.320*cos(33.3152 * pi / 180)) # km

# 20 km is equal to the following decimal degrees: [0.215] (inside edge)
# 20 / (111.320*cos(33.3152 * pi / 180))

# 200 km is equal to the following decimal degrees: [2.150] (outer edge)
# 200 / (111.320*cos(33.3152 * pi / 180))

# 50 km is equal to the following decimal degrees: [0.540] (edge cutoff)
# 25 / (111.320*cos(33.3152 * pi / 180))

# Create mesh
mesh <- inla.mesh.2d(loc      = cords,
                     max.edge = c(0.25, 2.150),
                     offset   = c(0.25, 2.150),
                     cutoff   = .75,
                     boundary = irq_bound)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SPDE                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# SPDE prior justification
# ----------------------------------- #
# Using the average (decimal degree) distance from Baghdad (location with
# most violence in conflict) to international border used to constrain mesh

# # --- Prior justification
# # Nb on priors - 1 degree at equator = 111 km
# # Baghdad saw most violence. Average distance from Baghdad to border:
# baghdad <- st_sfc(st_point(c(44.3661, 33.3152)),
#                   crs = st_crs("+proj=longlat")) %>%
#   st_transform(., st_crs(irq1))
#
# irq_border <- irq0 %>%
#   st_union() %>%
#   st_cast(., "POINT")
#
# res <- st_distance(irq_border, baghdad)
# summary(res / 1e3)
# # Mean = 331.7 km from border
# mean(res) %>% as.numeric
# # Mean = 331663.9 meters from border
# For unprojected data this is:
# # At Baghdad Latitude (33.3152), 1 degree of longitude is approximately:
# 331.7 / (111.320*cos(33.3152 * pi / 180))
# # ~ 3.56 decimal degrees # Unprojected
# ----------------------------------- #


# ----------------------------------- #
# Create spde model object with penalized complexity priors
# ----------------------------------- #
# [see: Fuglstad et al. 2018]
#  ~ prior.range = c(0.5, 0.01) ~ # p(range > 0.05) = 0.01
#  ~ prior.sigma = c(1, 0.01)   ~ # p(sigma > 1)    = 0.01
spde <- inla.spde2.pcmatern(mesh        = mesh,
                            alpha       = 2,  # Fractional operator, 2=default
                            prior.range = c(3.6, 0.1),
                            prior.sigma = c(001, 0.1))
# ----------------------------------- #


# ----------------------------------- #
# SPDE index object
# ----------------------------------- #
# Additional data preparation is required to build the space-time model. The
# index set is made taking into account the number of mesh points in the
# SPDE model and the number of groups, as:
dat$time_id  <- fct_drop(dat$time_id)
time_periods <- length(unique(dat$time_id))
iset <- inla.spde.make.index(name    = "i",
                             n.spde  = spde$n.spde,
                             n.group = time_periods)
# ----------------------------------- #


# ----------------------------------- #
# Projection matrix
# ----------------------------------- #
# The "A" matrix maps the GMRF from the n mesh nodes back to the observed
# data locations (A is a matrix with columns referring to observations in the
# data and rows corresponding to nodes in the mesh)

# The index set for the latent field does not depend on the data set
# locations. It only depends on the SPDE model size and on the time dimension.
# The projection matrix is defined using the coordinates of the observed data.
# We need to pass the time index to the group argument to build the projector
# matrix and the inla.spde.make.A() function:
A <- inla.spde.make.A(mesh  = mesh,
                      loc   = as.matrix(cords),
                      group = as.numeric(dat$time_id))
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# STACKS                                                                  ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Stack note
# ----------------------------------- #
# The A matrix should only be applied to spatial random effects, not fixed
# (i.e., observed covariates) which were observed at coordinate locations.
# The inla.stack() facilitates these restrictions.
#
# Below this is accomplished with the argument: A = list(A, 1) which provides
# "A" as the projector matrix for the latent field and 1 as a simple 1-to-1
# map for the covariates in a model
# ----------------------------------- #


# ----------------------------------- #
# Create Stacks
# ----------------------------------- #
data_stack <- inla.stack(
  data = list(sigacts = dat$p_s1_d
              # icews   = dat$p_icews_d,
              # ged     = dat$p_ged_d),
  ),
  A    = list(A, 1),
  effects = list(
    c(list(Intercept = 1), iset),
    list(phi        = dat$phi,
         # icews_l          = dat$p_icews_d_lag,
         # ged_l            = dat$p_ged_d_lag,
         p_spentcptotal_d = dat$p_spentcptotal_d,
         p_spentruzicka_d = dat$p_spentruzicka_d,
         coalitioncc_d    = dat$coalitioncc_d,
         insurgentcc_d    = dat$insurgentcc_d,
         p_spentcerpsmall_noncp_d  = dat$p_spentcerpsmall_noncp_d,
         p_spentusaid_nonruzicka_d = dat$p_spentusaid_nonruzicka_d,
         a_of_batt_d = dat$a_of_batt_d,
         cmoc        = dat$cmoc,
         dis_usprt   = dat$dis_usprt)),
  tag = "data_stack"
)
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# PRIORS                                                                  ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Temporal prior
# ----------------------------------- #
# PC temporal autoregressive prior [considers P(cor>0) = 0.9]
# h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))

# PC temporal prior on random walk:
h.spec <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))
# ----------------------------------- #

# ----------------------------------- #
# Gaussian precision prior
# ----------------------------------- #
# PC prior on gaussian precision
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# FORMULA                                                                 ----
#-----------------------------------------------------------------------------#
f <- . ~ - 1 + Intercept + p_spentcptotal_d + p_spentruzicka_d +
  coalitioncc_d + insurgentcc_d +
  p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
  a_of_batt_d + cmoc + dis_usprt +
  f(i,
    model         = spde,
    group         = i.group,
    control.group = list(model = "rw1", hyper = h.spec))

# change model to "ar1" for autoregressive
# Nb - for rw1 there will be no gmrf temporal parameter akin to Python et al
# 2018 Table 1, p10
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# INLA MODELS                                                             ----
#-----------------------------------------------------------------------------#

# With temporal lag
mod_list[["m4_1_spde"]] <- inla(formula = update(f, sigacts ~ phi + .),
                 family  = "gaussian",
                 data    = inla.stack.data(data_stack, spde = spde),
                 control.predictor = list(A       = inla.stack.A(data_stack),
                                          compute = TRUE),
                 control.family    = list(hyper   = list(theta = prec.prior),
                                          link    = "identity"),
                 control.compute   = list(config  = TRUE))

inla_prms <- inla_params(model_list = list("sigacts" = mod_list$m4_1_spde),
                         spde       = spde)

llk <- mod_list$m4_1_spde$mlik[2]
aic <- l_aic(llik = llk, k = mod_list$m4_1_spde$neffp[1])
bic <- l_bic(llik = llk, k = mod_list$m4_1_spde$neffp[1], n = 824)

mod_tex[["m4_1_spde"]] <- createTexreg(
  coef.names = inla_prms$sigacts$variable[1:14],
  coef       = inla_prms$sigacts$median[1:14],
  ci.low     = inla_prms$sigacts$lb[1:14],
  ci.up      = inla_prms$sigacts$ub[1:14],
  gof.names  = c("n", "LogLik", "AIC", "BIC"),
  gof        = c(824, llk, aic, bic),
  gof.decimal= c(F, T, T, T)
)

# Without temporal lag
mod_list[["m4_2_spde"]] <- inla(formula = update(f, sigacts ~ + .),
                family  = "gaussian",
                data    = inla.stack.data(data_stack, spde = spde),
                control.predictor = list(A       = inla.stack.A(data_stack),
                                         compute = TRUE),
                control.family    = list(hyper   = list(theta = prec.prior),
                                         link    = "identity"),
                control.compute   = list(config  = TRUE))

inla_prms <- inla_params(model_list = list("sigacts" = mod_list$m4_2_spde),
                         spde       = spde)

llk <- mod_list$m4_2_spde$mlik[2]
aic <- l_aic(llik = llk, k = mod_list$m4_2_spde$neffp[1])
bic <- l_bic(llik = llk, k = mod_list$m4_2_spde$neffp[1], n = 824)

mod_tex[["m4_2_spde"]] <- createTexreg(
  coef.names = inla_prms$sigacts$variable[1:13],
  coef       = inla_prms$sigacts$median[1:13],
  ci.low     = inla_prms$sigacts$lb[1:13],
  ci.up      = inla_prms$sigacts$ub[1:13],
  gof.names  = c("n", "LogLik", "AIC", "BIC"),
  gof        = c(824, llk, aic, bic),
  gof.decimal= c(F, T, T, T)
)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# INLA Model Fields
#-----------------------------------------------------------------------------#

# The inclusion or exclusion of the dynamic lag has no meaningful difference
# in regard to GRMF values.
m4_1_field <- inla_fields(model = mod_list$m4_1_spde,
                          mesh  = mesh,
                          index = iset,
                          time_periods = time_periods,
                          boundary = irq0)

# m4_2_field <- inla_fields(model = mod_list$m4_2_spde,
#                           mesh  = mesh,
#                           index = iset,
#                           time_periods = time_periods,
#                           boundary = irq0)

rw_1_mean <- ggplot(data = m4_1_field$mean) +
  geom_raster(aes(x = x, y = y, fill = val)) +
  geom_sf(data = irq1, fill = "transparent", color = "gray40", size = 0.1) +
  geom_sf(data = irq0, fill = "transparent", color = "black",  size = 0.4) +
  scale_fill_gradient2(low      = muted("blue"),
                       high     = muted("red"),
                       mid      = "white",
                       midpoint = 0,
                       limits   = c(-0.2, 0.85)) +
  # labs(title = expression(paste("GMRF, ", zeta))) +
  map_theme +
  theme(legend.text     = element_text(size = 16/.pt),
        legend.key.size = unit(5, "mm"))

windows(height = 6.5, width = 6.5);rw_1_mean
dev.off()

ggsave(plot     = rw_1_mean,
       filename = "Results/Figures/figure-gmrf_mean.png",
       width    = 5.5,
       height   = 5.0,
       units    = "in")

file.copy(
  from = "Results/Figures/figure-gmrf_mean.png",
  to   = "Drafts/Drafts/GAST20210906/figure-gmrf_mean.png",
  overwrite = TRUE
)

# m4_1_field$sd$val_ln <- log(m4_1_field$sd$val)
# rw_1_sigma <- ggplot(data = m4_1_field$sd) +
#   geom_raster(aes(x = x, y = y, fill = val_ln)) +
#   geom_sf(data = irq0, fill = "transparent", color = "black", size = 0.1) +
#   scale_fill_gradient2(low      = muted("blue"),
#                        high     = muted("red"),
#                        mid      = "white",
#                        midpoint = median(m4_1_field$sd$val_ln),
#                        limits   = NULL) +
#   labs(subtitle = "Sigma (logged)") +
#   map_theme

# rw_2_mean <- ggplot(data = m4_2_field$mean) +
#   geom_raster(aes(x = x, y = y, fill = val)) +
#   geom_sf(data = irq0, fill = "transparent", color = "black", size = 0.1) +
#   scale_fill_gradient2(low      = muted("blue"),
#                        high     = muted("red"),
#                        mid      = "white",
#                        midpoint = 0,
#                        limits   = NULL) +
#   labs(subtitle = "Mean") +
#   map_theme
#
# m4_2_field$sd$val_ln <- log(m4_2_field$sd$val)
#
# rw_2_sigma <- ggplot(data = m4_2_field$sd) +
#   geom_raster(aes(x = x, y = y, fill = val_ln)) +
#   geom_sf(data = irq0, fill = "transparent", color = "black", size = 0.1) +
#   scale_fill_gradient2(low      = muted("green"),
#                        high     = muted("red"),
#                        mid      = "white",
#                        midpoint = median(m4_2_field$sd$val_ln),
#                        limits   = NULL) +
#   labs(subtitle = "Sigma, logged") +
#   map_theme

# rw_field <- cowplot::plot_grid(rw_mean, rw_sigma, align = "hv", axis = "l",
#                                ncol = 2)
# ggsave(plot     = rw_field,
#        filename = "Results/Figures/figure-field-sigacts-rw.png",
#        width    = 6.5,
#        height   = 4.0,
#        units    = "in")
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# GMRF RANGE ESTIMATES                                                    ----
#-----------------------------------------------------------------------------#
matern <- function(lambda, kappa, dist){
  # Matern covariance function, used with matern_data()
  2^(1-lambda)/gamma(lambda) * (kappa*dist)^lambda* besselK(x=dist*kappa, nu=lambda)
}


matern_data <- function(model,
                        max_distance = NULL){
  # ----------------------------------- #
  # function description
  # ----------------------------------- #
  # matern_data takes an inla model where random effects estimated in a gaussian field are assumed
  # to follow a continuous decay process fit with a stochastic partial differential equation with a
  # matern covariance solution. This function extract model parameters to construct decay estimates
  # up to a max_distance argument (km) set by the user.
  # ----------------------------------- #


  # ----------------------------------- #
  # Extract spatial field and parameters from model input
  # ----------------------------------- #
  # Spatial parameters with nominal scale (time is aggregated if present)
  spde.resfinal <- inla.spde2.result(inla = model,
                                     name = "i",
                                     spde,do.transform = TRUE) # do.transform put in correct scale

  # Kappa / computed using Blangiardo & Cameletti (2015)
  Kappa    <-inla.emarginal(function(x) x, spde.resfinal$marginals.kappa[[1]]) # kappa (mean)
  Kappahpd <-inla.hpdmarginal(0.95, spde.resfinal$marginals.kappa[[1]])        # kappa (hpd 95%)

  # Variance of the random field
  variance    <- inla.emarginal(function(x) x, spde.resfinal$marginals.variance.nominal[[1]]) # variance (mean)
  variancehpd <- inla.hpdmarginal(0.95, spde.resfinal$marginals.variance.nominal[[1]])        # variance (hpd 95%)

  # Range + degree conversion
  range    <- inla.emarginal(function(x) x, spde.resfinal$marginals.range.nominal[[1]]) # range (mean)
  rangehpd <- inla.hpdmarginal(0.95, spde.resfinal$marginals.range.nominal[[1]])        # range (hpd 95%)
  deg      <- 2*pi*6371/360
  # ----------------------------------- #


  # ----------------------------------- #
  # Build plot data frame
  # ----------------------------------- #
  # Additional parameters:
  lambda <- 1

  if(is.null(max_distance)){
    dist.x <- seq(from = 0.01, to = rangehpd[,"high"], length.out = 100)
  } else{
    # Convert user-input kilometers "max_distance" to decimal degrees for
    # appropriate mapping to estimation space.
    max_distance <- max_distance / deg

    dist.x <- seq(from = 0.01, to = max_distance, length.out = 100)
  }

  plt_dat <- data.frame(
    "x"  = dist.x,
    "lb" = matern(lambda = lambda,
                  kappa  = Kappahpd[,"low"],
                  dist   = dist.x),
    "y"  = matern(lambda = lambda,
                  kappa  = Kappa,
                  dist   = dist.x),
    "ub" = matern(lambda = lambda,
                  kappa  = Kappahpd[,"high"],
                  dist   = dist.x))
  # ----------------------------------- #


  # ----------------------------------- #
  # Return statement
  # ----------------------------------- #
  res <- list("df"    = plt_dat,
              "range" = range)

  return(res)
  # ----------------------------------- #
}

range_est <- matern_data(model = mod_list$m4_1_spde, max_distance = 300)

# Breaks: 0 - 500km (converted to decimal degrees)
breaks <- c(0, 50, 100, 200, 300) / (2*pi*6371/360)
breaks <- sort(c(breaks, range_est$range))
labs   <- (breaks * (2*pi*6371/360)) %>% round(., 0)

rw_1_range <- ggplot(data = range_est$df, aes(x=x)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "gray50", alpha = 0.2) +
  geom_line(aes(y = y),  color = "black", size = 0.25) +

  geom_vline(aes(xintercept = range_est$range), linetype = "dashed", size = 0.2, color = "gray50") +
  geom_hline(aes(yintercept = 0.05),   linetype = "dashed", size = 0.2, color = "gray50") +

  scale_x_continuous(name   = "Distance [km]",
                     limits = c(0,max(breaks)),
                     breaks = breaks,
                     labels = labs,
                     expand = c(0,0)) +

  scale_y_continuous(name   = "Matern covariance function",
                     limits = c(0,1),
                     breaks = c(seq(0, 1, 0.25), 0.1)) +

  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title     = element_blank(),
        panel.grid       = element_blank(),
        axis.text        = element_text(size = 20/.pt),
        panel.background = element_rect(fill = NA, color = "black", size = 0.1),
        strip.background = element_rect(fill = "gray95", color = "black", size = 0.1),
        plot.margin      = unit(c(1,3,1,1), "mm"))
  # labs(title    = "Spatial correlation decay",
  #      subtitle = "SPDE - Insurgent Events")


windows(height = 3.5, width = 6.5);rw_1_range
dev.off(2)

ggsave(plot     = rw_1_range,
       filename = "Results/Figures/figure-spde_range.png",
       width    = 6.5,
       height   = 3.5,
       units    = "in")

file.copy(
  from = "Results/Figures/figure-spde_range.png",
  to   = "Drafts/Drafts/GAST20210906/figure-spde_range.png",
  overwrite = TRUE
)

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# PRINT MODELS                                                            ----
#-----------------------------------------------------------------------------#
screenreg(mod_tex,
          omit.coef = "district_|time_|su_vh|Intercept",
          custom.gof.rows = list("FE: District" = c("Yes", "Yes",
                                                    "No", "", "No", "No"),
                                 "FE: Time"     = c("Yes", "Yes",
                                                    "No", "", "No", "No"),
                                 "FE: Sunni VS" = c("Yes", "Yes",
                                                    "No", "", "No", "No")))
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
save(dat,m1_dat, m4_2_field, mod_list,mod_tex, spde, mesh,
     file = "Results/Estimates/results-apsa.Rdata")
#rm(list = ls())
#-----------------------------------------------------------------------------#
