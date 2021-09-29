#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          July 03, 2021
# Purpose:       Space-time dynamic modeling with R-INLA and SPDE
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# Notes:
#
# Code to simulate and estimate a "space-time coregionalizaton model"
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html
#
# Check - https://github.com/inlabru-org/inlabru [already have package] for
#         simpler code setup
#
# Still struggling - these all specify spatial random effects. The temporal
# dependence is structured on the fields - that is, the fields vary from
# time_i to time_j by dependece rho.
#
# See : https://ourcodingclub.github.io/tutorials/inla/
# for more on this spatial field temporal dependence.
#
# However - spatial fields are NOT the same as spatial dynamics.
#
# A tutorial - nb A. Python
#      https://punama.github.io/BDI_INLA/
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
library(forcats)
library(sf)
library(INLA)
#---------------------------#

#---------------------------#
# Load data
#---------------------------#
load("data/Oxford_HB_2021_APSA-Data.Rdata")
#---------------------------#

#---------------------------#
# Functions
#---------------------------#
inla_fit <- function(dv, bias = FALSE){
  cat(sprintf("\14Working on: %s", dv))
  dvl   <- paste(dv,"l", sep = "_")

  if(bias){
    tmp_f <- update(f, paste(dv, " ~ + sigacts2 + ."))
  } else{
    tmp_f <- update(f, paste(dv, " ~ + ", dvl, " + ."))
  }

  res   <- inla(formula = tmp_f,
                family  = "gaussian",
                data    = inla.stack.data(data_stack, spde = spde),
                control.predictor = list(A       = inla.stack.A(data_stack),
                                         compute = TRUE),
                control.family    = list(hyper   = list(theta = prec.prior),
                                         link    = "identity"),
                control.compute   = list(dic     = TRUE,
                                         waic    = TRUE,
                                         config  = TRUE))
  return(res)
}
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# DATA                                                                    ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# To facilitate easily switching between half-year or monthly, set here.
# ----------------------------------- #
dat <- irq_halfyr %>%
  rename(
    sigacts   = p_s1_d,
    sigacts_l = p_s1_d_lag,
    icews     = p_icews_d,
    icews_l   = p_icews_d_lag,
    ged       = p_ged_d,
    ged_l     = p_ged_d_lag
  )
rm(irq_halfyr, irq_month, sp_wts)
# ----------------------------------- #


# ----------------------------------- #
# Set time length as a variable to reuse
# ----------------------------------- #
time_periods <- length(unique(dat$time_id))
# ----------------------------------- #
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
                     max.edge = c(0.215, 2.150),
                     offset   = c(0.215, 2.150),
                     cutoff   = 0.27,
                     boundary = irq_bound)

png(filename = "results/figures/figure-mesh.png",
    width    = 6.5,
    height   = 6.5,
    res      = 350,
    units    = "in")
plot(mesh, main = "")
dev.off()
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

st_centroid(irq0)
prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +datum=NAD83 +units=m +no_defs"

x   <- irq1$district_name
res <- lapply(x, function(id){
  tmp <- irq1 %>% filter(district_name == id) %>% select(district_name)
  cnt <- st_centroid(tmp)

})


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
# Create DV vector
# ----------------------------------- #
dvs <- c("sigacts", "icews", "ged")
# ----------------------------------- #


# ----------------------------------- #
# Create Stacks
# ----------------------------------- #
data_stack <- inla.stack(
  data = list(sigacts = dat$sigacts,
              icews   = dat$icews,
              ged     = dat$ged),
  A    = list(A, 1),
  effects = list(
    c(list(Intercept = 1), iset),
    list(sigacts2         = dat$sigacts,
         sigacts_l        = dat$sigacts_l,
         icews_l          = dat$icews_l,
         ged_l            = dat$ged_l,
         p_spentcptotal_d = dat$p_spentcptotal_d,
         p_spentruzicka_d = dat$p_spentruzicka_d,
         coalitioncc_d    = dat$coalitioncc_d,
         insurgentcc_d    = dat$insurgentcc_d,
         p_spentcerpsmall_noncp_d  = dat$p_spentcerpsmall_noncp_d,
         p_spentusaid_nonruzicka_d = dat$p_spentusaid_nonruzicka_d,
         a_of_batt_d = dat$a_of_batt_d,
         cmoc        = dat$cmoc,
         dis_usprt   = dat$dis_usprt,
         baghdad     = dat$baghdad,
         time        = dat$time_id)),
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
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# MODELS                                                                  ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# Theoretical (conflict) models
# ----------------------------------- #
{conflict        <- lapply(dvs, inla_fit)
names(conflict) <- dvs
# ----------------------------------- #


# ----------------------------------- #
# Bias models
# ----------------------------------- #
bias        <- lapply(dvs[2:3], inla_fit, bias = TRUE)
names(bias) <- paste0(dvs[2:3],"_bias")
}# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
# Clean-up
rm(data_stack, h.spec, prec.prior, f, inla_fit)

# Save
save.image(file = "Results/Estimates/results-spde.Rdata")
#-----------------------------------------------------------------------------#

rm(list = ls())
