

library(SpatialNLS)
library(spatialreg)
library(texreg)

load("../Data/Data-All-Tidy-20210131.Rdata")
load("../Data/Data-Spatial-Weights-20210131.Rdata")

# ----------------------------------- #
# Spatial weights
# ----------------------------------- #
tmp <- as_Spatial(d_cs)
tmp <- spdep::poly2nb(tmp)
W_cs <- spdep::nb2mat(tmp, style = "W")
W_lw <- spdep::mat2listw(W_cs)
rm(tmp)
# ----------------------------------- #

# ----------------------------------- #
# DATA CS
# ----------------------------------- #
covars <- d_pn_mo %>% st_drop_geometry() %>%
  select(Year, Month, District, Pop_den, Unemp_rate, starts_with("CD")) %>%
  group_by(District) %>%
  summarize(Pop_den = mean(Pop_den, na.rm = T),
            Unemp_rate = mean(Unemp_rate, na.rm = T),
            across(starts_with("CD"), max))
dat <- left_join(d_cs, covars, by = "District")
# ----------------------------------- #


#-----------------------------------------------------------------------------#
# MODELS                                                                  ----
#-----------------------------------------------------------------------------#
form <- SIGACT_pc ~ Pop_den + Unemp_rate + CD_Coalition + CD_Insurgent + CD_Sectarian

# ----------------------------------- #
# OLS
# ----------------------------------- #
mod.ols <- lm(formula = form, data = dat)
# ----------------------------------- #

# ----------------------------------- #
# SPATIAL NLS
# ----------------------------------- #
mod.nls <- lagsarnls(formula = form,
                     data    = dat,
                     W       = W_cs)
# ----------------------------------- #

# ----------------------------------- #
# SPATIAL REG
# ----------------------------------- #
mod.spr <- lagsarlm(formula = form, data = dat, listw = W_lw)
# ----------------------------------- #

#-----------------------------------------------------------------------------#
# MODEL OUTPUT                                                            ----
#-----------------------------------------------------------------------------#
names(coefficients(mod.ols))
names(coefficients(mod.nls))
names(coefficients(mod.spr))

m_vnames = list("Pop_den"           = "Pop. Density",
                "b_Pop_den"         = "Pop. Density",

                "Pop_urban_pct"     = "Pop. Urban",
                "b_Pop_urban_pct"   = "Pop. Urban",

                "Unemp_rate"        = "Unemployment",
                "b_Unemp_rate"      = "Unemployment",

                "CD_Coalition"      = "CD Coalition",
                "CD_Insurgent"      = "CD Insurgent",
                "CD_Sectarian"      = "CD Sectarian",

                "b_CD_Coalition"      = "CD Coalition",
                "b_CD_Insurgent"      = "CD Insurgent",
                "b_CD_Sectarian"      = "CD Sectarian",

                "rho"               = "rho",

                "(Intercept)"       = "Intercept",
                "b_Intercept"       = "Intercept")


screenreg(l = list(mod.ols, mod.spr, mod.nls),
          custom.model.names = rep(c("OLS","SPR","NLS"),1),
          # custom.header      = list("No FEs" = 1:3, "FEs" = 4:6),
          custom.coef.map    = m_vnames)

