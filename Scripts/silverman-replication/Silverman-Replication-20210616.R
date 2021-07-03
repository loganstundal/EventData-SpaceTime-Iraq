#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          June 16, 2021
# Purpose:       Revised silverman analysis
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
library(texreg)
library(sandwich)
library(lmtest)

#---------------------------#
# Set working directory
#---------------------------#
# setwd()

#---------------------------#
# Load data
#---------------------------#
silver <- readstata13::read.dta13(
  paste0("Scripts/silverman-replication/",
         "silverman-replication-data/CP_Data_Final_Half-Year -- 3-29-20.dta"))

#---------------------------#
# Local Functions
#---------------------------#
clse_pvs <- function(model, data, cluster){
  # clse_pvs computes clustered standard errors and associated p-valuse
  # for a model object and returns these as lists to use with texreg
  clse <- sandwich::vcovCL(model, cluster = data[[cluster]])
  clse <- sqrt(diag(clse))

  bs   <- coef(model)[!is.na(coef(model))]

  pvls <- 2 * (1 - pnorm(abs(bs / clse)))

  res <- list("clse" = clse, "pvals" = pvls)
  return(res)
}

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# TIDY DATA                                                               ----
#-----------------------------------------------------------------------------#
"Task 1 - need to construct first differences and preserve original FEs in
          Silverman's data. His models include the following variables:

      - DV    ~ p_S1                            | differenced,
      - IVs   ~ p_spentcptotal & p_spentruzicka | differenced,
      - Ctrls ~ coalitioncc, insurgentcc,
                p_spentcerpsmall_noncp,
                p_spentusaid_nonruzicka,
                a_of_batt, cmoc, dis_usprt      | differenced,
      - FEs   ~ halfyr2-10, su_vh2-10, district"

silver <- silver %>%
  # Select relevant variables
  select(district, districtid, halfyr,
         p_S1, p_spentcptotal, p_spentruzicka,
         coalitioncc, insurgentcc,
         p_spentcerpsmall_noncp, p_spentusaid_nonruzicka,
         a_of_batt, cmoc, dis_usprt, POP,
         starts_with(c("halfyr", "su_vh"))) %>%

  # Generate appropriate first differences and drop NAs on DV
  group_by(district) %>%
  mutate(across(!starts_with(c("halfyr", "su_vh", "POP","cmoc","dis_usprt","districtid")),
                .fns = list(d = ~ . - dplyr::lag(., 1)))) %>%

  tidyr::drop_na(p_S1_d) %>%

  # Create a temporal lag of SIGACTS levels for future use
  mutate(p_S1_lag = dplyr::lag(p_S1, 1),
         p_S1_d_lag = dplyr::lag(p_S1_d)) %>%

  # Final variable select
  select(district, districtid, halfyr,
         p_S1, p_S1_lag,
         p_S1_d, p_S1_d_lag,
         p_spentcptotal, p_spentcptotal_d,
         p_spentruzicka, p_spentruzicka_d,
         coalitioncc, coalitioncc_d,
         insurgentcc, insurgentcc_d,
         p_spentcerpsmall_noncp, p_spentcerpsmall_noncp_d,
         p_spentusaid_nonruzicka, p_spentusaid_nonruzicka_d,
         a_of_batt, a_of_batt_d,
         cmoc, dis_usprt, POP,
         starts_with(c("halfyr", "su_vh")))
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# ANALYSIS - Table 2 replication                                          ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Base formula and data subset
# ----------------------------------- #
f <- as.formula(paste("p_S1_d ~ -1 + p_spentcptotal_d + p_spentruzicka_d + ",
                      paste(sprintf("halfyr%s",2:10), collapse = " + "), " + ",
                      paste(sprintf("su_vh%s", 2:10), collapse = " + "), " + ",
                      "district"))
f <- list("m1" = f,
          "m2" = update(f, ~ . + coalitioncc_d + insurgentcc_d),
          "m3" = update(f, ~ . + coalitioncc_d + insurgentcc_d +
                          p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d),
          "m4" = update(f, ~ . + coalitioncc_d + insurgentcc_d +
                          p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
                          a_of_batt_d + cmoc + dis_usprt))
silver2 <- silver %>% filter(district != "Karkh")
# ----------------------------------- #


# ----------------------------------- #
# Models
# ----------------------------------- #
m1 <- lm(formula = f$m1,
         data    = silver2,
         weights = silver2$POP)

m2 <- lm(formula = f$m2,
         data    = silver2,
         weights = silver2$POP)

m3 <- lm(formula = f$m3,
         data    = silver2,
         weights = silver2$POP)

m4 <- lm(formula = f$m4,
         data    = silver2,
         weights = silver2$POP)

# Confidence intervals for M4 primary variables
round(confint(lmtest::coeftest(m4, vcov. = vcovCL(m4, cluster = silver2$district))),3)[1:2,]
# ----------------------------------- #

# Testing sensitivity of m.x to other forms of population weights:
mx <- lm(formula = f$m1, data = silver2, weights = log(silver2$POP))
lmtest::coeftest(mx, vcov. = vcovCL(mx, cluster = silver2$district))[1:2,]
rm(mx)
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# ROBUSTNESS - ALTERNATIVE SPECIFICATIONS                                 ----
#-----------------------------------------------------------------------------#
# ----------------------------------- #
# No weights
# ----------------------------------- #
m1.nw <- lm(formula = f$m1,
            data    = silver2)

m2.nw <- lm(formula = f$m2,
            data    = silver2)

m3.nw <- lm(formula = f$m3,
            data    = silver2)

m4.nw <- lm(formula = f$m4,
            data    = silver2)
# ----------------------------------- #


# ----------------------------------- #
# Include Baghdad
# ----------------------------------- #
m1.bg <- lm(formula = f$m1,
            weights = silver$POP,
            data    = silver)

m2.bg <- lm(formula = f$m2,
            weights = silver$POP,
            data    = silver)

m3.bg <- lm(formula = f$m3,
            weights = silver$POP,
            data    = silver)

m4.bg <- lm(formula = f$m4,
            weights = silver$POP,
            data    = silver)
# ----------------------------------- #


# ----------------------------------- #
# No weights + Baghdad
# ----------------------------------- #
m1.00 <- lm(formula = f$m1,
            data    = silver)

m2.00 <- lm(formula = f$m2,
            data    = silver)

m3.00 <- lm(formula = f$m3,
            data    = silver)

m4.00 <- lm(formula = f$m4,
            data    = silver)

# Confidence intervals for M4 primary variables
round(confint(lmtest::coeftest(m4.00, vcov. = vcovCL(m4.00, cluster = silver$district))),3)[1:2,]
round(confint(lmtest::coeftest(m4.00)),3)[1:2,]
round(coef(m4.00)[1],3)
# ----------------------------------- #

# ----------------------------------- #
# Model 4 with no fixed effects
# ----------------------------------- #
f.nofe <- formula(p_S1_d ~ -1 + p_spentcptotal_d + p_spentruzicka_d + coalitioncc_d + insurgentcc_d +
                    p_spentcerpsmall_noncp_d + p_spentusaid_nonruzicka_d +
                    a_of_batt_d + cmoc + dis_usprt)

m4.no.unit <- lm(formula = update(f.nofe, ~ . +
                                    halfyr2 + halfyr3 + halfyr4 + halfyr5 + halfyr6 + halfyr7 + halfyr8 + halfyr9 + halfyr10 +
                                    su_vh2 + su_vh3 + su_vh4 + su_vh5 + su_vh6 + su_vh7 + su_vh8 + su_vh9 + su_vh10),
                 weights = silver$POP,
                 data = silver)
m4.no.unit <- coeftest(m4.no.unit, vcov = vcovCL(m4.no.unit, cluster = silver$district))[1:9,]
m4.no.unit <- as.data.frame(m4.no.unit)
m4.no.unit <- createTexreg(coef.names = rownames(m4.no.unit),
                           coef       = m4.no.unit$Estimate,
                           se         = m4.no.unit$`Std. Error`,
                           pvalues    = m4.no.unit$`Pr(>|t|)`,
                           gof.names  = c("Num. obs."),
                           gof        = 936)

m4.no.suvs <- lm(formula = update(f.nofe, ~ . + as.factor(district) +
                                    halfyr2 + halfyr3 + halfyr4 + halfyr5 + halfyr6 + halfyr7 + halfyr8 + halfyr9 + halfyr10),
                 weights = silver$POP,
                 data = silver)

m4.no.suvs <- coeftest(m4.no.suvs, vcov = vcovCL(m4.no.suvs, cluster = silver$district))[1:9,]
m4.no.suvs <- as.data.frame(m4.no.suvs)
m4.no.suvs <- createTexreg(coef.names = rownames(m4.no.suvs),
                           coef       = m4.no.suvs$Estimate,
                           se         = m4.no.suvs$`Std. Error`,
                           pvalues    = m4.no.suvs$`Pr(>|t|)`,
                           gof.names  = c("Num. obs."),
                           gof        = 936)

m4.no.time <- lm(formula = update(f.nofe, ~ . + as.factor(district) +
                                    su_vh2 + su_vh3 + su_vh4 + su_vh5 + su_vh6 + su_vh7 + su_vh8 + su_vh9 + su_vh10),
                 weights = silver$POP,
                 data = silver)
m4.no.time <- coeftest(m4.no.time, vcov = vcovCL(m4.no.time, cluster = silver$district))[1:9,]
m4.no.time <- as.data.frame(m4.no.time)
m4.no.time <- createTexreg(coef.names = rownames(m4.no.time),
                           coef       = m4.no.time$Estimate,
                           se         = m4.no.time$`Std. Error`,
                           pvalues    = m4.no.time$`Pr(>|t|)`,
                           gof.names  = c("Num. obs."),
                           gof        = 936)

m4.no.fes <-  lm(formula = f.nofe,
                 weights = silver$POP,
                 data = silver)
m4.no.fes <- coeftest(m4.no.fes, vcov = vcovCL(m4.no.fes, cluster = silver$district))[1:9,]
m4.no.fes <- as.data.frame(m4.no.fes)
m4.no.fes <- createTexreg(coef.names = rownames(m4.no.fes),
                          coef       = m4.no.fes$Estimate,
                          se         = m4.no.fes$`Std. Error`,
                          pvalues    = m4.no.fes$`Pr(>|t|)`,
                          gof.names  = c("Num. obs."),
                          gof        = 936)
# ----------------------------------- #
screenreg(l = list(m4.no.unit, m4.no.time, m4.no.suvs, m4.no.fes),
          custom.coef.map = vars,
          custom.model.names = c("No unit FEs", "No time FEs", "No Sunni VS FEs", "No FEs"))
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# TABLE MODELS                                                            ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Variable Names
# ----------------------------------- #
vars <- list(p_spentcptotal_d          = 'Condolence spending per cap',
             p_spentruzicka_d          = 'Ruzicka spending per cap',
             coalitioncc_d             = 'Coalition collateral damage',
             insurgentcc_d             = 'Insurgent collateral damage',
             p_spentcerpsmall_noncp_d  = 'Other small cerp spending',
             p_spentusaid_nonruzicka_d = 'Other USAID spending',
             a_of_batt_d               = 'Coalition troop strength',
             cmoc                      = 'CMOC presence',
             dis_usprt                 = 'PRT presence')

# mods <- list(m1, m2, m3, m4)
mods <- list(m1, m1.nw, m1.bg, m1.00,
             m2, m2.nw, m2.bg, m2.00,
             m3, m3.nw, m3.bg, m3.00,
             m4, m4.nw, m4.bg, m4.00)
# ----------------------------------- #


# ----------------------------------- #
# Standard errors and pvalues
# ----------------------------------- #
se_pvs <- sapply(mods, function(x){
  if(nrow(model.matrix(x)) == 936){
    dat <- silver
  } else{
    dat <- silver2
  }
  return(clse_pvs(model   = x,
                  data    = dat,
                  cluster = "districtid"))
  },
  simplify = F)
# ----------------------------------- #

# ----------------------------------- #
# Models
# ----------------------------------- #
# my_note <- paste0("%stars\\\\\n",
#                   "(a) - Silverman's model\\\\\n",
#                   "(b) - Excludes population analytic weights\\\\\n",
#                   "(c) - Includes Baghdad observations (retains weights)\\\\\n",
#                   "(d) - Excludes weights and includes Baghdad." )

my_note <- paste0("\n\\item %stars.",
                  "\\item (a) - Silverman's model",
                  "\\item (b) - Excludes population analytic weights",
                  "\\item (c) - Includes Baghdad observations (retains weights)",
                  "\\item (d) - Excludes weights and includes Baghdad.")

screenreg(l = mods,
          custom.coef.map    = vars,
          custom.header      = list("Model 1" = 1:4,
                                    "Model 2" = 5:8,
                                    "Model 3" = 9:12,
                                    "Model 4" = 13:16),
          custom.model.names = rep(sprintf("(%s)",letters[1:4]),4),
          override.se        = map(se_pvs, "clse"),
          override.pvalues   = map(se_pvs, "pvals"),
          custom.note        = my_note)
# ----------------------------------- #


# ----------------------------------- #
# Save for pdf table
# ----------------------------------- #
save(list = c("mods", "vars", "se_pvs", "my_note",
              "m4.no.unit","m4.no.time","m4.no.suvs","m4.no.fes"),
     file = "scripts/silverman-replication/silverman-writeup/silverman-models.Rdata")
# ----------------------------------- #
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# SAVE                                                                    ----
#-----------------------------------------------------------------------------#
save(list = c("silver", "f", "vars"),
     file = "scripts/silverman-replication/silverman_halfyr.Rdata")
#rm(list = ls())
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# OTHER CODE & NOTES                                                      ----
#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Notes
# ----------------------------------- #
# On why we need plm() to perform "between" unit regression ala 'areg' in stata
# https://www.statalist.org/forums/forum/general-stata-discussion/general/1475427-differences-between-reg-and-areg

# Information on these analytic weights:
# https://www.stata.com/support/faqs/statistics/analytical-weights-with-linear-regression/
# https://riptutorial.com/r/example/16441/weighting
# ----------------------------------- #


# ----------------------------------- #
# Good resource on why to cluster standard errors:
# https://www.r-bloggers.com/2021/05/clustered-standard-errors-with-r/
# ----------------------------------- #
# mclse <- function(model, cluster_name, data){
#   X <- model.matrix(model)
#   u <- residuals(model)
#   n <- nobs(model)
#   p <- ncol(X)
#
#   # sigma2 <- sum(u^2) / (n - p)
#   crossXinv <- solve(t(X) %*% X, diag(p))
#
#   omegaj <- lapply(levels(data[[cluster_name]]), function(cluster) {
#     j <- data[[cluster_name]] == cluster
#     # drop = FALSE: don't drop dimensions when we have only one obs.
#     X_j <- X[j, , drop = FALSE]
#     # tcrossprod is outer product x * x^T
#     t(X_j) %*% tcrossprod(u[j]) %*% X_j
#   })
#
#   n_cl <- length(levels(data[[cluster_name]]))
#   omega <- (n-1) / (n-p) * (n_cl / (n_cl-1)) * Reduce('+', omegaj)
#   # sandwich formula; extract diagonal and take square root to get SEs
#   m1clse <- sqrt(diag(crossXinv %*% omega %*% crossXinv))
#   return(m1clse)
# }
# ----------------------------------- #


# ----------------------------------- #
# PLM vs LM
# ----------------------------------- #
# Note: the following two models are essentially identical:
# library(plm)
# d       <- pdata.frame(silver2, index = c('districtid', 'halfyr'))
# m10 <- plm(formula = f,
#            data     = d,
#            weights  = d$POP,
#            model    = 'within')
#
# m11 <- lm(update(f, ~ . -1 + district), data = silver2, weights = silver2$POP)
# m11.se <- vcovCL(m11, cluster = ~districtid)
# m11.se <- sqrt(diag(m11.se))
#
# m11.bs <- coef(m11)
# m11.bs <- m11.bs[names(m11.bs) %in% names(m11.se)]
# m11.pv <- 2 * (1 - pnorm(abs(m11.bs / m11.se)))
# screenreg(list(m10,m11), custom.coef.map = vars)
# ----------------------------------- #

#-----------------------------------------------------------------------------#
