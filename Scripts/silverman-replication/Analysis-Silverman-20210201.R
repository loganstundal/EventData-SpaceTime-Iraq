#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          February 01, 2021
# Purpose:       Silverman (IO 2020) -- Replication Script
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
# ADMINISTRATIVE --------------------------------------------------------------

#---------------------------#
# Clear working environment
#---------------------------#
rm(list = ls())

#---------------------------#
# Load required packages
#---------------------------#
library(tidyverse)

#---------------------------#
# Set working directory
#---------------------------#
# setwd("c:/users/logan/googledrive/umn/research/ra_john/2021 - Time Series - RINLA")

#---------------------------#
# Load data
#---------------------------#
silver <- readstata13::read.dta13("Scripts/silverman-replication/silverman-replication-data/CP_Data_Final_Half-Year -- 3-29-20.dta")
d      <- read_csv("Data/IRQmonthly.csv")


# Need to use original "monthly.csv" from Ben because 3-cities are aggregaed in
# my tidied data resulting in 30 fewer observations.
# load("Data/Data-All-Tidy-20210131.Rdata")
# rm(d_ts, d_cs, d_pn_yr)

#-----------------------------------------------------------------------------#
# DATA TIDY -------------------------------------------------------------------

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

  # Final variable select
  select(district, districtid, halfyr,
         ends_with("_d"), cmoc, dis_usprt, POP,
         starts_with(c("halfyr", "su_vh")))


# ----------------------------------- #


"Task 2 - Aggregate our weekly data to half year and calculate first diffs."
d2 <- d %>%
  # Rectify names
  rename(district  = district_name,
         halfyrid  = half) %>%
  select(district, halfyrid,
         p_S1, p_spentcptotal, p_spentruzicka,
         coalitioncc, insurgentcc,
         p_spentcerpsmall_noncp, p_spentusaid_nonruzicka,
         a_of_batt, cmoc, dis_usprt, POP,
         starts_with(c("halfyr", "su_vh"))) %>%

  # Aggregate from monthly to half-yearly
  group_by(district, halfyrid) %>%

  summarize(across(starts_with(c("p_","a_of_batt")), mean),
            across(starts_with(c("cmoc","dis_usprt","halfyr","su_vh")), max),
            POP = POP[1],
            .groups = "keep") %>%
  ungroup() %>%

  # Drop last half-year (Silverman does not include halfyrid == 98)
  filter(halfyrid != 98) %>%

  # At this stage the monthly data are now equivalent to Silverman's
  # replication dta file

  # Generate appropriate first differences and drop NAs on DV
  group_by(district) %>%
  mutate(across(!starts_with(c("halfyr", "su_vh", "POP")),
                .fns = list(d = ~ . - lag(., 1)))) %>%
  tidyr::drop_na(p_S1_d) %>%

  # Final variable select
  select(district, halfyrid,
         ends_with("_d"), POP,
         starts_with(c("halfyr", "su_vh")))


# ----------------------------------- #

# Compare first 10 rows of data to ensure equality

vrs <- c("p_S1_d", "p_spentcptotal_d","p_spentruzicka_d")
# vrs <- c("a_of_batt_d","cmoc_d","dis_usprt_d")

cbind(d2[1:10,vrs], silver[1:10,vrs])




#-----------------------------------------------------------------------------#




"Attempting the weekly since I have no fucking idea how he aggregated these data..."

# d.p_S1 l.d.p_spentcptotal l.d.p_spentruzicka
#        l.d.coalitioncc    l.d.insurgentcc
#        l.d.p_spentcerpsmall_noncp l.d.p_spentusaid_nonruzicka
#        l.d.a_of_batt              l.dis_usprt l.cmoc
#        halfyr1-halfyr11 su_vh1-su_vh11 [aw=POP] if districtid~=69,
#        cluster(districtid) absorb(districtid)

# Going to try with plm
library(plm)

rm(list=ls())

silver <- readstata13::read.dta13("Data/Silverman-Replication-Data/CP_Data_Final_Half-Year -- 3-29-20.dta")

f <- as.formula(paste("diff(p_S1) ~ lag(diff(p_spentcptotal)) + lag(diff(p_spentruzicka)) +
           lag(diff(coalitioncc)) + lag(diff(insurgentcc)) + lag(diff(p_spentcerpsmall_noncp)) +
           lag(diff(p_spentusaid_nonruzicka)) + lag(diff(a_of_batt)) + lag(dis_usprt) + lag(cmoc) + ",
                      "as.factor(districtid) + ",
                      paste(sprintf("halfyr%s",1:11), collapse = " + "), " + ",
                      paste(sprintf("su_vh%s", 1:11), collapse = " + ")))



# silver <- silver %>% filter(district != "Karkh")

d <- silver
d2 <- d %>% filter(district != "Karkh")


d2 <- d %>%
  filter(district != "Karkh") %>%
  # group_by(district) %>%
  mutate(d_y  = p_S1 - dplyr::lag(p_S1, 1),
         d_x1 = p_spentcptotal - dplyr::lag(p_spentcptotal,1),
         d_x2 = p_spentruzicka - dplyr::lag(p_spentruzicka, 1)) %>%
  ungroup()

f <- as.formula(paste("d_y ~ d_x1 + d_x2 + ",
                      paste(sprintf("halfyr%s",2:10), collapse = " + "), " + ",
                      paste(sprintf("su_vh%s", 2:10), collapse = " + ")))

lm(formula = f, data = d2)

d <- pdata.frame(d, index = c("districtid","month"))
d2 <- pdata.frame(d2, index = c("districtid","month"))

m <- plm(formula = f,
         data    = d,
         weights  = d$POP,
         model   = "pooling")
nm = as.list(names(coef(m))[2:10])
names(nm) = names(coef(m))[2:10]


library(texreg)
screenreg(m, digits = 2, custom.coef.map = nm)

#-----------------------------------------------------------------------------#

# ----------------------------------- #
# Just going to try the model with districtid FEs
# ----------------------------------- #
library(texreg)

f0 <- d_p_S1 ~ d_p_spentcptotal + d_p_spentruzicka + as_factor(districtid)
f1 <- as.formula( paste("p_S1_d ~ p_spentcptotal_d +p_spentruzicka_d +",
                        "as.factor(districtid) + ",
                        # "coalitioncc_d + insurgentcc_d +",
                        paste(sprintf("halfyr%s",2:10), collapse = " + "), " + ",
                        paste(sprintf("su_vh%s",2:10), collapse = " + ")))

f2 <- as.formula( paste("d_p_S1 ~ d_p_spentcptotal +d_p_spentruzicka +",
                        "as.factor(districtid) + ",
                        # "coalitioncc_d + insurgentcc_d +",
                        paste(sprintf("halfyr%s_d",2:10), collapse = " + "), " + ",
                        paste(sprintf("su_vh%s",2:10), collapse = " + ")))

f3 <- as.formula( paste("d_p_S1 ~ d_p_spentcptotal +d_p_spentruzicka +",
                        "as.factor(districtid) + ",
                        # "coalitioncc_d + insurgentcc_d +",
                        paste(sprintf("halfyr%s",2:10), collapse = " + "), " + ",
                        paste(sprintf("su_vh%s_d",2:10), collapse = " + ")))

f4 <- as.formula( paste("d_p_S1 ~ d_p_spentcptotal +d_p_spentruzicka +",
                        "as.factor(districtid) + ",
                        # "coalitioncc_d + insurgentcc_d +",
                        paste(sprintf("halfyr%s_d",2:10), collapse = " + "), " + ",
                        paste(sprintf("su_vh%s_d",2:10), collapse = " + ")))

w <- d_diff %>% mutate(w = case_when(district == "Karkh" ~ 1,
                                     TRUE ~ POP)) %>% pull(w)

m <- lm(formula = f1,
        weights = w[w!=1],
        data    = d_diff %>% filter(district != "Karkh"))


"List of what does not work:
1. Model with just relevant variables + district FEs,
2. Above + undifferenced halfyr and su_vs vars,
3. 2 with halfyr only diffed,
4. 2 with su_vh only diffed,
5. 2 with both diffed,
6. repeated 1-5 with weights vector w"

screenreg(m,
          custom.coef.map = list("d_p_spentcptotal" = "d_p_spentcptotal",
                                 "d_p_spentruzicka" = "d_p_spentruzicka",
                                 "d_coalitioncc"    = "d_coalitioncc",
                                 "d_insurgentcc"    = "d_insurgentcc",
                                 "d_p_spentcptotal:d_coalitioncc" = "spendCPxCDcoal",
                                 "d_p_spentcptotal:d_insurgentcc" = "spendCPxCDinsur"))

"True values: cptotal = ~ -0.0640119; ruzicka = ~ -0.6340512"
round(coefficients(m)[2:3],8)
# IDENTICAL

#-----------------------------------------------------------------------------#


"I think the author's theory is actually:

D.SIGACTS ~ D.COMPENSATION * D.COLLATERAL_DAMAGE,

That is: compensation dampens (i.e., moderates) the effect of collateral damage
on insurgent attacks.

But the author did not test this."


f5 <- update(f4, ~ . + d_coalitioncc + d_insurgentcc)

m <- lm(formula = f5,
        weights = w[w!=1],
        data = d_diff  %>% filter(district != "Karkh"))

screenreg(m,
          custom.coef.map = list("d_p_spentcptotal" = "d_p_spentcptotal",
                                 "d_p_spentruzicka" = "d_p_spentruzicka",
                                 "d_coalitioncc"    = "d_coalitioncc",
                                 "d_insurgentcc"    = "d_insurgentcc",
                                 "d_p_spentcptotal:d_coalitioncc" = "spendCPxCDcoal",
                                 "d_p_spentcptotal:d_insurgentcc" = "spendCPxCDinsur"))

" ^^^ Identical to T2M2"

f6 <- update(f4, ~ . + d_p_spentcptotal * d_coalitioncc + d_p_spentcptotal * d_insurgentcc)



m <- lm(formula = f6,
        weights = w[w!=1],
        data = d_diff  %>% filter(district != "Karkh"))

screenreg(m,
          custom.coef.map = list("d_p_spentcptotal" = "d_p_spentcptotal",
                                 "d_p_spentruzicka" = "d_p_spentruzicka",
                                 "d_coalitioncc"    = "d_coalitioncc",
                                 "d_insurgentcc"    = "d_insurgentcc",
                                 "d_p_spentcptotal:d_coalitioncc" = "spendCPxCDcoal",
                                 "d_p_spentcptotal:d_insurgentcc" = "spendCPxCDinsur"))
# SOME SUPPORT FOR MODERATION HYPOTHESES - NEED TO CLUSTER ERRORS TO DOUBLE-CHECK.



# All the fes really put a strain on this estimation so temp form:
f7 <- as.formula( paste("d_p_S1 ~ d_p_spentcptotal +d_p_spentruzicka +",
                        "d_p_spentcptotal * d_coalitioncc + d_p_spentcptotal * d_insurgentcc"))

dd<- d_diff  %>% filter(district != "Karkh") # margins can't use tidyverse internally
m <- lm(formula = f7,
        weights = w[w!=1],
        data = dd)


library(margins)

mm <- margins(m, variables = "d_p_spentcptotal")
summary(mm)


cplot(m, what = "effect", x = "d_p_spentcptotal", dx = "d_coalitioncc")


#-----------------------------------------------------------------------------#

















#-----------------------------------------------------------------------------#
# SAVE ------------------------------------------------------------------------
# save.image()
# rm(list = ls())



#-----------------------------------------------------------------------------#
# My notes attempting to figure out how to replicate Silverman
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# What is going on with Stata's "agre" - Silverman do file?
#-----------------------------------------------------------------------------#
# https://www.stata.com/manuals/rareg.pdf (p. 3)
"If you have more variables than can fit in a Stata matrix, 'regress' will not work.
'areg' provides a way of obtaining estimates of beta - but not the (dummy parameter)
gamma[i's]- in these cases. The effects of the dummy variables are said to be 'absorbed"

# So it is quite literally a regression with more dummy variables than stata's very
# limited can handle...

# This document:
# http://www.princeton.edu/~otorres/Panel101.pdf (p. 22)
# further validates this - "areg" is just reg automating the FEs


# ----------------------------------- #
# Silverman's Table 2, Model 1 fn call:
# ----------------------------------- #
# eststo: quietly areg d.p_S1 d.p_spentcptotal d.p_spentruzicka
#                      halfyr2-halfyr10
#                      su_vh2-su_vh10
#                      [aw=POP] if district~="Karkh",
#                      cluster(districtid) absorb(districtid)

# Breakdown:
"eststo - just stores model estimates for later formatting"
"quietly - just suppresses the output on the console"
"areg - command as specified above for reg with more dummies than stata can
        handle in a matrix"

"Then the inclusion of EXACTLY the variables called with 'd ~ data' above with
exception of the d. operator - more in a bit.

cluster() and absorb() are stata default se clustering and 'absorption' of the
excess dummies per the areg command.

So any problems I'm encountering are on the bahavior of 'd.'!"



#-----------------------------------------------------------------------------#
# So, what is going on with Stata's "d." operator?
#-----------------------------------------------------------------------------#

"nothing... it's just a first differences operator. Presumably when this dta is
loaded into stata it recalls xtset."

#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Weights
#-----------------------------------------------------------------------------#

"[aw=POP] if district~='Karkh'"

"This assigns analytic weights (sigma/omega ~ omega = weights) to population
 for all observations where district does NOT equal Karkh.

This choice DEFINITELY impacts the results flipping coefficient signs and
significance which is really surprising since it is NOT discussed anywhere
in the article or appendix and is entirely unmotivated theoretically."


"Silverman's models have 927 observations. There are no missing values remaining
in d meaning he drops a variable whithout saying how/why. likely related to
[aw=POP] if district~='Karkh'"
# table(d_diff$district == "Karkh")

#-----------------------------------------------------------------------------#






















