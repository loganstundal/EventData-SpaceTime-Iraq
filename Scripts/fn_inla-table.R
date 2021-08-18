#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          August 06, 2021
# Purpose:       functions to format INLA model parameters into tables
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# inla_params() ~ Extract parameters of interest from a fit inla model object
#-----------------------------------------------------------------------------#
inla_params <- function(model_list, spde){
  # ----------------------------------- #
  # function description
  # ----------------------------------- #
  # qoi() [quantity of interest] takes an INLA model list and returns a list
  # containing quantities of interest: posterior parameter estimates as well as
  # credibility intervals on included regressors and parameters for the GMRF
  # (range, kappa, and sigma).
  # To extract these quantities for one model, input your model as a named list:
  # result <- qoi(model_list = list("my_model" = estimated_model))

  # ----------------------------------- #

  params <- lapply(model_list, function(mod){
    # ----------------------------------- #
    # Fixed effects
    # ----------------------------------- #
    betas <- mod[["summary.fixed"]][,c("0.5quant",
                                       "0.025quant","0.975quant")]

    betas <- betas %>% rename(median = `0.5quant`,
                      lb     = `0.025quant`,
                      ub     = `0.975quant`) %>%
      mutate(variable = rownames(.), type = "fixed") %>%
      dplyr::select(variable, median, lb, ub, type)
    # ----------------------------------- #

    # ----------------------------------- #
    # Hyper-parameters
    # ----------------------------------- #
    spde_pars <- inla.spde2.result(inla = mod,
                                   name = "i",
                                   spde = spde,
                                   do.transform = TRUE)
    # Medians
    Kappa <- inla.qmarginal(0.50, spde_pars$marginals.kappa[[1]])
    Sigma <- inla.qmarginal(0.50, spde_pars$marginals.variance.nominal[[1]])
    Range <- inla.qmarginal(0.50, spde_pars$marginals.range.nominal[[1]])
    # Rho   <- inla.qmarginal(0.50, mod$marginals.hyperpar$`GroupRho for i`)

    # Compute 95% HPD intervals
    Kappahpd <- inla.hpdmarginal(0.95, spde_pars$marginals.kappa[[1]])
    Sigmahpd <- inla.hpdmarginal(0.95, spde_pars$marginals.variance.nominal[[1]])
    Rangehpd <- inla.hpdmarginal(0.95, spde_pars$marginals.range.nominal[[1]])
    # Rhohpd   <- inla.hpdmarginal(0.95, mod$marginals.hyperpar$`GroupRho for i`)

    # Convert range to km (degrees = 2*pi*6371/360)
    Range    <- Range * 2*pi*6371/360
    Rangehpd <- Rangehpd * 2*pi*6371/360

    # Collect hyper-parameters
    hyper <- rbind(c(Kappa,  Kappahpd),
                   c(Sigma,  Sigmahpd),
                   c(Range,  Rangehpd)
                   # c(Rho,    Rhohpd)
                   ) %>%
      as.data.frame() %>%
      rename(median = 1, lb = 2, ub = 3) %>%
      mutate(variable = c("Kappa","Sigma^2","Range"),
             type     = "hyper")
    # ----------------------------------- #

    # ----------------------------------- #
    # Model log likelihoods and n
    # ----------------------------------- #
    llik  <- mod$mlik[2]

    mod_n <- mod$size.linear.predictor$Ntotal -
        mod$size.linear.predictor$n

    str <- data.frame(
      variable = c("N", "LogLik"),
      median   = c(mod_n, llik),
      type     = "str"
    )
    # ----------------------------------- #

    # Return tidy model results
    res <- bind_rows(betas, hyper, str)
    rownames(res) <- 1:nrow(res)
    return(res)
  })

  # ----------------------------------- #
  # fn Return statement
  # ----------------------------------- #
  return(params)
  # ----------------------------------- #
}
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
# inla_table() ~ print a table using data returned from inla_params()
#-----------------------------------------------------------------------------#
inla_table <- function(params,
                       caption = NULL, model_names = NULL,
                       header  = NULL,
                       ...){
  # ----------------------------------- #
  # Description
  # ----------------------------------- #
  # Takes model parameter data exported from inla_params() and produces a
  # formatted regression table.
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(kableExtra)
  # ----------------------------------- #

  # ----------------------------------- #
  # Tidy parameters for table
  # ----------------------------------- #
  tidy_params <- lapply(params, function(mod){
    mod %>%
      mutate(median = format(round(median, digits = 3), nsmall = 3),
             lb     = format(round(lb,     digits = 3), nsmall = 3),
             ub     = format(round(ub,     digits = 3), nsmall = 3)) %>%
      mutate(across(.cols = c(median, lb, ub), .fns = ~ str_trim(.x))) %>%
      mutate(hpd    = sprintf("[%s, %s]", lb, ub)) %>%
      pivot_longer(.,
                   cols = c(median,hpd),
                   values_to = "val",
                   names_to = "stat") %>%
      filter(!(type == "str" & stat == "hpd")) %>%
      dplyr::select(variable, val, stat, type)
  }) %>%
    bind_rows(., .id = "model") %>%
    mutate(variable = case_when(variable %in%
                                  c("sigacts_l","icews_l","ged_l") ~ "Phi",
                                TRUE ~ variable)) %>%
    pivot_wider(.,
                id_cols     = c(variable, stat, type),
                names_from  = model,
                values_from = val) %>%
    mutate(type = factor(type, levels = c("fixed","hyper","str")))

  tidy_params <- tidy_params %>%
    group_by(type) %>%
    mutate(z = 1:n()) %>%
    ungroup %>%
    mutate(z = case_when(variable == "Phi" ~ as.integer(100),
                         variable == "Intercept" ~ as.integer(101),
                         TRUE ~ z)) %>%
    arrange(type, z)

  # TO DO - RENAME WITH LIST, DROP WITH VECTOR

  par_lengths <- table(tidy_params$type)

  tidy_params <- tidy_params %>%
    dplyr::select(-stat, -type, -z) %>%
    mutate(variable = str_remove_all(variable, "_")) %>%
    mutate(across(.cols = everything(),
                  .fns  = ~replace_na(.x, ""))) %>%
    mutate(across(.cols = !variable,
                  ~cell_spec(.x, font_size = ifelse(stringr::str_detect(.x, "\\["), 8, 10))))
  # ----------------------------------- #


  # ----------------------------------- #
  # Create table
  # ----------------------------------- #
  # Number of models
  md_ct <- length(params)

  # Model names
  if(is.null(model_names)){
    model_names <- c("", sprintf("(%s)", 1:md_ct))
  }

  kbl_table <-  kbl(x         = tidy_params,
                    caption   = caption,
                    escape    = FALSE,
                    col.names = model_names,
                    booktabs  = TRUE,
                    align     = c("l",rep("c", md_ct))) %>%
    kable_styling(full_width = TRUE,
                  font_size  = 10) %>%
    collapse_rows(columns     = 1,
                  valign      = "middle",
                  latex_hline = "none") %>%
    row_spec(as.numeric(par_lengths["fixed"]),
             extra_latex_after = "\\midrule") %>%
    row_spec(as.numeric(par_lengths["fixed"] + par_lengths["hyper"]),
             extra_latex_after = "\\midrule")
    # footnote(general = footnote, footnote_as_chunk = T)

  if(!is.null(header)){
    kbl_table <- kbl_table %>%
      add_header_above(header)
  }
  # ----------------------------------- #

  # ----------------------------------- #
  # Function return
  # ----------------------------------- #
  return(kbl_table)
  # ----------------------------------- #
}
#-----------------------------------------------------------------------------#
