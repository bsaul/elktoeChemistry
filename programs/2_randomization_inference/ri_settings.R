#-----------------------------------------------------------------------------#
#   Title: Configuration for analysis
#  Author: B Saul
# Purpose: 
#-----------------------------------------------------------------------------#

NSIMS = 1000
# .elements <- c("Pb_ppm_m209", "U_ppm_m238")
 .elements <- "all"
RI_CONFIG <- list(
  
  # # A : Are sites (including baseline) comparable in past?
  # 
  # A_gam = list(
  #   label     = "A",
  #   sublabel  = "ri[gam]",
  #   desc      = "Are sites (including baseline) comparable in past?",
  #   filters   = list(
  #     contrast  = "all",
  #     agrp      = "agrp_first_transect_with_AB",
  #     elements  = .elements,
  #     signals   = c("base", "avg5_trunc_3sd"),
  #     group_by_valve = FALSE,
  #     transect_opts  = list(
  #       .layers    = "ncr",
  #       .annuli    = LETTERS[c(2:14)], # B-N
  #       .min_n_obs = 10L
  #     )
  #   ),
  #   # RI functions
  #   prep_FUN      = prep_for_gam_ri,
  #   dec_FUN       = define_multiarm_cluster_declaration,
  #   test_data_FUN = identity,
  #   test_statistic_FUN  = make_gam_ts(
  #     m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + factor(drawer) +  s(analysis_id, bs = "re"),
  #     m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + factor(drawer) + s(analysis_id, bs = "re")
  #   ),
  #   ri_FUN = do_ri_gam,
  #   nsims      = NSIMS
  # ),
  A_wls = list(
    label     = "A",
    sublabel  = "ri[wls]",
    desc      = "Are experimental specimens similar to baseline specimens in past?",
    filters   = list(
      contrast  = "baseline_v_sites",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[2:14], # A-N
        .min_n_obs = 5L
      )
    ),
    # RI functions
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = create_wls_data,
    test_statistic_FUN  = make_wls_ts(
      m1_rhs = ~ Z + factor(drawer) + I(baseline_volume/1000),
      m2_rhs = ~     factor(drawer) + I(baseline_volume/1000)
    ),
    ri_FUN = do_ri_gam,
    nsims  = NSIMS
  ),

  # # B: Are sites (excluding baseline) comparable in past?
  # 
  # B_gam = list(
  #   label     = "B",
  #   sublabel  = "ri[gam]",
  #   desc      = "Are sites (excluding baseline) comparable in past?",
  #   filters   = list(
  #     contrast  = "nobaseline",
  #     agrp      = "agrp_first_transect_with_AB",
  #     elements  = .elements,
  #     signals   = c("base", "avg5_trunc_3sd"),
  #     group_by_valve = FALSE,
  #     transect_opts  = list(
  #       .layers    = "ncr",
  #       .annuli    = LETTERS[c(2:14)], # B-N
  #       .min_n_obs = 10L
  #     )
  #   ),
  #   # RI functions
  #   prep_FUN      = prep_for_gam_ri,
  #   dec_FUN       = define_multiarm_cluster_declaration,
  #   test_data_FUN = identity,
  #   test_statistic_FUN  = make_gam_ts(
  #     m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + factor(drawer) + baseline_volume + s(analysis_id, bs = "re"),
  #     m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + factor(drawer) + baseline_volume + s(analysis_id, bs = "re")
  #   ),
  #   ri_FUN = do_ri_gam,
  #   nsims      = NSIMS
  # ),
  B_wls = list(
    label     = "B",
    sublabel  = "ri[wls]",
    desc      = "Are sites (excluding baseline) comparable in past?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[2:14], # A-N
        .min_n_obs = 5L
      )
    ),
    # RI functions
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = create_wls_data,
    test_statistic_FUN  = make_wls_ts(
      m1_rhs = ~ Z + factor(drawer) + I(baseline_volume/1000),
      m2_rhs = ~     factor(drawer) + I(baseline_volume/1000)
    ),
    ri_FUN = do_ri_gam,
    nsims  = NSIMS
  ),

  # # C: Are sites (excluding baseline) different?
  
  #####
  # C_gam = list(
  #   label     = "C",
  #   sublabel  = "ri[gam]",
  #   desc      = "Are sites (excluding baseline) different?",
  #   filters   = list(
  #     contrast  = "nobaseline",
  #     agrp      = "agrp_transect_most_A",
  #     elements  = .elements,
  #     signals   = c("base", "avg5_trunc_3sd"),
  #     group_by_valve = FALSE,
  #     transect_opts  = list(
  #       .layers    = "ncr",
  #       .annuli    = LETTERS[c(1:14)], # A-N
  #       .min_n_obs = 10L
  #     )
  #   ),
  #   # RI functions
  #   prep_FUN      = prep_for_gam_ri,
  #   dec_FUN       = define_multiarm_cluster_declaration,
  #   test_data_FUN = identity,
  #   test_statistic_FUN  = make_gam_ts(
  #     m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + factor(drawer) + baseline_volume + s(analysis_id, bs = "re"),
  #     m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + factor(drawer) + baseline_volume + s(analysis_id, bs = "re")
  #   ),
  #   ri_FUN = do_ri_gam,
  #   nsims      = NSIMS
  # )
  
  C_wls = list(
    label     = "C",
    sublabel  = "ri[wls]",
    desc      = "Are sites (excluding baseline) different?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[1], # A-N
        .min_n_obs = 5L
      )
    ),
    # RI functions
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = create_wls_data,
    test_statistic_FUN  = make_wls_ts(
      m1_rhs = ~ Z + factor(drawer) + I(baseline_volume/1000),
      m2_rhs = ~     factor(drawer) + I(baseline_volume/1000)
    ),
    ri_FUN = do_ri_gam,
    nsims  = NSIMS
  ),
  
  
  D_wls = list(
    label     = "D",
    sublabel  = "ri[wls]",
    desc      = "Are sites (excluding baseline) different in pio?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "pio",
        .min_n_obs = 5L
      )
    ),
    # RI functions
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = create_wls_data,
    test_statistic_FUN  = make_wls_ts(
      m1_rhs = ~ Z + factor(drawer) + I(baseline_volume/1000),
      m2_rhs = ~     factor(drawer) + I(baseline_volume/1000)
    ),
    ri_FUN = do_ri_gam,
    nsims  = NSIMS
  ),
    
  E_wls = list(
    label     = "E",
    sublabel  = "ri[wls]",
    desc      = "Are sites (excluding baseline) different in psm?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "psm",
        .min_n_obs = 5L
      )
    ),
    # RI functions
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = create_wls_data,
    test_statistic_FUN  = make_wls_ts(
      m1_rhs = ~ Z + factor(drawer) + I(baseline_volume/1000),
      m2_rhs = ~     factor(drawer) + I(baseline_volume/1000)
    ),
    ri_FUN = do_ri_gam,
    nsims  = NSIMS
  )
)
