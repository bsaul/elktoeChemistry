NSIMS = 5

RI_CONFIG <- list(
  
  A_gam = list(
    label = "A",
    sublabel = "ri[gam]",
    desc  = "Are sites (including baseline) comparable in past?",
    filters   = list(
      contrast  = "all",
      agrp      = "agrp_first_transect_with_AB",
      group_by_valve = FALSE,
      # elements  = "all",
      elements  = c("As_ppm_m75", "Co_ppm_m59"),
      signals   = c("base", "avg5_trunc_3sd"),
      transect_opts = list(
        .layers    = "ncr",
        .annuli    = LETTERS[c(2:14)], # B-N
        .min_n_obs = 10L
      )
    ),
    # RI functions
    prep_FUN      = prep_for_gam_ri,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = identity,
    test_statistic_FUN  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    ri_FUN = do_ri_gam,
    nsims      = NSIMS
  ),
  
  A_mom = list(
    label = "A",
    sublabel = "ri[mom]",
    desc  = "Are sites (including baseline) comparable in past?",
    filters   = list(
      contrast  = "all",
      group_by_valve = TRUE,
      agrp      = "all",
      # elements  = "all",
      elements  = c("As_ppm_m75", "Co_ppm_m59"),
      signals   = c("base", "avg5_trunc_3sd"),
      transect_opts = list(
        .layers    = "ncr",
        .annuli    = LETTERS[c(2:14)], # B-N
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_A_mom,
    prep_FUN      = identity,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN = do_ri_summary_stats,
    nsims      = NSIMS
  )
  
)
