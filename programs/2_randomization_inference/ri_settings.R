#-----------------------------------------------------------------------------#
#   Title: Configuration for analysis
#  Author: B Saul
# Purpose: 
#-----------------------------------------------------------------------------#

NSIMS = 1000
.elements <- "all"

RI_CONFIG <- list(
  
  # A : Are sites (including baseline) comparable in past?
  
  A_gam = list(
    label     = "A",
    sublabel  = "ri[gam]",
    desc      = "Are sites (including baseline) comparable in past?",
    filters   = list(
      contrast  = "all",
      agrp      = "agrp_first_transect_with_AB",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = FALSE,
      transect_opts  = list(
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
      m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    ri_FUN = do_ri_gam,
    nsims      = NSIMS
  ),
  
  A_mom = list(
    label     = "A",
    sublabel  = "ri[mom]",
    desc      = "Are sites (including baseline) comparable in past?",
    filters   = list(
      contrast  = "all",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[2:14], # B-N
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_A_mom,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN = do_ri_summary_stats,
    nsims      = NSIMS
  ),
  
  A_momd = list(
    label     = "A",
    sublabel  = "ri[momd]",
    desc      = "Are sites (including baseline) comparable in past?",
    filters   = list(
      contrast  = "all",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[2:14], # B-N
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_A_mom_diff,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN      = do_ri_summary_stats,
    nsims       = NSIMS
  ),
  
  
  # B: Are sites (excluding baseline) comparable in past?
  
  B_gam = list(
    label     = "B",
    sublabel  = "ri[gam]",
    desc      = "Are sites (excluding baseline) comparable in past?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "agrp_first_transect_with_AB",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = FALSE,
      transect_opts  = list(
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
      m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    ri_FUN = do_ri_gam,
    nsims      = NSIMS
  ),
  
  B_mom = list(
    label     = "B",
    sublabel  = "ri[mom]",
    desc      = "Are sites (excluding baseline) comparable in past?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[2:14], # B-N
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_A_mom,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN = do_ri_summary_stats,
    nsims      = NSIMS
  ),
  
  B_momd = list(
    label     = "B",
    sublabel  = "ri[momd]",
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
    test_data_FUN = create_summary_stats_data_A_mom_diff,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN      = do_ri_summary_stats,
    nsims       = NSIMS
  ),
 
  
  # C: Are sites (excluding baseline) different?
  
  C_gam = list(
    label     = "C",
    sublabel  = "ri[gam]",
    desc      = "Are sites (excluding baseline) different?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "agrp_first_transect_with_A",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = FALSE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[c(1:14)], # B-N
        .min_n_obs = 10L
      )
    ),
    # RI functions
    prep_FUN      = prep_for_gam_ri,
    dec_FUN       = define_multiarm_cluster_declaration,
    test_data_FUN = identity,
    test_statistic_FUN  = make_gam_ts(
      m1_rhs = ~ s(d, bs = "ts") + d*Z + pd*annuli*Z + s(pd, bs = "ts") + s(analysis_id, bs = "re"),
      m2_rhs = ~ s(d, bs = "ts") + d   + pd*annuli   + s(pd, bs = "ts") + s(analysis_id, bs = "re")
    ),
    ri_FUN = do_ri_gam,
    nsims      = NSIMS
  ),
  
  C_mom = list(
    label     = "C",
    sublabel  = "ri[mom]",
    desc      = "Are sites (excluding baseline) different?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[1:14], 
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_A_mom,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN = do_ri_summary_stats,
    nsims      = NSIMS
  ),
  
  C_momd = list(
    label     = "C",
    sublabel  = "ri[momd]",
    desc      = "Are sites (excluding baseline) different?",
    filters   = list(
      contrast  = "nobaseline",
      agrp      = "all",
      elements  = .elements,
      signals   = c("base", "avg5_trunc_3sd"),
      group_by_valve = TRUE,
      transect_opts  = list(
        .layers    = "ncr",
        .annuli    = LETTERS[1:14], # A-N
        .min_n_obs = 5L
      )
    ),
    
    # RI functions
    test_data_FUN = create_summary_stats_data_C_mom_diff,
    prep_FUN      = replace_na_median,
    dec_FUN       = define_multiarm_declaration,
    test_statistic_FUN = kw_test_fun,
    ri_FUN      = do_ri_summary_stats,
    nsims       = NSIMS
  ) 
)
