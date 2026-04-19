
plot_traces <- function(model, pars, label_mapping = NULL) {
  df_chains <- extract_chains(model, pars)
  if (!is.null(label_mapping)) {
    df_chains <- (
      df_chains
      %>% mutate(
        parameter = recode(
          parameter,
          !!!label_mapping
        )
      )
    )
  }
  (
    df_chains
    %>% ggplot()
    + geom_line(
      aes(x = step, y = value, color = chain), 
      alpha = 1
    )
    + facet_wrap(
      ~ parameter,
      ncol = 3,
      scales = "free",
      labeller = label_parsed
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      plot.background = element_blank(),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 18),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 15),
      axis.ticks.length = unit(4, "pt"),
      axis.ticks = element_line(linewidth = 1),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15)
    )
    +scale_color_brewer(
      name = "Chains",
      palette = "Reds",
      direction = -1
    )
    + scale_size(range=c(1, 2), guide=FALSE)
    + scale_y_continuous(
      breaks = pretty_breaks(n = 2)
    )
    + scale_x_continuous(
      breaks = pretty_breaks(n = 3)
    )
    + labs(
      x = "MCMC Iteration",
      y = "Value"
    )
    + guides(
      color = guide_legend(override.aes = list(linewidth = 4)),
    )
  )
}

plot_post_pairs <- function(model, pars, label_mapping = NULL) {
  df_post <- extract_post(model, pars)
  if (!is.null(label_mapping)) {
    names(df_post) <- label_mapping[names(df_post)]
  }
  (
    ggpairs(
      df_post,
      lower = list(continuous = wrap("points", alpha = 0.1, pch=21, fill = col_palette[6]), color = col_palette[5]),
      upper = list(continuous = wrap("density", color = col_palette[5])),
      diag = list(continuous = wrap("barDiag", fill = col_palette[5], alpha = 0.7, bins = 20)),
      progress = FALSE,
      labeller = label_parsed
    )
    + theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      axis.ticks.length = unit(4, "pt"),
      axis.ticks = element_line(linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 14),
    )
    + scale_y_continuous(
      breaks = pretty_breaks(n = 3)
    )
    + scale_x_continuous(
      breaks = pretty_breaks(n = 3)
    )
  )
}

plot_ppc_overlay <- function(fit, data, outcome_var, n_samples = 500, common_legend = TRUE, title = NA) {
  strip_theme <- theme(
    plot.title = element_text(size = 14),
    legend.position.inside = c(0.75, 0.85)
  )
  df_ppc <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      re_formula = NA
    )
    %>% ungroup()
    %>% select(
      country,
      cluster,
      .row,
      .draw,
      sv_pred = .epred,
      sv_obs = !!sym(outcome_var)
    )
    %>% group_by(country, cluster, .draw)
    %>% summarise(
      mean_sv_pred = mean(sv_pred, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
  )
  df_obs <- (
    df_observed
    %>% group_by(country, cluster)
    %>% summarise(
      mean_sv_obs = mean(!!sym(outcome_var), na.rm = TRUE)
    )
  )
  country_list <- unique(df_observed$country)
  plot_list <- list()
  for (i in seq_along(country_list)) {
    plot_list[[i]] <- (
      ggplot()
      + geom_density(
        data = df_ppc %>% filter(country == country_list[[i]]),
        mapping = aes(
          x = mean_sv_pred,
          group = .draw
        ),
        alpha = 0.2, 
        color = "#81c9db"
      )
      + geom_density(
        data = df_obs %>% filter(country == country_list[[i]]),
        mapping = aes(
          x = mean_sv_obs
        ),
        color = "#0292b7",
        linewidth = 1.5
      )
      # + scale_color_manual(
      #   name = "",
      #   values = c(
      #     "observed" = col_palette[1],
      #     "simulated" = col_palette[2]
      #   )
      # )
      # + scale_x_continuous(lim = c(-0.01, 0.2))
      + theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        plot.background = element_blank(),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 14),
        axis.ticks.length = unit(4, "pt"),
        axis.ticks = element_line(linewidth = 1),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.7),
        strip.text = element_text(size = 15),
        strip.background = element_blank(),
      )
      + strip_theme
      + guides(
        color = "none",
      )
    )
  }

  p <- (
    (plot_list[[1]]+plot_list[[2]]) /
    (plot_list[[3]]+plot_list[[4]]) /
    (plot_list[[5]]+plot_list[[6]])
    + plot_layout(
      axes = "collect",
      axis_titles = "collect",
    )
  )

  p
}

plot_psis_diagnostics <- function(loo, df_observed, title = "") {
  ks <- loo$pointwise[, "influence_pareto_k"]
  (
    tibble(
      degurba = df_observed$degurba,
      ks = ks
    )
    %>% arrange(degurba)
    %>% mutate(
      degurba = as.factor(degurba),
      n = 1:n()
    )
    %>% ggplot()
    + plot_theme
    + theme(
      legend.position = "right"
    )
    + geom_point(
      aes(
        x = n,
        y = ks,
        color = degurba
      )
    )
    + geom_hline(yintercept = 1.0, linetype = "dashed", color = "firebrick", linewidth = 1)
    + geom_hline(yintercept = 0.7, linetype = "dotted", color = "lightcoral", linewidth = 1)
    + scale_color_manual(
      labels = ur_labels,
      values = ur_colors
    )
    + labs (
      x = "Data point",
      y = "Pareto shape k",
      color = "Degree of urbanisation", 
      title = title
    )
    + guides(
      color = guide_legend(override.aes = list(size = 4))
    )
  )
}

plot_ppc_rootogram <- function(fit, data, outcome_var, cred1 = 0.5, cred2 = 0.95) {
  df_ppc <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      re_formula = NA
    )
    %>% ungroup()
    %>% select(
      country,
      adm_name,
      .row,
      .draw,
      sv_pred = .epred,
      sv_obs = !!sym(outcome_var)
    )
    %>% pivot_longer(
      cols = c(sv_pred, sv_obs),
      names_to = "group",
      values_to = "value"
    )
    %>% group_by(country, adm_name, group)
    %>% summarise(
      hdi50low = bayestestR::hdi(value, ci = cred1)$CI_low,
      hdi50high = bayestestR::hdi(value, ci = cred1)$CI_high,
      hdi95low = bayestestR::hdi(value, ci = cred2)$CI_low,
      hdi95high = bayestestR::hdi(value, ci = cred2)$CI_high,
      # q025 = quantile(value, 0.025, na.rm = TRUE),
      # q25 = quantile(value, 0.1, na.rm = TRUE),
      mean_sv = mean(value, na.rm = TRUE),
      median_sv = median(value, na.rm = TRUE),
      # q75 = quantile(value, 0.9, na.rm = TRUE),
      # q975 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
    %>% group_by(country, adm_name)
    %>% mutate(
      mean_sv_obs = mean_sv[group == "sv_obs"][1]
    )
    %>% ungroup()
    %>% mutate(
      adm_name = reorder_within(adm_name, desc(mean_sv_obs), country)
    )
  )
  p <- (
    ggplot()
    + geom_segment(
      data = df_ppc %>% filter(group == "sv_pred"),
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = hdi95low,
        # y = q025,
        yend = hdi95high
        # yend = q975
      ),
      color = "#DEEBF7FF",
      linewidth = 1
    )
    + geom_segment(
      data = df_ppc %>% filter(group == "sv_pred"),
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = hdi50low,
        # y = q25,
        yend = hdi50high
        # yend = q75
      ),
      color = "#9ECAE1FF",
      linewidth = 1
    )
    + geom_point(
      data = df_ppc %>% filter(group == "sv_pred"),
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        # y = median_sv,
        y = mean_sv,
        color = "predicted"
      ),
      size = 1.5,
    )
    + geom_point(
      data = df_ppc %>% filter(group == "sv_obs"),
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = mean_sv,
        color = "observed"
      ),
      size = 1.5,
    )
    + facet_wrap(
      ~ country,
      nrow = 1,
      scales = "free_x"
    )
    + scale_color_manual(
      labels = c(
        "observed" = "Observed",
        "predicted" = "Replicated"
      ),
      values = c(
        "observed" = "#2171B5FF",
        "predicted" = "#9ECAE1FF"
      )
    )
    + scale_x_reordered()
    + scale_y_continuous(
      labels = function(x) format(x * 100, nsmall = 0)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 13),
      axis.ticks.length = unit(4, "pt"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 1),
      axis.title = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13),
      legend.title = element_blank(),
      legend.text = element_text(size = 13, margin = margin()),
      legend.background = element_blank(),
      strip.text = element_text(size = 13),
      strip.background = element_blank()
    )
    + labs(
      x = "Administrative area",
      y = "NP violence (%)"
    )
  )
  p
}

plot_ppc_rootogram_weighted <- function(fit, data, outcome_var, cred1 = 0.5, cred2 = 0.95) {
  df_obs <- (
    data
    %>% group_by(country, adm_name)
    %>% summarise(
      mean_sv_obs := sum(!!sym(outcome_var) * sample_weight, na.rm = TRUE) / sum(sample_weight, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
    # %>% mutate(
    #   adm_name = reorder_within(adm_name, desc(mean_sv), country)
    # )
  )
  df_ppc <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      re_formula = NA
    )
    %>% ungroup()
    %>% select(
      country,
      adm_name,
      .row,
      .draw,
      sv_pred = .epred,
    )
    %>% group_by(country, adm_name)
    %>% summarise(
      hdi50low = bayestestR::hdi(sv_pred, ci = cred1)$CI_low,
      hdi50high = bayestestR::hdi(sv_pred, ci = cred1)$CI_high,
      hdi95low = bayestestR::hdi(sv_pred, ci = cred2)$CI_low,
      hdi95high = bayestestR::hdi(sv_pred, ci = cred2)$CI_high,
      mean_sv_pred = mean(sv_pred, na.rm = TRUE),
      median_sv_pred = median(sv_pred, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
    %>% merge(
      df_obs,
      by = c("country", "adm_name"),
      all.x = TRUE
    )
    # %>% mutate(
    #   adm_name = reorder_within(adm_name, desc(mean_sv_obs), country)
    # )
  )
  p <- (
    ggplot()
    + geom_segment(
      data = df_ppc,
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = hdi95low,
        yend = hdi95high
      ),
      color = "#DEEBF7FF",
      linewidth = 1
    )
    + geom_segment(
      data = df_ppc ,
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = hdi50low,
        yend = hdi50high
      ),
      color = "#9ECAE1FF",
      linewidth = 1
    )
    + geom_point(
      data = df_ppc,
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = mean_sv_pred,
        color = "predicted"
      ),
      size = 1.5,
    )
    + geom_point(
      data = df_ppc,
      mapping = aes(
        x = reorder(adm_name, desc(mean_sv_obs)),
        y = mean_sv_obs,
        color = "observed"
      ),
      size = 1.5,
    )
    + facet_wrap(
      ~ country,
      nrow = 1,
      scales = "free_x"
    )
    + scale_color_manual(
      labels = c(
        "observed" = "Observed",
        "predicted" = "Replicated"
      ),
      values = c(
        "observed" = "#2171B5FF",
        "predicted" = "#9ECAE1FF"
      )
    )
    + scale_x_reordered()
    + scale_y_continuous(
      labels = function(x) format(x * 100, nsmall = 0)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 13),
      axis.ticks.length = unit(4, "pt"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 1),
      axis.title = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13),
      legend.title = element_blank(),
      legend.text = element_text(size = 13, margin = margin()),
      legend.background = element_blank(),
      strip.text = element_text(size = 13),
      strip.background = element_blank()
    )
    + labs(
      x = "Administrative area",
      y = "NP violence (%)"
    )
  )
  p
}


plot_baseline_levels <- function(fit, lims = NULL) {
  baseline <- (
    fit
    %>% as_draws_rvars()
    %>% spread_rvars(b_Intercept, r_country[country, 1])
    %>% filter(`1` == "Intercept")
    %>% select(-c(`1`))
    %>% mutate(
      country_mean_lat = (b_Intercept + r_country),
      country_mean = 1 / (1 + exp(-country_mean_lat))
    )
  )
  p <- (
    ggplot(
      data = baseline,
      aes(
        xdist=country_mean,
        y = fct_reorder(country, country_mean, .fun = "median")
      ),
    )
    + stat_halfeye(
      color = "#A63603FF",
      fill = "#FDD0A2FF",
      adjust = 1.0,
      p_limits = c(0.1, 0.9),
      normalize = "panels",
      trim = TRUE
    )
    + scale_x_continuous(
      labels = scales::percent,
      limits = lims
    )
    + labs(
      x = "Baseline Sexual Violence (%)",
      y = "Country"
    )
    + guides(
      fill = "none",
      color = "none"
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
    )
  )
  p
}

plot_bin_marginal <- function(fit, data, covs, cred1 = 0.5, cred2 = 0.95, lims = NULL) {
  marg_diff <- create_df_marginal_bin_differences(fit, data, covs)
  marg_diff_stats <- create_df_marginal_bin_differences_stats(marg_diff, cred1 = cred1, cred2 = cred2)
  p <- (
    ggplot(
      data = marg_diff,
      mapping = aes(
        y = fct_reorder(covariate, expected_diff, .fun = ~ -median(.)),
        x = expected_diff,
      )
    )
    + ggridges::geom_density_ridges(
      aes(height = after_stat(scaled)),
      stat = "density",
      orientation = "x",
      scale = 0.6,
      rel_min_height = 0.01,
      alpha = 0.7,
      color = NA,
      fill = "#FDD0A2FF"
    )
    + geom_segment(
      data = marg_diff_stats,
      # aes(y = covariate, yend = covariate, x = lower95_exp_diff, xend = upper95_exp_diff),
      aes(y = covariate, yend = covariate, x = lower2_exp_diff, xend = upper2_exp_diff),
      linewidth = 1,
      alpha = 0.4,
      color = "#A63603FF",
      inherit.aes = FALSE
    )
    + geom_segment(
      data = marg_diff_stats,
      # aes(y = covariate, yend = covariate, x = lower95_exp_diff, xend = upper95_exp_diff),
      aes(y = covariate, yend = covariate, x = lower1_exp_diff, xend = upper1_exp_diff),
      linewidth = 1.5,
      alpha = 0.7,
      color = "#A63603FF",
      inherit.aes = FALSE
    )
    + geom_point(
      data = marg_diff_stats,
      aes(y = covariate, x = median_exp_diff),
      color = "#A63603FF",
      size = 3,
      inherit.aes = FALSE
    )
    + geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
    )
    + labs(
      x = "Average marginal effect\n(p(violence|1) − p(violence|0))",
    )
    + scale_y_discrete(
      labels = label_map,
      expand = expansion(mult = c(0.1, 0.25))
    )
    + scale_x_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
  )
  p
}

plot_bin_marginal_by_group <- function(fit, data, covs, group = "country", cred1 = 0.5, cred2 = 0.95, lims = NULL) {
  n_groups <- length(unique(data[[group]]))
  marg_diff <- create_df_marginal_bin_differences_by_group(fit, data, covs, group = group)
  marg_diff_stats <- create_df_marginal_bin_differences_stats_by_group(marg_diff, cred1 = cred1, cred2 = cred2)
  p <- (
    ggplot(
      data = marg_diff,
      mapping = aes(
        y = fct_reorder(covariate, expected_diff, .fun = ~ -median(.)),
        x = expected_diff,
      )
    )
    + ggridges::geom_density_ridges(
      aes(
        height = after_stat(scaled),
        color = group,
        fill = group
      ),
      stat = "density",
      orientation = "x",
      scale = 0.6,
      rel_min_height = 0.01,
      alpha = 0.3,
      color = NA
    )
    + geom_segment(
      data = marg_diff_stats,
      # aes(y = covariate, yend = covariate, x = lower95_exp_diff, xend = upper95_exp_diff),
      aes(y = covariate, yend = covariate, x = lower2_exp_diff, xend = upper2_exp_diff, color = group),
      linewidth = 1,
      alpha = 0.6,
      inherit.aes = FALSE
    )
    + geom_segment(
      data = marg_diff_stats,
      # aes(y = covariate, yend = covariate, x = lower95_exp_diff, xend = upper95_exp_diff),
      aes(y = covariate, yend = covariate, x = lower1_exp_diff, xend = upper1_exp_diff, color = group),
      linewidth = 1.5,
      alpha = 0.9,
      inherit.aes = FALSE
    )
    + geom_point(
      data = marg_diff_stats,
      aes(y = covariate, x = median_exp_diff, color = group),
      size = 3,
      inherit.aes = FALSE
    )
    + geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
    + facet_wrap(
      ~ group,
      ncol = n_groups
    )
    + paletteer::scale_color_paletteer_d(
      "wesanderson::Zissou1"
    )
    + paletteer::scale_fill_paletteer_d(
      "wesanderson::Zissou1"
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      strip.text = element_text(size = 15),
      strip.background = element_blank()
    )
    + labs(
      x = "Average marginal effect\n(p(violence|1) − p(violence|0))",
    )
    + scale_y_discrete(
      labels = label_map,
      expand = expansion(mult = c(0.1, 0.35))
    )
    + scale_x_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
  )
  p
}

plot_single_bin_marginal_by_group <- function(fit, data, cov, group = "country", cred1 = 0.5, cred2 = 0.95, lims = NULL) {
  marg_diff <- create_df_marginal_bin_differences_by_group(fit, data, cov, group = group)
  marg_diff_stats <- create_df_marginal_bin_differences_stats_by_group(marg_diff, cred1 = cred1, cred2 = cred2)
  p <- (
    ggplot(
      data = marg_diff,
      mapping = aes(
        y = fct_reorder(group, expected_diff, .fun = ~ -median(.)),
        x = expected_diff
      )
    )
    + ggridges::geom_density_ridges(
      aes(
        height = after_stat(scaled),
      ),
      stat = "density",
      orientation = "x",
      scale = 0.6,
      rel_min_height = 0.01,
      alpha = 0.7,
      fill = "#FDD0A2FF",
      color = NA
    )
    + geom_segment(
      data = marg_diff_stats,
      aes(y = group, yend = group, x = lower2_exp_diff, xend = upper2_exp_diff),
      linewidth = 1,
      alpha = 0.4,
      color = "#A63603FF",
      inherit.aes = FALSE
    )
    + geom_segment(
      data = marg_diff_stats,
      aes(y = group, yend = group, x = lower1_exp_diff, xend = upper1_exp_diff),
      linewidth = 1.5,
      alpha = 0.7,
      color = "#A63603FF",
      inherit.aes = FALSE
    )
    + geom_point(
      data = marg_diff_stats,
      aes(y = group, x = median_exp_diff, color = group),
      size = 3,
      color = "#A63603FF",
      inherit.aes = FALSE
    )
    + geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      strip.text = element_text(size = 15),
      strip.background = element_blank()
    )
    + labs(
      x = "Average marginal effect\n(p(violence|1) − p(violence|0))",
    )
    + scale_y_discrete(
      expand = expansion(mult = c(0.15, 0.25))
    )
    + scale_x_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
  )
  p
}

plot_mono_marginal <- function(fit, data, covs, lims = NULL) {
  p <- (
    generate_mono_marginal_df(fit, data, covs)
    %>% ggplot() 
    + geom_segment(
      aes(
        x = covariate_value,
        y = lower2,
        yend = upper2
      ),
      color = "#A63603FF",
      alpha = 0.3,
      linewidth = 1
    )
    + geom_segment(
      aes(
        x = covariate_value,
        y = lower1,
        yend = upper1
      ),
      color = "#A63603FF",
      alpha = 0.6,
      linewidth = 1.3
    )
    + geom_point(
      aes(
        x = covariate_value, 
        y = median, 
      ),
      color = "#A63603FF",
      size = 2
    )
    + facet_wrap(
      ~ factor(covariate, levels = covs),
      nrow = 1,
      scales = "free_x",
      labeller = as_labeller(label_map),
      strip.position = "bottom"
    )
    + labs(
      y = "Marginal expected\np(violence)",
    ) 
    + scale_y_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),,
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "inside",
      legend.position.inside = c(0.15, 0.8),
      strip.placement = "outside",
      strip.text = element_text(size = 15),
      strip.background = element_blank()
    )
  )
  p
}

plot_mono_marginal_by_group <- function(fit, data, covs, group = "country", dodge = 0.1, lims = NULL) {
  pd <- position_dodge(dodge)
  p <- (
    generate_mono_marginal_by_group_df(fit, data, covs, group = group)
    %>% ggplot() 
    + geom_linerange(
      aes(
        x = covariate_value,
        ymin = lower2,
        ymax = upper2,
        color = group
      ),
      position = pd,
      alpha = 0.3,
      linewidth = 1
    )
    + geom_linerange(
      aes(
        x = covariate_value,
        ymin = lower1,
        ymax = upper1,
        color = group
      ),
      position = pd,
      alpha = 0.6,
      linewidth = 1.3
    )
    + geom_point(
      aes(
        x = covariate_value, 
        y = median, 
        color = group
      ),
      position = pd,
      size = 2
    )
    + facet_wrap(
      ~ factor(covariate, levels = covs),
      nrow = 1,
      scales = "free_x",
      labeller = as_labeller(label_map),
      strip.position = "bottom"
    )
    + labs(
      y = "Marginal expected\np(violence)",
    ) 
    + paletteer::scale_color_paletteer_d(
      "wesanderson::Zissou1"
    )
    + scale_y_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),,
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      strip.placement = "outside",
      strip.text = element_text(size = 15),
      strip.background = element_blank()
    )
  )
  p
}

plot_mono_marginal_by_group_2 <- function(fit, data, covs, group = "country", lims = NULL) {
  n_groups <- length(unique(data[[group]]))
  p <- (
    generate_mono_marginal_by_group_df(fit, data, covs, group = group)
    %>% ggplot() 
    + geom_segment(
      aes(
        x = covariate_value,
        y = lower2,
        yend = upper2,
        color = group
      ),
      alpha = 0.3,
      linewidth = 1
    )
    + geom_segment(
      aes(
        x = covariate_value,
        y = lower1,
        yend = upper1,
        color = group
      ),
      alpha = 0.6,
      linewidth = 1.3
    )
    + geom_point(
      aes(
        x = covariate_value, 
        y = median, 
        color = group
      ),
      size = 2
    )
    + facet_wrap(
      ~ factor(covariate, levels = covs) + group,
      nrow = n_groups,
      scales = "free_x",
      labeller = as_labeller(label_map),
      strip.position = "bottom"
    )
    + labs(
      y = "Marginal expected\np(violence)",
    ) 
    + paletteer::scale_color_paletteer_d(
      "wesanderson::Zissou1"
    )
    + scale_y_continuous(
      limits = lims,
      labels = scales::percent_format(accuracy = 1)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),,
      axis.title.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "inside",
      legend.position.inside = c(0.15, 0.8),
      strip.placement = "outside",
      strip.text = element_text(size = 15),
      strip.background = element_blank()
    )
  )
  p
}

plot_continuous_marginal <- function(fit, data, covs, grid_points = 20, ndraws = 4000, lims = NULL) {
  p <- (
    generate_cont_marginal_df(fit, data, covs, grid_points, ndraws = ndraws)
    %>% ggplot(
      aes(
        x = covariate_value, 
        y = median, 
      )
    ) 
    + geom_ribbon(
      aes(ymin = lower2, ymax = upper2), 
      alpha = 0.3, 
      fill = col_palette[5],
      color = NA
    ) 
    + geom_ribbon(
      aes(ymin = lower1, ymax = upper1), 
      alpha = 0.6, 
      fill = col_palette[5],
      color = NA
    ) 
    + geom_line(
      linewidth = 1,
      color = col_palette[5]
    )
    + facet_wrap(
      ~ covariate,
      nrow = 1,
      scales = "free_x",
      labeller = as_labeller(label_map),
      strip.position = "bottom"
    )
    + labs(
      y = "Marginal expected\nviolence (%)",
    ) 
    + scale_y_continuous(
      limits = lims,
      labels = function(x) format(x * 100, nsmall = 0)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),,
      axis.title.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 12),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.position = "inside",
      legend.position.inside = c(0.15, 0.8),
      strip.placement = "outside",
      strip.text = element_text(size = 13),
      strip.background = element_blank()
    )
  )
  p
}


plot_single_continuous_marginal_by_group <- function(fit, data, cov, grid_points = 20, ndraws = 4000, group = "country", lims = NULL) {
  n_groups <- length(unique(data[[group]]))
  p <- (
    generate_single_cont_marginal_by_group_df(fit, data, cov, grid_points, ndraws = ndraws, group = group)
    %>% ggplot(
      mapping = aes(
        x = covariate_value, 
        y = median
      )
    ) 
    + geom_ribbon(
      aes(
        ymin = lower2, 
        ymax = upper2, 
        # fill = group
      ), 
      fill = col_palette[5],
      alpha = 0.3, 
      color = NA
    ) 
    + geom_ribbon(
      aes(
        ymin = lower1, 
        ymax = upper1, 
        # fill = group
      ), 
      fill = col_palette[5],
      alpha = 0.6, 
      color = NA
    ) 
    + geom_line(
      # mapping = aes(
      #   color = group
      # ),
      color = col_palette[5],
      linewidth = 1
    )
    + ggh4x::facet_nested_wrap(
      ~ group,
      ncol = n_groups,
      # scales = "free_y",
      remove_labels = "all",
    )
    + labs(
      x = label_map[[cov]],
      y = "Marginal expected\nviolence (%)"
    ) 
    + scale_y_continuous(
      limits = lims,
      labels = function(x) format(x * 100, nsmall = 0)
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),,
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      strip.text = element_text(size = 13),
      strip.background = element_blank()
    )
  )
  p
}

plot_calibration <- function(fit, data, outcome_var) {
  epreds <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      seed = seed
    )
    %>% group_by(across(-c(.draw, .epred)))
    %>% summarise(
      median_epred = median(.epred, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
  )
  
  rel_diag <- (
    reliabilitydiag(
      epreds$median_epred, 
      y = epreds[[outcome_var]]
    )
  )
  
  rel_diag_est <- (
    rel_diag$...1$bins
    %>% pivot_longer(
      cols = c(x_min, x_max),
      names_to = "var",
      values_to = "x"
    )
  )
  
  cases <- rel_diag$...1$cases
  breaks <- seq(0, max(cases$x) + 0.02, by = 0.02)
  hist_counts <- (
    cases
    %>% mutate(
      bins = cut(
        x = x,
        breaks = breaks,
        include.lowest = TRUE
      )
    )
    %>% count(bins, .drop = FALSE, name = "n")
    %>% ungroup()
    %>% mutate(
      prop = n / sum(n)
    )
    %>% pull(prop)
  )
  breaks_full <- rep(breaks, each = 3)
  breaks_full <- breaks_full[2:(length(breaks_full) - 1)]
  hist_counts_full <- rep(hist_counts, each = 3)
  hist_counts_full[seq(1, length(hist_counts_full), by = 3)] <- 0
  hist_counts_full <- c(hist_counts_full, 0)
  max_x <- max(
    max(breaks_full),
    max(rel_diag_est$x)
  )
  cons_reg <- (
    rel_diag$...1$regions
    %>% filter(x <= max_x)
  )
  # upper_lvls <- rel_diag$...1$regions$upper
  # cep <- rel_diag_est$CEP_pav
  # max_y <- max(
  #   max(hist_counts_full),
  #   max(upper_lvls[rel_diag_est$x < max(rel_diag_est$x)]),
  #   max(cep[rel_diag_est$x < max(rel_diag_est$x)])
  # )
  
  cal_plot <- (
    ggplot()
    # + geom_point(
    #   data = sample_frac(rel_diag$...1$cases, 0.2),
    #   mapping = aes(
    #     x = x,
    #     y = y
    #   ),
    #   alpha = 0.1,
    #   color = col_palette[5],
    #   fill = col_palette[6],
    #   position = position_jitter(height = 0.05)
    # )
    # + geom_histogram(
    #   data = rel_diag$...1$cases,
    #   mapping = aes(
    #     x = x,
    #     y = after_stat(width * density)
    #   ),
    #   fill = NA,
    #   color = col_palette[5],
    #   linewidth = 1.2,
    #   binwidth = 0.02
    # )
    + geom_line(
      mapping = aes(
        x = breaks_full,
        y = hist_counts_full
      ),
      color = col_palette[5],
      linewidth = 1.1
    )
    + geom_ribbon(
      data = cons_reg,
      mapping = aes(
        x = x,
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.5,
      fill = col_palette[3]
    )
    + geom_line(
      data = rel_diag_est,
      mapping = aes(
        x = x,
        y = CEP_pav
      ),
      color = col_palette[1],
      linewidth = 1.2
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
    )
    + geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linetype = "dashed"
    )
    + scale_x_continuous(
      limits = c(0, max(rel_diag_est$x))
    )
    # + scale_y_continuous(
    #   limits = c(0, max_y)
    # )
    + coord_cartesian(clip = "off")
    + labs(
      x = "Median Posterior\nPredictive",
      y = "Observed"
    )
  )
  
  cal_plot
}

plot_calibration_weighted_new <- function(fit, data, outcome_var) {
  epreds <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      seed = seed
    )
    %>% mutate(
      weighted_pred = (.epred * sample_weight) / sum(sample_weight, na.rm = TRUE)
    )
    %>% group_by(across(-c(.draw, .epred, weighted_pred)))
    %>% summarise(
      median_epred = quantile(weighted_pred, p = 0.5),
      .groups = "drop"
    )
    %>% ungroup()
  )
  
  rel_diag <- (
    reliabilitydiag(
      epreds$median_epred, 
      y = epreds[[outcome_var]]
    )
  )
  
  rel_diag_est <- (
    rel_diag$...1$bins
    %>% pivot_longer(
      cols = c(x_min, x_max),
      names_to = "var",
      values_to = "x"
    )
  )
  
  cases <- rel_diag$...1$cases
  breaks <- seq(0, max(cases$x) + 0.02, by = 0.02)
  hist_counts <- (
    cases
    %>% mutate(
      bins = cut(
        x = x,
        breaks = breaks
      )
    )
    %>% group_by(bins)
    %>% summarise(
      prop = n() / nrow(.)
    )
    %>% pull(prop)
  )
  breaks_full <- rep(breaks, each = 3)
  breaks_full <- breaks_full[2:(length(breaks_full) - 1)]
  hist_counts_full <- rep(hist_counts, each = 3)
  hist_counts_full[seq(1, length(hist_counts_full), by = 3)] <- 0
  hist_counts_full <- c(hist_counts_full, 0)
  max_x <- max(
    max(breaks_full),
    max(rel_diag_est$x)
  )
  cons_reg <- (
    rel_diag$...1$regions
    %>% filter(cons_reg$x <= max_x)
  )
  # upper_lvls <- rel_diag$...1$regions$upper
  # cep <- rel_diag_est$CEP_pav
  # max_y <- max(
  #   max(hist_counts_full),
  #   max(upper_lvls[rel_diag_est$x < max(rel_diag_est$x)]),
  #   max(cep[rel_diag_est$x < max(rel_diag_est$x)])
  # )
  
  cal_plot <- (
    ggplot()
    # + geom_point(
    #   data = sample_frac(rel_diag$...1$cases, 0.2),
    #   mapping = aes(
    #     x = x,
    #     y = y
    #   ),
    #   alpha = 0.1,
    #   color = col_palette[5],
    #   fill = col_palette[6],
    #   position = position_jitter(height = 0.05)
    # )
    # + geom_histogram(
    #   data = rel_diag$...1$cases,
    #   mapping = aes(
    #     x = x,
    #     y = after_stat(width * density)
    #   ),
    #   fill = NA,
    #   color = col_palette[5],
    #   linewidth = 1.2,
    #   binwidth = 0.02
    # )
    + geom_line(
      mapping = aes(
        x = breaks_full,
        y = hist_counts_full
      ),
      color = col_palette[5],
      linewidth = 1.1
    )
    + geom_ribbon(
      data = cons_reg,
      mapping = aes(
        x = x,
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.5,
      fill = col_palette[3]
    )
    + geom_line(
      data = rel_diag_est,
      mapping = aes(
        x = x,
        y = CEP_pav
      ),
      color = col_palette[1],
      linewidth = 1.2
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
    )
    + geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linetype = "dashed"
    )
    + scale_x_continuous(
      limits = c(0, max(rel_diag_est$x))
    )
    # + scale_y_continuous(
    #   limits = c(0, max_y)
    # )
    + coord_cartesian(clip = "off")
    + labs(
      x = "Median Posterior\nPredictive",
      y = "Observed"
    )
  )
  
  cal_plot
}

plot_calibration_weighted <- function(fit, data, outcome_var) {
  
  epreds <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = 100,
      seed = seed
    )
    %>% group_by(across(-c(.draw, .epred)))
    %>% summarise(
      median_epred = median(.epred, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
    %>% arrange(median_epred)
  )
  
  pava <- isotone::gpava(
    z = epreds$median_epred,
    y = epreds[[outcome_var]],
    weights = epreds$sample_weight,
    solver = weighted.mean
  )
  
  cal_plot <- (
    ggplot()
    + geom_point(
      data = epreds,
      mapping = aes(
        x = median_epred,
        y = !!sym(outcome_var)
      ),
      alpha = 0.05,
      color = col_palette[5],
      fill = col_palette[6],
      position = position_jitter(height = 0.05)
    )
    + geom_line(
      mapping = aes(
        x = pava$z,
        y = pava$x
      ),
      color = col_palette[1],
      linewidth = 1.2
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
    )
    + geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linetype = "dashed"
    )
    + scale_x_continuous(
      limits = c(min(pava$z) - 0.002, max(pava$z) + 0.005)
    )
    + labs(
      x = "Median Posterior Predictive",
      y = "Observed"
    )
  )
  cal_plot
}


plot_moran_residuals <- function(fit, data, shapes) {

  df_resid <- (
    data
    %>% add_residual_draws(fit)
    %>% group_by(country, adm_name)
    %>% summarise(
      residual = mean(.residual, na.rm = TRUE),
      .groups = "drop"
    )
  )

  moran_plot <- (
    spdep::moran.plot(
      df_resid$residual,
      spdep::nb2listw(spdep::poly2nb(shapes, row.names = shapes$adm_name)),
      # labels = as.character(df_resid$country),
      labels = FALSE,
      pch = 19,
      tck = -0.05,
      xlab = "Probability of Violence (residuals)",
      ylab = "Spatial Lagged Values",
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      axis.title = element_text(size = 16, family = "sans"),
      axis.text.y = element_text(size = 15, family = "sans"),
      axis.text.x = element_text(size = 15, family = "sans"),
      axis.ticks.length = unit(1, "pt"),
      legend.background = element_blank()
    )
  )
  moran_plot
}

plot_moran <- function(data, shapes, fit = NULL, outcome_var = NULL, residuals = FALSE, ndraws = 4000) {
  
  if (residuals) {
    xlab <- "p(viol, res.)"
    df <- (
      data
      %>% add_residual_draws(fit, ndraws = ndraws)
      %>% group_by(country, adm_name)
      %>% summarise(
        value = mean(.residual, na.rm = TRUE),
        .groups = "drop"
      )
    )
  } else {
    xlab <- "p(viol, obs.)"
    df <- (
      data
      %>% group_by(country, adm_name)
      %>% summarise(
        value = mean(!!sym(outcome_var), na.rm = TRUE),
        .groups = "drop"
      )
    )
  }
  
  moran_data <- (
    spdep::moran.plot(
      df$value,
      spdep::nb2listw(spdep::poly2nb(shapes, row.names = shapes$adm_name)),
      labels = FALSE,
      plot = FALSE,
      return_df = TRUE
    )
  )
  
  moran_plot <- (
    ggplot(
      data = moran_data,
      mapping = aes(x = x, y = wx)
    )
    + geom_hline(
      yintercept = mean(moran_data$wx),
      linetype = "dotted"
    )
    + geom_vline(
      xintercept = mean(moran_data$x),
      linetype = "dotted"
    )
    + geom_point(
      shape = 19,
      color = col_palette[5]
    )
    + geom_smooth(
      method = "lm",
      se = TRUE,
      color = col_palette[7],
      fill = col_palette[8],
      alpha = 0.5,
      linetype = "dashed"
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.background = element_blank(),
      axis.title = element_text(size = 16, family = "sans"),
      axis.text = element_text(size = 15, family = "sans"),
    )
    + coord_cartesian(clip = "off")
    + labs(
      x = xlab,
      y = "Lagged Values"
    )
  )
  moran_plot
}

plot_moran_raw <- function(data, shapes, outcome_var) {

  df_obs_prob <- (
    data
    %>% group_by(country, adm_name)
    %>% summarise(
      viol_prob = mean(!!sym(outcome_var), na.rm = TRUE),
      .groups = "drop"
    )
  )

  moran_data <- (
    spdep::moran.plot(
      df_obs_prob$viol_prob,
      spdep::nb2listw(spdep::poly2nb(shapes, row.names = shapes$adm_name)),
      labels = FALSE,
      plot = FALSE,
      return_df = TRUE
    )
  )

  moran_plot <- (
    ggplot(
      data = moran_data,
      mapping = aes(x = x, y = wx)
    )
    + geom_hline(
      yintercept = mean(moran_data$wx),
      linetype = "dotted"
    )
    + geom_vline(
      xintercept = mean(moran_data$x),
      linetype = "dotted"
    )
    + geom_point(
      shape = 19,
      color = col_palette[5]
    )
    + geom_smooth(
      method = "lm",
      se = TRUE,
      color = col_palette[7],
      fill = col_palette[8],
      alpha = 0.5,
      linetype = "dashed"
    )
    + theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      plot.background = element_blank(),
      axis.title = element_text(size = 16, family = "sans"),
      axis.text = element_text(size = 15, family = "sans"),
    )
    + labs(
      x = "Probability of Violence (observed)",
      y = "Spatial Lagged Values"
    )
  )
  moran_plot
}