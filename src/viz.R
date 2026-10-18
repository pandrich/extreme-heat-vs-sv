plot_ppc_rootogram <- function(fit, data, outcome_var, cred1 = 0.5, cred2 = 0.95) {
  df_obs <- (
    data
    %>% group_by(country, adm_name)
    %>% summarise(
      mean_sv_obs := sum(!!sym(outcome_var) * sample_weight, na.rm = TRUE) / sum(sample_weight, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
  )
  ndraws = min(ndraws(fit), 100)
  df_ppc <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = ndraws,
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
      ), 
      fill = col_palette[5],
      alpha = 0.3, 
      color = NA
    ) 
    + geom_ribbon(
      aes(
        ymin = lower1, 
        ymax = upper1, 
      ), 
      fill = col_palette[5],
      alpha = 0.6, 
      color = NA
    ) 
    + geom_line(
      color = col_palette[5],
      linewidth = 1
    )
    + ggh4x::facet_nested_wrap(
      ~ group,
      ncol = n_groups,
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
  ndraws = min(ndraws(fit), 100)
  epreds <- (
    data
    %>% add_epred_draws(
      fit,
      ndraws = ndraws,
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
  
  cal_plot <- (
    ggplot()
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
    + coord_cartesian(clip = "off")
    + labs(
      x = "Median Posterior\nPredictive",
      y = "Observed"
    )
  )
  
  cal_plot
}