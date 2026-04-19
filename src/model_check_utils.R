extract_post <- function(model, pars) {
  dfs = list()
  post <- rstan::extract(model, pars = pars)
  for (i in seq(1:length(pars))) {
    ncols <- ncol(post[[pars[i]]])
    if (is.na(ncols)) {
      dfs[[i]] <- (
        as.data.frame(post[[pars[i]]])
        %>% set_names(pars[i])
      )
    } else {
      dfs[[i]] <- (
        as.data.frame(post[[pars[i]]])
        %>% set_names(paste0(pars[i], "[", 1:ncols, "]"))
      )
    }
  }
  bind_cols(dfs)
}

extract_chains <- function(model, pars) {
  dfs <- list()
  chains_data <- rstan::extract(model, pars = pars, permuted = FALSE)
  for (i in 1:dim(chains_data)[2]) {
    chain <- chains_data[, i, ]
    dfs[[i]] <- (
      as_tibble(
        matrix(
          chain,
          nrow = dim(chain)[1],
          ncol = dim(chain)[2],
          dimnames = list(
            NULL,
            names(chain[1, ])
          )
        )
      )
    )
    param_levels <- names(dfs[[i]])
    dfs[[i]] <- (
      dfs[[i]]
      %>% mutate(
        chain = as.factor(i),
        step = row_number()
      )
      %>% pivot_longer(cols = -c(step, chain), names_to = "parameter", values_to = "value")
      %>% mutate(
        parameter = factor(parameter, levels = param_levels),
      )
    )
  }
  bind_rows(dfs)
}

extract_ppc <- function(model, n_samples = 500, df_observed) {
  ppc <- rstan::extract(model)$sv_pred[1:n_samples,]
  rownames(ppc) <- 1:nrow(ppc)
  df_ppc <- (
    as_tibble(t(ppc))
    %>% bind_cols(
      country = df_observed$country,
      adm_name = df_observed$adm_name,
      rural = df_observed$rural,
      over_18 = df_observed$over_18,
      school_level = df_observed$school_level,
      live_with_male_relative = df_observed$live_with_male_relative,
      clim = df_observed$clim,
      sv_true = df_observed$sv
    )
    %>% pivot_longer(
      -(country:sv_true),
      names_to = "sample",
      values_to = "sv_pred"
    )
  )
}


create_ranef_df <- function(ranef_arr, group_name) {
  df <- as.data.frame.table(ranef_arr, responseName = "value", stringsAsFactors = FALSE)
  colnames(df)[1:3] <- c("iter", "level", "term")
  df$group <- group_name
  df$value <- as.numeric(df$value)
  df
}


calculate_marginal_bin_difference <- function(fit, data, covariate) {
  newdata0 <- data
  newdata1 <- data
  newdata0[[covariate]] <- 0
  newdata1[[covariate]] <- 1
  
  epreds0 <- posterior_epred(fit, newdata = newdata0, re_formula = NA)
  epreds1 <- posterior_epred(fit, newdata = newdata1, re_formula = NA)
  
  diff_means <- rowMeans(epreds1 - epreds0)
  
  tibble(
    covariate = covariate,
    expected_diff = diff_means
  )
}

calculate_marginal_bin_difference_by_group <- function(fit, data, covariate, group = "country") {
  
  newdata0 <- data
  newdata1 <- data
  newdata0[[covariate]] <- 0
  newdata1[[covariate]] <- 1
  
  group_vec <- data[[group]]
  groups <- unique(group_vec)
  group_indices <- split(seq_len(nrow(newdata0)), group_vec)
  
  epreds0 <- posterior_epred(fit, newdata = newdata0, re_formula = NULL)
  epreds1 <- posterior_epred(fit, newdata = newdata1, re_formula = NULL)
  epred_diff <- epreds1 - epreds0
  
  diff_means <- lapply(
    group_indices,
    function(country_idxs) {
      if (length(country_idxs) == 1) {
        epred_diff[, country_idxs]
      } else {
        rowMeans(epred_diff[, country_idxs])
      }
    }
  )
  
  group_names <- names(diff_means)
  
  df_list <- lapply(
    seq_along(diff_means),
    function(i) {
      tibble(
        group = group_names[i],
        covariate = covariate,
        expected_diff = diff_means[[group_names[i]]]
      )
    }
  )
  
  bind_rows(df_list)

}

create_df_marginal_bin_differences <- function(fit, data, covs) {
  (
    lapply(covs, function(cov) {
      calculate_marginal_bin_difference(fit, data, cov)
    })
    %>% bind_rows()
  )
}

create_df_marginal_bin_differences_by_group <- function(fit, data, covs, group = "country") {
  (
    lapply(covs, function(cov) {
      calculate_marginal_bin_difference_by_group(fit, data, cov, group = group)
    })
    %>% bind_rows()
  )
}

create_df_marginal_bin_differences_stats <- function(marg_diff, cred1 = 0.5, cred2 = 0.95) {
  (
    marg_diff
    %>% group_by(covariate)
    %>% summarise(
      median_exp_diff = median(expected_diff),
      mean_exp_diff = mean(expected_diff),
      lower1_exp_diff = bayestestR::hdi(expected_diff, ci = cred1)$CI_low,
      upper1_exp_diff = bayestestR::hdi(expected_diff, ci = cred1)$CI_high,
      lower2_exp_diff = bayestestR::hdi(expected_diff, ci = cred2)$CI_low,
      upper2_exp_diff = bayestestR::hdi(expected_diff, ci = cred2)$CI_high,
      # lower95_exp_diff = quantile(expected_diff, 0.025),
      # upper95_exp_diff = quantile(expected_diff, 0.975),
      .groups = "drop"
    )
    %>% arrange(median_exp_diff)
  )
}

create_df_marginal_bin_differences_stats_by_group <- function(marg_diff, cred1 = 0.5, cred2 = 0.95) {
  (
    marg_diff
    %>% group_by(group, covariate)
    %>% summarise(
      median_exp_diff = median(expected_diff),
      mean_exp_diff = mean(expected_diff),
      lower1_exp_diff = bayestestR::hdi(expected_diff, ci = cred1)$CI_low,
      upper1_exp_diff = bayestestR::hdi(expected_diff, ci = cred1)$CI_high,
      lower2_exp_diff = bayestestR::hdi(expected_diff, ci = cred2)$CI_low,
      upper2_exp_diff = bayestestR::hdi(expected_diff, ci = cred2)$CI_high,
      # lower95_exp_diff = quantile(expected_diff, 0.025),
      # upper95_exp_diff = quantile(expected_diff, 0.975),
      .groups = "drop"
    )
    %>% ungroup()
    %>% arrange(median_exp_diff)
  )
}

calculate_marginal <- function(fit, data, cov, grid_vals, ndraws = 1000) {

  draws_list <- vector("list", length(grid_vals))
  for (i in seq_along(grid_vals)) {
    newdata <- data
    newdata[[cov]] <- grid_vals[i]
    ep <- posterior_epred(fit, newdata = newdata, re_formula = NA, ndraws=ndraws)
    draws_list[[i]] <- rowMeans(ep)
  }
  
  scale <- attr(data[[cov]], "scaled:scale")
  if (!is.null(scale)) {
    center <- attr(data[[cov]], "scaled:center")
    grid_vals <- (grid_vals * scale) + center
  }
  
  list(
    covariate_value = grid_vals,
    expected_viol = draws_list
  )
}

calculate_marginal_by_group <- function(
    fit, 
    data, 
    cov, 
    grid_vals,
    ndraws = 1000,
    group = "country"
  ) {
  
  group_vec <- data[[group]]
  groups <- unique(group_vec)
  
  draws_list <- vector("list", length(grid_vals))
  
  for (i in seq_along(grid_vals)) {
    newdata <- data
    newdata[[cov]] <- grid_vals[i]
    group_indices <- split(seq_len(nrow(newdata)), group_vec)
    ep <- posterior_epred(fit, newdata = newdata, re_formula = NULL, ndraws = ndraws)
    draws_list[[i]] <- sapply(
      group_indices,
      function(country_idxs) {
        if (length(country_idxs) == 1) {
          ep[, country_idxs]
        } else {
          rowMeans(ep[, country_idxs])
        }
      }
    )
  }
  
  scale <- attr(data[[cov]], "scaled:scale")
  if (!is.null(scale)) {
    center <- attr(data[[cov]], "scaled:center")
    grid_vals <- (grid_vals * scale) + center
  }
  
  list(
    covariate_value = grid_vals,
    expected_viol = draws_list
  )
}

summarize_draws <- function(draws_list, cred1 = 0.5, cred2 = 0.95) {
  grid <- draws_list$covariate_value
  draws <- draws_list$expected_viol
  df_list <- lapply(seq_along(grid), function(i) {
    draws_i <- draws[[i]]
    hdi1_i <- bayestestR::hdi(draws_i, ci = cred1)
    hdi2_i <- bayestestR::hdi(draws_i, ci = cred2)
    tibble(
      covariate_value = grid[i],
      median = median(draws_i),
      mean = mean(draws_i),
      lower1 = hdi1_i$CI_low,
      upper1 = hdi1_i$CI_high,
      lower2 = hdi2_i$CI_low,
      upper2 = hdi2_i$CI_high
    )
  })
  bind_rows(df_list)
}

summarize_draws_by_group <- function(draws_list, cred1 = 0.5, cred2 = 0.95) {
  groups <- dimnames(draws_list$expected_viol[[1]])[[2]]
  grid <- draws_list$covariate_value
  draws <- draws_list$expected_viol
  df_list_grid <- lapply(seq_along(grid), function(i) {
    draws_i <- draws[[i]]
    hdi1_i <- lapply(seq_len(ncol(draws_i)), function(j) bayestestR::hdi(draws_i[, j], ci = cred1))
    names(hdi1_i) <- colnames(draws_i)
    hdi2_i <- lapply(seq_len(ncol(draws_i)), function(j) bayestestR::hdi(draws_i[, j], ci = cred2))
    names(hdi2_i) <- colnames(draws_i)
    df_list_group <- lapply(seq_along(groups), function(group_idx) {
      group_name <- groups[group_idx]
      tibble(
        group = groups[group_idx],
        covariate_value = grid[i],
        median = median(draws_i[, group_idx]),
        mean = mean(draws_i[, group_idx]),
        lower1 = hdi1_i[[group_name]]$CI_low,
        upper1 = hdi1_i[[group_name]]$CI_high,
        lower2 = hdi2_i[[group_name]]$CI_low,
        upper2 = hdi2_i[[group_name]]$CI_high
      )
    })
    bind_rows(df_list_group)
  })
  bind_rows(df_list_grid)
}

generate_mono_marginal_df <- function(fit, data, covs, ndraws = 1000) {
  results <- list()
  for (cov in covs) {
    grid_vals <- sort(as.integer(as.character(unique(data[[cov]]))))
    mg <- calculate_marginal(fit, data, cov, grid_vals, ndraws = ndraws)
    s <- (
      summarize_draws(mg)
      %>% mutate(
        covariate = cov
      )
    )
    results[[cov]] <- s
  }
  bind_rows(results) 
}

generate_mono_marginal_by_group_df <- function(fit, data, covs, group = "country") {
  results <- list()
  for (cov in covs) {
    grid_vals <- as.integer(levels(data[[cov]]))
    mg <- calculate_marginal_by_group(fit, data, cov, grid_vals, group = group)
    s <- (
      summarize_draws_by_group(mg)
      %>% mutate(
        covariate = cov
      )
    )
    results[[cov]] <- s
  }
  bind_rows(results) 
}

generate_cont_marginal_df <- function(fit, data, covs, grid_points = 20, ndraws = 1000) {
  results <- list()
  for (cov in covs) {
    rng <- range(data[[cov]], na.rm = TRUE)
    grid_vals <- seq(rng[1], rng[2], length.out = grid_points)
    mg <- calculate_marginal(fit, data, cov, grid_vals, ndraws = ndraws)
    s <- (
      summarize_draws(mg)
      %>% mutate(
        covariate = cov
      )
    )
    results[[cov]] <- s
  }
  bind_rows(results) 
}

generate_single_cont_marginal_by_group_df <- function(fit, data, cov, grid_points = 20, ndraws = 1000, group = "country") {
  rng <- range(data[[cov]], na.rm = TRUE)
  grid_vals <- seq(rng[1], rng[2], length.out = grid_points)
  mg <- calculate_marginal_by_group(fit, data, cov, grid_vals, ndraws = ndraws, group = group)
  s <- (
    summarize_draws_by_group(mg)
    %>% mutate(
      covariate = cov
    )
  )
  s
}

get_calibration_scores <- function(fit, data, outcome_var, model_name) {
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
    reliabilitydiag::reliabilitydiag(
      epreds$median_epred, 
      y = epreds[[outcome_var]]
    )
  )
  
  summary <- (
    utils::getS3method("summary", "reliabilitydiag")(rel_diag, score = "brier")
    %>% mutate(
      forecast = model_name
    )
  )
  summary
}
