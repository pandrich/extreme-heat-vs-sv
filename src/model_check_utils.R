calculate_marginal <- function(fit, data, cov, grid_vals, ndraws = 1000) {

  draws_list <- vector("list", length(grid_vals))
  for (i in seq_along(grid_vals)) {
    newdata <- data
    newdata[[cov]] <- grid_vals[i]
    ep <- brms::posterior_epred(fit, newdata = newdata, re_formula = NA, ndraws=ndraws)
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
    ep <- brms::posterior_epred(fit, newdata = newdata, re_formula = NULL, ndraws = ndraws)
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
