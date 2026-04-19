build_temp_daily <- function(rast_temp, adm_shapes) {
  df_temp <- (
    rast_temp
    %>% exactextractr::exact_extract(
      adm_shapes,
      fun = "max",
      # fun = "weighted_mean",
      # weights = "area",
      max_cells_in_memory = 176067000,
      append_cols = c("country", "adm_code", "adm_name")
    )
  )
  colnames(df_temp) <- c("country", "adm_code", "adm_name", as.character(terra::time(rast_temp)))
  
  df_temp <- (
    df_temp
    %>% pivot_longer(
      cols = -c("country", "adm_code", "adm_name"),
      names_to = "date",
      values_to = "temp_2m_max"
    )
    %>% mutate(
      date = as.Date(date, format = "%Y-%m-%d"),
      temp_2m_max = temp_2m_max - 273.15
    )
  )
  df_temp
}


build_utci_daily <- function(utci_data_dir, country_map, adm_shapes) {
  print(paste0("Starting UTCI data processing"))
  country_dirs <- list.dirs(utci_data_dir, recursive = FALSE)
  country_isos <- sapply(str_split(country_dirs, "/"), tail, 1)
  country_dirs <- setNames(country_dirs, country_isos)
  dfs_utci <- list()
  for (c in country_isos) {
    print(paste0("Processing UTCI data for: ", country_map[[c]]))
    shape <- adm_shapes %>% filter(country == country_map[[c]])
    rast_files <- list.files(country_dirs[[c]], pattern = "^ECMWF_utci_daily_stats", full.names = TRUE, recursive = TRUE)
    thresholds <- round(quantile(seq_along(rast_files), probs = c(0, 0.25, 0.5, 0.75, 1)))
    dfs_utci_country <- list()
    for (i in seq_along(rast_files)) {
      if (i %in% thresholds) {
        thresh_idx <- which(thresholds == i)
        print(paste0(names(thresh_idx), " complete"))
      }
      rast_utci <- terra::rast(rast_files[[i]])["utci_daily_max"]
      date <- str_split_i(terra::time(rast_utci)[[1]], " ", 1)
      dfs_utci_country[[i]] <- (
        rast_utci
        %>% exactextractr::exact_extract(
          shape,
          # fun = "max",
          fun = "weighted_mean",
          weights = "area",
          max_cells_in_memory = 176067000,
          append_cols = c("country", "adm_code", "adm_name"),
          progress = FALSE
        )
        %>% mutate(
          date = as.Date(date, format = "%Y-%m-%d"),
          utci_max = max - 273.15
        )
      )
    }
    dfs_utci[[c]] <- (
      bind_rows(dfs_utci_country)
      %>% select(-max)
      %>% mutate(
        country = country_map[[c]]
      )
      %>% relocate(country)
    )
  }
  df_utci <- bind_rows(dfs_utci)
}

calculate_temp_12m_stats <- function(df_temp_daily, roll = FALSE, roll_days = 3) {
  variable_name = "temp"
  if (roll) {
    df_temp_daily <- (
      df_temp_daily
      %>% arrange(country, adm_code, adm_name, date)
      %>% group_by(country, adm_code, adm_name)
      %>% mutate(
        temp_2m_max = slider::slide_index_dbl(
          .x = temp_2m_max,
          .i = date,
          .f = ~ mean(.x, na.rm = TRUE),
          .before = ~ .x %m-% days(roll_days),
          .complete = TRUE
        )
      )
      %>% ungroup()
    )
    variable_name = paste0(variable_name, "_roll")
  }
  df_temp_12m <- (
    df_temp_daily
    %>% group_by(country, adm_code, adm_name)
    %>% mutate(
      hist_mean = mean(temp_2m_max, na.rm = TRUE),
      hist_sd = sd(temp_2m_max, na.rm = TRUE),
    )
    %>% ungroup()
    %>% mutate(
      sd_temp = (temp_2m_max - hist_mean) / hist_sd,
      above_1sd = as.numeric(sd_temp > 1)
    )
    %>% select(-c(hist_mean, hist_sd, sd_temp))
    %>% arrange(country, adm_code, adm_name, date)
    %>% group_by(country, adm_code, adm_name)
    %>% mutate(
      across(
        -c(date, temp_2m_max),
        ~ (
          slider::slide_index_dbl(
            .x = (.),
            .i = date,
            .f = ~ sum(.x, na.rm = TRUE),
            .before = ~ .x %m-% months(12),
            .complete = TRUE
          )
        ),
        .names = "{.col}_12m"
      ),
      max_12m = slider::slide_index_dbl(
        .x = temp_2m_max,
        .i = date,
        .f = ~ max(.x, na.rm = TRUE),
        .before = ~ .x %m-% months(12),
        .complete = TRUE
      )
    )
    %>% ungroup()
    %>% select(country, adm_code, adm_name, date, contains("12m"))
    %>% filter(!is.na(above_1sd_12m))
    %>% pivot_longer(
      cols = -c(country, adm_code, adm_name, date),
      names_to = "variable",
      values_to = "value"
    ) 
    %>% group_by(country, adm_code, adm_name, variable)
    %>% mutate(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    )
    %>% ungroup()
    %>% mutate(
      dev_value = if_else(
        sd_value == 0,
        0,
        (value - mean_value) / sd_value
      ),
      variable = paste0(variable_name, "_", variable)
    )
    %>% filter(!is.na(value))
    %>% pivot_wider(
      id_cols = c(country, adm_code, adm_name, date),
      names_from = "variable",
      values_from = c(value, dev_value)
    )
  )
  df_temp_12m
}


calculate_utci_12m_stats <- function(df_utci_daily, roll = FALSE, roll_days = 3) {
  variable_name = "utci"
  if (roll) {
    df_utci_daily <- (
      df_utci_daily
      %>% arrange(country, adm_code, adm_name, date)
      %>% group_by(country, adm_code, adm_name)
      %>% mutate(
        utci_max = slider::slide_index_dbl(
          .x = utci_max,
          .i = date,
          .f = ~ mean(.x, na.rm = TRUE),
          .before = ~ .x %m-% days(roll_days),
          .complete = TRUE
        )
      )
      %>% ungroup()
    )
    variable_name = paste0(variable_name, "_roll")
  } 
  df_utci_12m <- (
    df_utci_daily
    %>% mutate(
      above_32 = as.numeric(utci_max > 32)
    )
    %>% arrange(country, adm_code, adm_name, date)
    %>% group_by(country, adm_code, adm_name)
    %>% mutate(
      across(
        -c(date, utci_max),
        ~ (
          slider::slide_index_dbl(
            .x = (.),
            .i = date,
            .f = ~ sum(.x, na.rm = TRUE),
            .before = ~ .x %m-% months(12),
            .complete = TRUE
          )
        ),
        .names = "{.col}_12m"
      ),
      max_12m = slider::slide_index_dbl(
        .x = utci_max,
        .i = date,
        .f = ~ max(.x, na.rm = TRUE),
        .before = ~ .x %m-% months(12),
        .complete = TRUE
      )
    )
    %>% ungroup()
    %>% select(country, adm_code, adm_name, date, contains("12m"))
    %>% filter(!is.na(above_32_12m))
    %>% pivot_longer(
      cols = -c(country, adm_code, adm_name, date),
      names_to = "variable",
      values_to = "value"
    ) 
    %>% group_by(country, adm_code, adm_name, variable)
    %>% mutate(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE)
    )
    %>% ungroup()
    %>% mutate(
      dev_value = if_else(
        sd_value == 0,
        0,
        (value - mean_value) / sd_value
      ),
      variable = paste0(variable_name, "_", variable)
    )
    %>% filter(!is.na(value))
    %>% pivot_wider(
      id_cols = c(country, adm_code, adm_name, date),
      names_from = "variable",
      values_from = c(value, dev_value)
    )
  )
  df_utci_12m
}
