make_plant_data <- function(df, global_df) {
  years <- unique(df$year)
  list(N= nrow(df),
       y= df$abundance,
       Nyr= length(years),
       year_id= as.integer(factor(df$year)),
       DOYs= df$DOYs,
       year= df$year,
       year_vec = (years-mean(years))/sd(years),
       Nplots  = length(unique(df$plot_id)),
       plot_id = as.integer(factor(df$PlotID)),
       DOY_sd  = sd(global_df$DOY),
       DOY_mean = mean(global_df$DOY))
}

fit_plant_model <- function(
    species_name,
    data,
    model,
    out_dir = "models",
    csv_dir = "cmdstan_csv") {

  message("Fit model for: ", species_name)

  sp_df <- data %>% filter(taxon == species_name)

  plant_data <- make_plant_data(sp_df, data)

  # Ensure output directories exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  species_csv_dir <- file.path(csv_dir, species_name)
  if (!dir.exists(species_csv_dir)) dir.create(species_csv_dir, recursive = TRUE)

  fit <- model$sample(
    data = plant_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 2000,
    iter_warmup = 500,
    output_dir = species_csv_dir,
    save_warmup = TRUE)

  saveRDS(fit, file = file.path(out_dir, paste0("phenology_", species_name, ".rds")))

  return(fit)
}

plot_plant_preds <- function(species_name, data,
                             model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models\\",
                             out_dir = "predictions") {

  message("Plotting predictions of: ", species_name)

  # Open model
  model_file <- file.path(model_dir, paste0("phenology_", species_name, ".rds"))
  if (!file.exists(model_file)) stop("Model file not found: ", model_file)
  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Subset data
  sp_df <- dplyr::filter(data, taxon == species_name)

  mu         <- fit$draws("mu", format="draws_matrix")
  width      <- fit$draws("width", format="draws_matrix")
  alpha_year <- fit$draws("alpha_year", format="draws_matrix")
  u_plot     <- fit$draws("u_plot", format="draws_matrix")
  u_plot_mu  <- fit$draws("u_plot_mu", format="draws_matrix")

  # Setup
  #log_obs_days <- log(mean(sp_df$TrapDays, na.rm = TRUE))

  DOY <- seq(-2, 2, length.out=100)

  # PDF
  pdf(file.path(out_dir, paste0(species_name, "_preds.pdf")))
  par(mfrow=c(3, 2))

  plot_ids <- sort(unique(sp_df$plot_id))
  year_ids <- sort(unique(sp_df$year))  # actual year numbers
  P <- length(plot_ids)
  cols <- rainbow(P)

  for (y_idx in seq_along(year_ids)) {
    y <- year_ids[y_idx]  # actual year

    for (pl in seq_along(plot_ids)) {
      p_id <- plot_ids[pl]
      pl_idx <- which(plot_ids == p_id)  # correct column index

      eta_hat <- numeric(length(DOY))

      for (nd in seq_along(DOY)) {
        eta_post <- exp(alpha_year[, y_idx] + u_plot[, pl_idx] -
                          ((DOY[nd] - (mu[, y_idx] + u_plot_mu[, pl_idx]))^2 / (width[, y_idx]^2)))

        eta_hat[nd] <- median(eta_post)
      }

      if (pl == 1) {
        plot(DOY, eta_hat, type="l", col=cols[pl],
             main=paste(species_name, "Year", y),
             xlab="Scaled DOY", ylab="Abundance")
      } else {
        lines(DOY, eta_hat, col=cols[pl])
      }

      # raw points
      ind <- sp_df$plot_id == p_id & sp_df$year == y
      points(sp_df$DOYs[ind], sp_df$abundance[ind],
             col=cols[pl], pch=16)
    }
  }

  dev.off()
}

diagnose_peak_parameters <- function(
    species_name,
    data,
    model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models\\",
    out_dir = "diagnostics",
    doy_min = -2,
    doy_max = 2
) {
  message("Diagnosing peak parameters for: ", species_name)

  # -----------------------------
  # Load model
  # -----------------------------
  model_file <- file.path(model_dir, paste0("phenology_", species_name, ".rds"))

  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }

  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # -----------------------------
  # Subset species data
  # -----------------------------
  sp_df <- dplyr::filter(data, taxon == species_name)

  if (nrow(sp_df) == 0) {
    stop("No data found for species: ", species_name)
  }

  # -----------------------------
  # Extract posterior draws
  # -----------------------------
  mu        <- fit$draws("mu", format = "draws_matrix")
  u_plot_mu <- fit$draws("u_plot_mu", format = "draws_matrix")
  width     <- fit$draws("width", format = "draws_matrix")

  # Optional, useful for checking abundance scale
  alpha_year <- fit$draws("alpha_year", format = "draws_matrix")
  u_plot     <- fit$draws("u_plot", format = "draws_matrix")

  # -----------------------------
  # Setup IDs
  # -----------------------------
  plot_ids <- sort(unique(sp_df$plot_id))
  year_ids <- sort(unique(sp_df$year))

  P <- length(plot_ids)
  Y <- length(year_ids)

  if (ncol(mu) != Y) {
    warning(
      "Number of mu columns does not match number of years in data.\n",
      "ncol(mu) = ", ncol(mu), ", number of year_ids = ", Y, "\n",
      "Check that year indexing in Stan matches sorted unique(data$year)."
    )
  }

  if (ncol(u_plot_mu) != P) {
    warning(
      "Number of u_plot_mu columns does not match number of plots in data.\n",
      "ncol(u_plot_mu) = ", ncol(u_plot_mu), ", number of plot_ids = ", P, "\n",
      "Check that plot indexing in Stan matches sorted unique(data$plot_id)."
    )
  }

  # -----------------------------
  # Storage for summaries
  # -----------------------------
  peak_summary <- list()

  row_id <- 1

  for (y_idx in seq_along(year_ids)) {
    year_val <- year_ids[y_idx]

    for (pl_idx in seq_along(plot_ids)) {
      plot_val <- plot_ids[pl_idx]

      peak_draws <- mu[, y_idx] + u_plot_mu[, pl_idx]

      observed_doy <- sp_df$DOYs[
        sp_df$year == year_val &
          sp_df$plot_id == plot_val
      ]

      observed_abundance <- sp_df$abundance[
        sp_df$year == year_val &
          sp_df$plot_id == plot_val
      ]

      peak_summary[[row_id]] <- data.frame(
        species = species_name,
        year = year_val,
        year_index = y_idx,
        plot_id = plot_val,
        plot_index = pl_idx,

        n_obs = length(observed_doy),
        observed_doy_min = ifelse(length(observed_doy) > 0, min(observed_doy, na.rm = TRUE), NA_real_),
        observed_doy_max = ifelse(length(observed_doy) > 0, max(observed_doy, na.rm = TRUE), NA_real_),
        observed_abundance_max = ifelse(length(observed_abundance) > 0, max(observed_abundance, na.rm = TRUE), NA_real_),

        peak_mean = mean(peak_draws),
        peak_median = median(peak_draws),
        peak_sd = sd(peak_draws),
        peak_q025 = unname(quantile(peak_draws, 0.025)),
        peak_q10 = unname(quantile(peak_draws, 0.10)),
        peak_q90 = unname(quantile(peak_draws, 0.90)),
        peak_q975 = unname(quantile(peak_draws, 0.975)),

        prop_peak_less_than_doy_min = mean(peak_draws < doy_min),
        prop_peak_greater_than_doy_max = mean(peak_draws > doy_max),
        prop_peak_inside_doy_range = mean(peak_draws >= doy_min & peak_draws <= doy_max),

        width_median = median(width[, y_idx]),
        width_q025 = unname(quantile(width[, y_idx], 0.025)),
        width_q975 = unname(quantile(width[, y_idx], 0.975)),

        alpha_year_median = median(alpha_year[, y_idx]),
        u_plot_median = median(u_plot[, pl_idx])
      )

      row_id <- row_id + 1
    }
  }

  peak_summary <- dplyr::bind_rows(peak_summary)

  # -----------------------------
  # Flag suspicious combinations
  # -----------------------------
  peak_summary <- peak_summary |>
    dplyr::mutate(
      peak_problem = dplyr::case_when(
        prop_peak_inside_doy_range < 0.50 ~ "Most posterior peak draws outside DOY range",
        peak_q025 < doy_min & peak_q975 > doy_max ~ "Peak highly uncertain across full DOY range",
        width_median < 0.20 ~ "Very narrow width",
        width_median > 2.00 ~ "Very wide width",
        TRUE ~ "OK"
      )
    )
}
