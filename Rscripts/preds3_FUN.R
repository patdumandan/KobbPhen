plot_plant_preds3 <- function(species_name, data,
                             model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models3\\",
                             out_dir = "predictions3",
                             doy_min = -2,
                             doy_max = 2,
                             n_doy = 100) {

  message("Plotting predictions of: ", species_name)

  # Open model
  model_file <- file.path(model_dir, paste0("phenology_", species_name, ".rds"))

  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }

  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Subset data
  sp_df <- dplyr::filter(data, taxon == species_name)

  if (nrow(sp_df) == 0) {
    stop("No data found for species: ", species_name)
  }

  # Use abundance column if present, otherwise use y
  if ("abundance" %in% names(sp_df)) {
    obs_y <- sp_df$abundance
  } else if ("y" %in% names(sp_df)) {
    obs_y <- sp_df$y
  } else {
    stop("Data must contain either an 'abundance' column or a 'y' column.")
  }

  # -----------------------------
  # Posterior draws
  # -----------------------------

  # Scalar parameters
  alpha_draws <- as.numeric(fit$draws("alpha", format = "draws_matrix"))
  width_draws <- as.numeric(fit$draws("width", format = "draws_matrix"))

  # Indexed parameters
  alpha_year <- as.matrix(fit$draws("alpha_year", format = "draws_matrix"))
  u_plot     <- as.matrix(fit$draws("u_plot", format = "draws_matrix"))
  mu_year    <- as.matrix(fit$draws("mu_year", format = "draws_matrix"))

  # -----------------------------
  # IDs and prediction grid
  # -----------------------------

  DOY <- seq(doy_min, doy_max, length.out = n_doy)

  plot_ids <- sort(unique(sp_df$plot_id))

  if ("year_id" %in% names(sp_df)) {
    year_lookup <- unique(sp_df[, c("year", "year_id")])
    year_lookup <- year_lookup[order(year_lookup$year_id), ]

    year_ids <- year_lookup$year
    year_indices <- year_lookup$year_id
  } else {
    year_ids <- sort(unique(sp_df$year))
    year_indices <- seq_along(year_ids)
  }

  P <- length(plot_ids)
  cols <- rainbow(P)

  # -----------------------------
  # PDF setup
  # -----------------------------

  out_file <- file.path(out_dir, paste0(species_name, "_preds.pdf"))

  pdf(out_file)

  n_years <- length(year_ids)
  par(mfrow = c(ceiling(n_years / 2), 2))

  # -----------------------------
  # Plot by year and plot
  # -----------------------------

  for (yy in seq_along(year_ids)) {
    y_actual <- year_ids[yy]
    y_idx <- year_indices[yy]

    year_obs <- obs_y[sp_df$year == y_actual]

    plot_started <- FALSE

    for (pl in seq_along(plot_ids)) {
      p_id <- plot_ids[pl]

      # If plot_id is already 1:Nplots, this is correct.
      # If your plot_id is not 1:Nplots, you should create a separate plot_index column.
      pl_idx <- p_id

      abundance_hat <- numeric(length(DOY))

      for (nd in seq_along(DOY)) {

        eta_post <-
          alpha_draws +
          alpha_year[, y_idx] +
          u_plot[, pl_idx] -
          ((DOY[nd] - mu_year[, y_idx])^2 / (width_draws^2))

        abundance_post <- exp(eta_post)

        abundance_hat[nd] <- median(abundance_post)
      }

      ind <- sp_df$plot_id == p_id & sp_df$year == y_actual

      if (!plot_started) {
        plot(
          DOY,
          abundance_hat,
          type = "l",
          col = cols[pl],
          main = paste(species_name, "Year", y_actual),
          xlab = "Scaled DOY",
          ylab = "Expected abundance",
          ylim = range(
            0,
            year_obs,
            abundance_hat,
            na.rm = TRUE
          )
        )

        plot_started <- TRUE
      } else {
        lines(DOY, abundance_hat, col = cols[pl])
      }

      points(
        sp_df$DOYs[ind],
        obs_y[ind],
        col = cols[pl],
        pch = 16
      )
    }

    legend(
      "topright",
      legend = plot_ids,
      col = cols,
      lty = 1,
      pch = 16,
      title = "Plot",
      cex = 0.7,
      bty = "n"
    )
  }

  dev.off()
}
