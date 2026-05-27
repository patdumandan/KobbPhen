plot_plant_preds2 <- function(species_name, data,
                             model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models2\\",
                             out_dir = "predictions2") {

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

  # Posterior draws
  # Keep vector parameters as numeric vectors
  alpha_draws <- as.numeric(fit$draws("alpha", format = "draws_matrix"))
  mu_draws    <- as.numeric(fit$draws("mu", format = "draws_matrix"))
  width_draws <- as.numeric(fit$draws("width", format = "draws_matrix"))

  # Keep indexed parameters as matrices
  alpha_year <- fit$draws("alpha_year", format = "draws_matrix")
  u_plot     <- fit$draws("u_plot", format = "draws_matrix")

  # Prediction grid on standardized DOY scale
  DOY <- seq(-2, 2, length.out = 100)

  # PDF
  pdf(file.path(out_dir, paste0(species_name, "_preds.pdf")))

  plot_ids <- sort(unique(sp_df$plot_id))
  year_ids <- sort(unique(sp_df$year))

  P <- length(plot_ids)
  cols <- rainbow(P)

  n_years <- length(year_ids)
  par(mfrow = c(ceiling(n_years / 2), 2))

  for (y_idx in seq_along(year_ids)) {
    y <- year_ids[y_idx]

    for (pl_idx in seq_along(plot_ids)) {
      p_id <- plot_ids[pl_idx]

      abundance_hat <- numeric(length(DOY))

      for (nd in seq_along(DOY)) {
        eta_post <-
          alpha_draws +
          alpha_year[, y_idx] +
          u_plot[, pl_idx] -
          ((DOY[nd] - mu_draws)^2 / (width_draws^2))

        abundance_post <- exp(eta_post)

        abundance_hat[nd] <- median(abundance_post)
      }

      ind <- sp_df$plot_id == p_id & sp_df$year == y

      if (pl_idx == 1) {
        plot(
          DOY,
          abundance_hat,
          type = "l",
          col = cols[pl_idx],
          main = paste(species_name, "Year", y),
          xlab = "Scaled DOY",
          ylab = "Expected abundance",
          ylim = range(
            0,
            sp_df$abundance[ind],
            abundance_hat,
            na.rm = TRUE
          )
        )
      } else {
        lines(DOY, abundance_hat, col = cols[pl_idx])
      }

      points(
        sp_df$DOYs[ind],
        sp_df$abundance[ind],
        col = cols[pl_idx],
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
