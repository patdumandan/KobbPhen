plot_plant_preds1 <- function(species_name, data,
                             model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models1\\",
                             out_dir = "predictions1") {

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

  # alpha_draws <- as.numeric(fit$draws("alpha", format = "draws_matrix"))
  # mu_draws    <- as.numeric(fit$draws("mu", format = "draws_matrix"))
  # width_draws <- as.numeric(fit$draws("width", format = "draws_matrix"))
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
         eta_post <- alpha_year[, y_idx] + u_plot[, pl_idx] -
                           ((DOY[nd] - (mu[, y_idx] + u_plot_mu[, pl_idx]))^2 / (width[, y_idx]^2))
        # eta_post <-
        #   alpha_draws +
        #   alpha_year[, y_idx] +
        #   u_plot[, pl_idx] -
        #   ((DOY[nd] - mu_draws)^2 / (width_draws^2))

        abundance_post <- exp(eta_post)
        eta_hat[nd] <- median(abundance_post)
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
