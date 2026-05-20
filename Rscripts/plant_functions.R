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
    iter_sampling = 200,
    iter_warmup = 50,
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

        eta_hat[nd] <- mean(eta_post)
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

