plot_taxon_preds_by_plot <- function(species_names, data,
                                     model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models4\\",
                                     out_dir = "predictions4",
                                     out_file = "taxon_preds_by_plot.pdf",
                                     doy_min = -2,
                                     doy_max = 2,
                                     n_doy = 100) {

  message("Plotting taxon-specific predictions by plot.")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # -----------------------------
  # Helper: load model draws
  # -----------------------------

  load_fit_draws <- function(species_name, sp_df) {
    model_file <- file.path(model_dir, paste0("phenology_", species_name, ".rds"))

    if (!file.exists(model_file)) {
      stop("Model file not found: ", model_file)
    }

    fit <- readRDS(model_file)

    species_years <- sort(unique(sp_df$year))

    list(
      species_name = species_name,
      species_years = species_years,

      alpha_draws = as.numeric(fit$draws("alpha", format = "draws_matrix")),

      alpha_year = as.matrix(fit$draws("alpha_year", format = "draws_matrix")),
      u_plot = as.matrix(fit$draws("u_plot", format = "draws_matrix")),
      mu_year = as.matrix(fit$draws("mu_year", format = "draws_matrix")),
      width_year = as.matrix(fit$draws("width_year", format = "draws_matrix")),
      u_plot_year = as.matrix(fit$draws("u_plot_year", format = "draws_matrix"))
    )
  }

  # -----------------------------
  # Helper: predict curve
  # -----------------------------

  predict_curve <- function(draws, DOY, pl_idx, year_actual) {

    # Species-specific year index
    y_idx <- match(year_actual, draws$species_years)

    if (is.na(y_idx)) {
      return(rep(NA_real_, length(DOY)))
    }

    u_plot_year_name <- paste0("u_plot_year[", pl_idx, ",", y_idx, "]")

    if (!u_plot_year_name %in% colnames(draws$u_plot_year)) {
      stop(
        "Could not find column in posterior draws: ", u_plot_year_name,
        "\nSpecies: ", draws$species_name,
        "\nYear: ", year_actual,
        "\nAvailable u_plot_year columns include:\n",
        paste(head(colnames(draws$u_plot_year), 20), collapse = ", ")
      )
    }

    u_plot_year_draws <- draws$u_plot_year[, u_plot_year_name]

    abundance_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {
      eta_post <-
        draws$alpha_draws +
        draws$alpha_year[, y_idx] +
        draws$u_plot[, pl_idx] +
        u_plot_year_draws -
        ((DOY[nd] - draws$mu_year[, y_idx])^2 /
           (draws$width_year[, y_idx]^2))

      abundance_post <- exp(eta_post)

      abundance_hat[nd] <- median(abundance_post, na.rm = TRUE)
    }

    abundance_hat
  }

  # -----------------------------
  # Subset data
  # -----------------------------

  plot_data <- dplyr::filter(data, taxon %in% species_names)

  if (nrow(plot_data) == 0) {
    stop("No data found for requested species.")
  }

  if ("abundance" %in% names(plot_data)) {
    plot_data$obs_y <- plot_data$abundance
  } else if ("y" %in% names(plot_data)) {
    plot_data$obs_y <- plot_data$y
  } else {
    stop("Data must contain either an 'abundance' column or a 'y' column.")
  }

  required_cols <- c("taxon", "plot_id", "year", "DOYs", "obs_y")

  if (!all(required_cols %in% names(plot_data))) {
    stop("Data must contain columns: taxon, plot_id, year, DOYs, and abundance or y.")
  }

  DOY <- seq(doy_min, doy_max, length.out = n_doy)

  plot_ids <- sort(unique(plot_data$plot_id))
  year_ids <- sort(unique(plot_data$year))

  taxon_cols <- rainbow(length(species_names))
  names(taxon_cols) <- species_names

  # -----------------------------
  # Load all model draws
  # -----------------------------

  fits <- list()

  for (sp in species_names) {
    message("Loading model for: ", sp)

    sp_df <- dplyr::filter(plot_data, taxon == sp)

    if (nrow(sp_df) == 0) {
      warning("No data found for species: ", sp)
      next
    }

    fits[[sp]] <- load_fit_draws(sp, sp_df)
  }

  # -----------------------------
  # PDF setup
  # -----------------------------

  pdf_file <- file.path(out_dir, out_file)

  pdf(pdf_file)
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(ceiling(length(plot_ids) / 2), 2))

  # -----------------------------
  # Plot by year, then plot_id
  # -----------------------------

  for (year_actual in year_ids) {
    message("Plotting year: ", year_actual)

    for (pl in seq_along(plot_ids)) {
      p_id <- plot_ids[pl]
      pl_idx <- p_id

      curves <- list()

      for (sp in species_names) {
        if (!sp %in% names(fits)) {
          next
        }

        curves[[sp]] <- predict_curve(
          draws = fits[[sp]],
          DOY = DOY,
          pl_idx = pl_idx,
          year_actual = year_actual
        )
      }

      obs_ind <- plot_data$plot_id == p_id & plot_data$year == year_actual

      curve_values <- unlist(curves)
      curve_values <- curve_values[is.finite(curve_values)]

      ylim_vals <- range(
        0,
        plot_data$obs_y[obs_ind],
        curve_values,
        na.rm = TRUE,
        finite = TRUE
      )

      plot(
        DOY,
        rep(NA_real_, length(DOY)),
        type = "n",
        main = paste("Plot", p_id, "Year", year_actual),
        xlab = "Scaled DOY",
        ylab = "Expected abundance",
        ylim = ylim_vals
      )

      for (sp in species_names) {
        if (!sp %in% names(curves)) {
          next
        }

        if (all(is.na(curves[[sp]]))) {
          next
        }

        lines(
          DOY,
          curves[[sp]],
          col = taxon_cols[sp],
          lwd = 2
        )

        sp_obs_ind <-
          plot_data$taxon == sp &
          plot_data$plot_id == p_id &
          plot_data$year == year_actual

        points(
          plot_data$DOYs[sp_obs_ind],
          plot_data$obs_y[sp_obs_ind],
          col = taxon_cols[sp],
          pch = 16
        )
      }

      legend(
        "topright",
        legend = species_names,
        col = taxon_cols,
        lty = 1,
        lwd = 2,
        pch = 16,
        title = "Taxon",
        cex = 0.65,
        bty = "n"
      )
    }
  }

  message("Saved plot to: ", pdf_file)

  invisible(pdf_file)
}

plot_taxon_preds_by_plot(
  species_names = c(
    "Betula_nana",
    "Rhododendron_groenlandicum",
    "Salix_glauca",
    "Salix_herbacea",
    "Vaccinium_uliginosum"
  ),
  data = kobb_dat2
)
