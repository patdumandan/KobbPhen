plot_taxon_community_preds_by_plot <- function(species_names, data,
                                               model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen\\models4\\",
                                               out_dir = "predictions4",
                                               out_file = "taxon_community_preds_by_plot.pdf",
                                               doy_min = -2,
                                               doy_max = 2,
                                               n_doy = 100) {

  message("Plotting taxon-specific and community predictions by plot.")

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

    # If a species model does not have this plot, return NA curve
    if (pl_idx > ncol(draws$u_plot)) {
      return(rep(NA_real_, length(DOY)))
    }

    if (y_idx > ncol(draws$alpha_year) ||
        y_idx > ncol(draws$mu_year) ||
        y_idx > ncol(draws$width_year)) {
      return(rep(NA_real_, length(DOY)))
    }

    u_plot_year_name <- paste0("u_plot_year[", pl_idx, ",", y_idx, "]")

    if (!u_plot_year_name %in% colnames(draws$u_plot_year)) {
      return(rep(NA_real_, length(DOY)))
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

  year_cols <- grDevices::colorRampPalette(
    c("#c6dbef", "#08306b")
  )(length(year_ids))

  names(year_cols) <- year_ids

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

  if (length(fits) == 0) {
    stop("No fitted models were loaded.")
  }

  # -----------------------------
  # PDF setup
  # -----------------------------

  pdf_file <- file.path(out_dir, out_file)

  pdf(pdf_file)
  on.exit(dev.off(), add = TRUE)

  n_panels <- length(species_names) + 1
  n_cols <- 2
  n_rows <- ceiling(n_panels / n_cols)

  # -----------------------------
  # Plot one page per plot_id
  # -----------------------------

  for (pl in seq_along(plot_ids)) {
    p_id <- plot_ids[pl]

    # If plot_id is already 1:Nplots, this is correct.
    # If your model uses a different plot index, replace this with that index.
    pl_idx <- p_id

    message("Plotting plot_id: ", p_id)

    par(mfrow = c(n_rows, n_cols))
    par(mar = c(4, 4, 3, 1))

    # Store curves so we can sum them for the community panel
    all_curves <- list()

    # -----------------------------
    # Species panels
    # -----------------------------

    for (sp in species_names) {

      species_curves <- list()

      for (year_actual in year_ids) {
        if (!sp %in% names(fits)) {
          species_curves[[as.character(year_actual)]] <- rep(NA_real_, length(DOY))
          next
        }

        species_curves[[as.character(year_actual)]] <- predict_curve(
          draws = fits[[sp]],
          DOY = DOY,
          pl_idx = pl_idx,
          year_actual = year_actual
        )
      }

      all_curves[[sp]] <- species_curves

      curve_values <- unlist(species_curves)
      curve_values <- curve_values[is.finite(curve_values)]

      obs_ind_species <-
        plot_data$taxon == sp &
        plot_data$plot_id == p_id

      ylim_vals <- range(
        0,
        plot_data$obs_y[obs_ind_species],
        curve_values,
        na.rm = TRUE,
        finite = TRUE
      )

      if (!all(is.finite(ylim_vals))) {
        ylim_vals <- c(0, 1)
      }

      plot(
        DOY,
        rep(NA_real_, length(DOY)),
        type = "n",
        main = paste(sp, "- Plot", p_id),
        xlab = "Scaled DOY",
        ylab = "Expected abundance",
        ylim = ylim_vals
      )

      for (year_actual in year_ids) {
        year_chr <- as.character(year_actual)
        curve <- species_curves[[year_chr]]

        if (!all(is.na(curve))) {
          lines(
            DOY,
            curve,
            col = year_cols[year_chr],
            lwd = 2
          )
        }

        obs_ind <-
          plot_data$taxon == sp &
          plot_data$plot_id == p_id &
          plot_data$year == year_actual

        points(
          plot_data$DOYs[obs_ind],
          plot_data$obs_y[obs_ind],
          col = year_cols[year_chr],
          pch = 16
        )
      }

      legend(
        "topright",
        legend = year_ids,
        col = year_cols,
        lty = 1,
        lwd = 2,
        pch = 16,
        title = "Year",
        cex = 0.65,
        bty = "n"
      )
    }

    # -----------------------------
    # Community total panel
    # -----------------------------

    community_curves <- list()

    for (year_actual in year_ids) {
      year_chr <- as.character(year_actual)

      curves_this_year <- list()

      for (sp in species_names) {
        if (!sp %in% names(all_curves)) {
          next
        }

        curve <- all_curves[[sp]][[year_chr]]

        if (is.null(curve) || all(is.na(curve))) {
          next
        }

        curves_this_year[[sp]] <- curve
      }

      if (length(curves_this_year) == 0) {
        community_curves[[year_chr]] <- rep(NA_real_, length(DOY))
      } else {
        community_curves[[year_chr]] <- Reduce("+", curves_this_year)
      }
    }

    community_values <- unlist(community_curves)
    community_values <- community_values[is.finite(community_values)]

    community_obs <- plot_data[
      plot_data$plot_id == p_id,
      c("year", "DOYs", "obs_y")
    ]

    if (nrow(community_obs) > 0) {
      community_obs_sum <- stats::aggregate(
        obs_y ~ year + DOYs,
        data = community_obs,
        FUN = sum,
        na.rm = TRUE
      )
    } else {
      community_obs_sum <- data.frame(
        year = numeric(),
        DOYs = numeric(),
        obs_y = numeric()
      )
    }

    ylim_vals <- range(
      0,
      community_obs_sum$obs_y,
      community_values,
      na.rm = TRUE,
      finite = TRUE
    )

    if (!all(is.finite(ylim_vals))) {
      ylim_vals <- c(0, 1)
    }

    plot(
      DOY,
      rep(NA_real_, length(DOY)),
      type = "n",
      main = paste("Community total - Plot", p_id),
      xlab = "Scaled DOY",
      ylab = "Summed expected abundance",
      ylim = ylim_vals
    )

    for (year_actual in year_ids) {
      year_chr <- as.character(year_actual)
      curve <- community_curves[[year_chr]]

      if (!all(is.na(curve))) {
        lines(
          DOY,
          curve,
          col = year_cols[year_chr],
          lwd = 3
        )
      }

      obs_ind <- community_obs_sum$year == year_actual

      points(
        community_obs_sum$DOYs[obs_ind],
        community_obs_sum$obs_y[obs_ind],
        col = year_cols[year_chr],
        pch = 16
      )
    }

    legend(
      "topright",
      legend = year_ids,
      col = year_cols,
      lty = 1,
      lwd = 2,
      pch = 16,
      title = "Year",
      cex = 0.65,
      bty = "n"
    )
  }

  message("Saved plot to: ", pdf_file)

  invisible(pdf_file)
}

plot_taxon_community_preds_by_plot(
  species_names = c(
    "Betula_nana",
    "Rhododendron_groenlandicum",
    "Salix_glauca",
    "Salix_herbacea",
    "Vaccinium_uliginosum"
  ),
  data = kobb_dat2
)
