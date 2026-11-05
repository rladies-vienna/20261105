# functions for building the mesh and plot the variogram
library(fields)
library(spam)

################# MESH FUNCTIONS ##################
estimate_spatial_range <- function(coords, response, unit = "units") {
  
  # data frame from inputs, gstat reads the coordinates from it (no sp object needed)
  df <- data.frame(x = coords[, 1], y = coords[, 2], z = response)
  
  # empirical variogram of the response (constant mean)
  vgm_emp <- variogram(z ~ 1, locations = ~ x + y, data = df)
  
  # fit both Exponential (nu = 0.5) and Matérn (nu = 1) to get direct estimates
  fit_exp <- fit.variogram(vgm_emp, model = vgm("Exp"))
  a_exp   <- fit_exp$range[2]                       # gstat scale a of the Exp fit
  
  # Matérn start: same practical range as the Exp fit (2 a_exp = sqrt(8) a_mat)
  a_mat_start <- 2 * a_exp / sqrt(8)
  init_mat <- vgm(psill = var(df$z) * 0.9, model = "Mat", range = a_mat_start,
                  kappa = 1, nugget = fit_exp$psill[1])
  fit_mat  <- fit.variogram(vgm_emp, model = init_mat, fit.kappa = FALSE)
  
  # gstat scale a -> INLA practical range rho = sqrt(8 * nu) * a
  rho_exp_est <- fit_exp$range[2] * sqrt(8 * 0.5)
  rho_mat_est <- fit_mat$range[2] * sqrt(8 * 1)
  
  # sill metrics: row 1 = nugget, row 2 = structured part
  nugget_exp <- fit_exp$psill[1]
  nugget_mat <- fit_mat$psill[1]
  psill_exp  <- fit_exp$psill[2]
  psill_mat  <- fit_mat$psill[2]
  
  # partial sill -> INLA marginal sd of the field
  sigma_exp_est <- sqrt(psill_exp)
  sigma_mat_est <- sqrt(psill_mat)
  
  cat("\n--- DETECTED EMPIRICAL ESTIMATES (INLA METRICS) ---\n")
  cat(sprintf("Exponential (nu=0.5) -> Implied INLA rho: %.3g %s | sigma: %.3f\n",
              rho_exp_est, unit, sigma_exp_est))
  cat(sprintf("Matérn      (nu=1.0) -> Implied INLA rho: %.3g %s | sigma: %.3f\n",
              rho_mat_est, unit, sigma_mat_est))
  
  list(
    vgm_emp    = vgm_emp,
    fit_exp    = fit_exp,
    fit_mat    = fit_mat,
    nugget_exp = nugget_exp,
    nugget_mat = nugget_mat,
    psill_exp  = psill_exp,
    psill_mat  = psill_mat,
    rho_exp    = rho_exp_est,
    rho_mat    = rho_mat_est,
    sigma_exp  = sigma_exp_est,
    sigma_mat  = sigma_mat_est,
    unit       = unit)
}


build_mesh_ensemble <- function(coords, boundary, response,
                                n_levels = 6, unit = "units") {
  
  est       <- estimate_spatial_range(coords, response, unit = unit)
  emp_range <- mean(c(est$rho_exp, est$rho_mat))
  
  inner_edge_vals <- seq(emp_range / 3, emp_range, length.out = n_levels)
  
  # signif() instead of round(): normalised values are below 1
  cat(sprintf("Effective range: %.3g %s\n", emp_range, unit))
  cat("Inner edge levels:", signif(inner_edge_vals / 2, 3), unit, "\n")
  
  # build the meshes; all settings are fractions of the range, so no unit is fixed here
  mesh_list <- lapply(inner_edge_vals, function(inner) {
    fm_mesh_2d(
      #loc.domain = coords,    #your coords if you do not set boundary (but then only rough polygon)
      boundary = boundary,     #sets the inner domain, the outer ring is built from offset[2]
      #Inner edges run from range/6 to range/2. The usual rule is range/5 to range/10,
      #so only your finest level meets it.
      max.edge = c(inner / 2, inner * 1.5),   # longest triangle edge inside and in the outer ring
      cutoff   = inner / 5,                   # merges input points closer than this
      offset   = c(inner, emp_range * 2)      # width of the inner and outer extension
    )
  })
  
  cat("Built", length(mesh_list), "meshes | nodes:",
      sapply(mesh_list, function(m) m$n), "\n")
  
  list(mesh_list = mesh_list, inner_edge_vals = inner_edge_vals,
       emp_range = emp_range, est = est)
}


################# VARIOGRAM PLOTS ##################

plot_variogram_fit <- function(est, unit = "m",
                               ann_x = c(0.45, 0.75),  # box positions, fraction of max distance
                               ann_y = 0.05) {         # box bottom, fraction of max semivariance
  
  ## gstat range a = rho / sqrt(8 * nu)
  a_exp <- est$rho_exp / sqrt(8 * 0.5)
  a_mat <- est$rho_mat / sqrt(8 * 1)
  
  ## rebuild the two fitted models from the returned values
  fit_exp <- vgm(psill = est$psill_exp, model = "Exp", range = a_exp,
                 nugget = est$nugget_exp)
  fit_mat <- vgm(psill = est$psill_mat, model = "Mat", range = a_mat,
                 nugget = est$nugget_mat, kappa = 1)
  
  h    <- seq(0, max(est$vgm_emp$dist), length.out = 300)
  mods <- c("Exponential (nu = 0.5)", "Matern (nu = 1)")
  curves <- rbind(
    data.frame(dist = h, gamma = variogramLine(fit_exp, dist_vector = h)$gamma, model = mods[1]),
    data.frame(dist = h, gamma = variogramLine(fit_mat, dist_vector = h)$gamma, model = mods[2]))
  
  ## text blocks: gstat scale, INLA range, partial sill, INLA sd, nugget
  txt <- function(name, a, rho, psill, sigma, nugget) sprintf(
    "%s\ngstat scale a = %.0f %s\nINLA range rho = %.0f %s\npartial sill = %.4f\nINLA sd sigma = %.3f\nnugget = %.4f",
    name, a, unit, rho, unit, psill, sigma, nugget)
  
  ymax <- max(c(est$vgm_emp$gamma, curves$gamma))
  ann  <- data.frame(
    model = mods,
    dist  = max(h) * ann_x,
    gamma = ymax * ann_y,
    label = c(txt(mods[1], a_exp, est$rho_exp, est$psill_exp, est$sigma_exp, est$nugget_exp),
              txt(mods[2], a_mat, est$rho_mat, est$psill_mat, est$sigma_mat, est$nugget_mat)))
  
  ggplot(as.data.frame(est$vgm_emp), aes(dist, gamma)) +
    geom_point(aes(size = np), colour = "#13262F", alpha = 0.8) +
    geom_line(data = curves, aes(colour = model), linewidth = 1) +
    geom_label(data = ann, aes(label = label, colour = model),
               hjust = 0, vjust = 0, size = 3.6, lineheight = 1.1,
               fill = "white", show.legend = FALSE) +
    scale_colour_manual(values = c("#2F6F8F", "#D9772B")) +
    scale_size_continuous(range = c(1.5, 4)) +
    scale_y_continuous(limits = c(0, ymax * 1.05)) +
    labs(x = paste0("distance (", unit, ")"), y = "semivariance",
         colour = NULL, size = "pairs",
         title = "Empirical variogram and fitted models") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}


######################## A MATRIX PLOTS ###################
plot_A_row <- function(mesh, A, coords, i) {
  a   <- A[i, ]
  idx <- which(a > 0)                      # mesh nodes used by station i
  psi <- a[idx]                            # psi_k(s_i) = A[i, k], the hat function values
  v   <- mesh$loc[idx, 1:2, drop = FALSE]  # node coordinates
  s <- as.numeric(coords[i, ])
  
  ## neighbouring triangles for context
  tv  <- mesh$graph$tv
  nb  <- tv[rowSums(matrix(tv %in% idx, ncol = 3)) > 0, , drop = FALSE]
  e   <- rbind(nb[, 1:2], nb[, 2:3], nb[, c(3, 1)])
  edges <- data.frame(x = mesh$loc[e[, 1], 1], y = mesh$loc[e[, 1], 2],
                      xend = mesh$loc[e[, 2], 1], yend = mesh$loc[e[, 2], 2])
  
  tri   <- data.frame(x = v[, 1], y = v[, 2])
  spoke <- data.frame(x = s[1], y = s[2], xend = v[, 1], yend = v[, 2])
  lab   <- data.frame(x = v[, 1], y = v[, 2],
                      label = sprintf("node %d\n\u03c8 = %.2f", idx, psi))
  
  ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
                 colour = "grey80", linewidth = 0.4) +
    geom_polygon(data = tri, aes(x, y), fill = "#E4EEF2",
                 colour = "#2F6F8F", linewidth = 1.2) +
    geom_segment(data = spoke, aes(x, y, xend = xend, yend = yend),
                 colour = "#D9772B", linetype = 2, linewidth = 0.9) +
    geom_point(data = tri, aes(x, y), colour = "#13262F", size = 5) +
    geom_point(aes(x = s[1], y = s[2]), colour = "#D9772B", size = 5) +
    geom_label(data = lab, aes(x, y, label = label), size = 3.5,
               vjust = -0.4, label.size = 0, fill = "white") +
    coord_equal() +
    labs(title = sprintf("Station %d = row %d of A", i, i),
         subtitle = sprintf("u(s%d) = %s   (sum = %.2f)", i,
                            paste(sprintf("%.2f w%d", psi, idx), collapse = " + "),
                            sum(psi)),
         x = "x", y = "y") +
    theme_minimal(base_size = 13)
}

######################## SCORING FUNCTION ###################
#Predictive distribution of a new measurement: N(mu, sd^2), sd includes the noise

score_gauss <- function(y, mu, sd, dens = NULL) {
  z  <- (y - mu) / sd
  ls <- if (is.null(dens)) -dnorm(y, mu, sd, log = TRUE) else -log(dens)
  crps <- sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
  c(n        = length(y),
    logscore = mean(ls),                       # lower = better
    RMSE     = sqrt(mean((y - mu)^2)),
    MAE      = mean(abs(y - mu)),
    CRPS     = mean(crps),                     # lower = better, log10 units
    cover95  = mean(abs(z) <= qnorm(0.975)))   # should be close to 0.95
}



# for inla.group.cv output. $mean and $sd are for the linear predictor,
# group_size counts stations only (no grid rows).
cv_scores <- function(g) c(
  score_gauss(y_obs, g$mean[idx_est], g$sd[idx_est]),     # sd of eta only, no noise
  group_size = mean(sapply(g$groups[seq_along(idx_est)],
                           function(x) sum(x$idx %in% idx_est))))


## One row of the final table. Scores a method does not give stay NA.
as_row <- function(x, method, rad = NA) {
  cols <- c("n", "group_size", "logscore", "RMSE", "MAE", "CRPS", "cover95")
  out  <- setNames(rep(NA_real_, length(cols)), cols)
  hit  <- intersect(names(x), cols)
  out[hit] <- x[hit]
  data.frame(method = method, rad_m = rad, t(out))
}


plot_cv_group <- function(i, removed, rad = NULL, title = NULL) {
  pts <- data.frame(x = coords_stations[, 1], y = coords_stations[, 2], role = "kept")
  pts$role[removed] <- "left out"
  pts$role[i]       <- "test station"
  p <- ggplot(pts, aes(x, y, colour = role)) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c(kept = "grey75", `left out` = "#2F6F8F",
                                   `test station` = "#D9772B"))
  if (!is.null(rad) && rad > 0) {
    a    <- seq(0, 2 * pi, length.out = 200)
    circ <- data.frame(x = coords_stations[i, 1] + rad * cos(a),
                       y = coords_stations[i, 2] + rad * sin(a))
    p <- p + geom_path(data = circ, aes(x, y), inherit.aes = FALSE,
                       colour = "#D9772B", linetype = 2)
  }
  p + coord_equal() +
    labs(title = title, x = "x", y = "y", colour = NULL) +
    theme_minimal(base_size = 13)
}




