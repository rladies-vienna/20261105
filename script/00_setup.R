## =============================================================================
## 00_setup.R  |  Run this once BEFORE the tutorial
## =============================================================================
## INLA is not on CRAN. It comes from its own repository.
## Everything else comes from CRAN, except INLAutils (GitHub).

## ---- 1. Install what is missing --------------------------------------------
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(
    "INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/stable"),
    dependencies = TRUE
  )
}

cran_pkgs <- c("fmesher", "sp", "ggplot2", "patchwork", "remotes", "lubridate",
               "raster", "cowplot", "reshape2", "tidyr", "assertthat", "dplyr", "sf",
               "RColorBrewer", "scales", "gstat", "terra", "patchwork","fields",
               "Matrix", "usethis", "devtools", "spam"
               )


missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1),
                             quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

if (!requireNamespace("INLAutils", quietly = TRUE)) {
  remotes::install_github("timcdlucas/INLAutils", upgrade = "never")
}

## ---- 2. Load packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(INLA)      # model fitting
  library(fmesher)   # meshes and projection matrices
  library(ggplot2)   # plotting
  library(patchwork) # combine ggplots with + and |
})

cat("INLA   ", as.character(packageVersion("INLA")), "\n")
cat("fmesher", as.character(packageVersion("fmesher")), "\n")

## ---- 3. Smoke test: can R talk to the inla binary? ---------------------------
smoke <- inla(y ~ 1, data = data.frame(y = rnorm(20)))
stopifnot(inherits(smoke, "inla"))
cat("INLA smoke test OK\n")

## ---- 4. Data: Meuse river zinc ---------------------------------------------
## 155 topsoil samples along the Meuse (NL). Coordinates are in metres
## (Dutch RD grid). We rescale to km: priors and mesh settings read easier,
## and the numerics behave better with smaller coordinate values.
data(meuse, package = "sp")       # observations
data(meuse.grid, package = "sp")  # 3103 prediction cells, 40 m x 40 m

obs <- data.frame(
  x        = meuse$x / 1000,   # km
  y        = meuse$y / 1000,   # km
  log_zinc = log(meuse$zinc),  # response: log zinc concentration (ppm)
  dist     = meuse$dist        # covariate: scaled distance to the river [0, 1]
)
grid <- data.frame(
  x    = meuse.grid$x / 1000,
  y    = meuse.grid$y / 1000,
  dist = meuse.grid$dist
)

coords      <- as.matrix(obs[, c("x", "y")])   # n x 2 matrix, needed by fmesher
grid_coords <- as.matrix(grid[, c("x", "y")])

## A first look at the data
p_data <- ggplot(obs, aes(x, y, colour = log_zinc)) +
  geom_point(size = 2) +
  scale_colour_viridis_c() +
  coord_equal() +
  labs(x = "Easting (km)", y = "Northing (km)", colour = "log(zinc)",
       title = "Point data: what does the surface look like in between?")
print(p_data)
