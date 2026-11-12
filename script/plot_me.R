library(ggplot2)
library(scales)
library(patchwork)

load("113AK_tutorial_data.RData")

################################ data plot ##############
# scale the legend for the plot
brks <- c(0, 25, 50, 100, 150, 200, 300, 500, 600)
labs <- c("0-25", "25-50", "50-100", "100-150", "150-200", "200-300", "300-500", "500-600")
stations$concentration_cat <- cut(stations$concentration, breaks = brks, labels = labs,
                                  include.lowest = TRUE)

sulfate_palette <- c(
  "0-25"    = "#000000", "25-50"   = "#7F7F7F", "50-100"  = "#6BAED6",
  "100-150" = "#2171B5", "150-200" = "#74C476", "200-300" = "#E31A1C",
  "300-500" = "#CC0033", "500-600" = "#FF9800")

p_stations <- ggplot() +
  # domain: prediction grid in light grey, aquifer outline on top
  geom_tile(data = pred_grid, aes(x, y), fill = "lightgrey") +
  geom_sf(data = boundary, fill = NA, colour = "grey40") +
  # quantified stations: filled points coloured by category
  geom_point(data = subset(stations, !censored),
             aes(x, y, fill = concentration_cat),
             shape = 21, size = 3, stroke = 0.3, colour = "black") +
  # censored stations: crosses, with their own legend entry
  geom_point(data = subset(stations, censored),
             aes(x, y, shape = "below detection"),
             size = 3, stroke = 1, colour = "black") +
  scale_fill_manual(
    values = sulfate_palette,
    name   = "Sulfate levels (mg/L)",
    drop   = FALSE,
    guide  = guide_legend(title.position = "top", order = 1,
                          override.aes = list(size = 4))) +
  scale_shape_manual(values = c("below detection" = 3), name = NULL,
                     guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = seq(-1, 2, by = 0.2)) +
  scale_y_continuous(breaks = seq(-1, 2, by = 0.2)) +
  theme_minimal() +
  theme(
    axis.title        = element_text(size = 16, face = "bold"),
    axis.text         = element_text(size = 16, face = "bold", colour = "black"),
    legend.text       = element_text(size = 16, face = "bold"),
    legend.title      = element_text(size = 16, face = "bold"),
    legend.key.height = unit(0.8, "cm"),
    plot.title        = element_text(size = 22, face = "bold")) +
  labs(x = "x", y = "y")

p_stations


################################ covariate plots ##############
theme_norm <- theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())

colourbar <- function() guide_colorbar(
  barheight = unit(6, "cm"), barwidth = unit(0.5, "cm"),
  ticks.colour = "black", frame.colour = "black", title.position = "top")

## One map per covariate: raster, boundary, colour scale, title
plot_covariate <- function(var, fill_scale, title) {
  ggplot() +
    geom_raster(data = pred_grid, aes(x, y, fill = .data[[var]])) +
    geom_sf(data = boundary, fill = NA, colour = "black") +
    fill_scale +
    labs(x = "x ", y = "y ", title = title) +
    theme_norm
}

# altitude plot
my_terrain_colors <- c("#1B5E20", "#4CAF50", "#FFEB3B", "#FFC107", "#8D6E63", "#CFD8DC")

p_altitude <- plot_covariate("altitude", scale_fill_gradientn(
  colors = my_terrain_colors,
  breaks = seq(0, 300, by = 50),
  labels = paste0(seq(0, 300, by = 50), " m"),
  name   = "Altitude",
  # map colours to quantiles of the data
  values = rescale(quantile(pred_grid$altitude, c(0, 0.2, 0.5, 0.7, 0.9, 1), na.rm = TRUE)),
  guide  = colourbar()),
  title = "Altitude")

# IDPR plot
idpr_colors <- c("#08306B",                          # very high infiltration
                 "#2171B5", "#6BAED6", "#F7F7F7",
                 "#FEB24C", "#E31A1C",
                 "#800026")                          # very high runoff

p_idpr <- plot_covariate("IDPR", scale_fill_gradientn(
  colors = idpr_colors,
  name   = "IDPR",
  breaks = seq(0, 2000, by = 250),
  guide  = colourbar()),
  title = "IDPR")

# distance to river plot
river_colors <- c("#084594", "#2171B5", "#6BAED6", "#C6DBEF", "#F7F7F7")

p_river <- plot_covariate("river_dist", scale_fill_gradientn(
  colors = river_colors,
  name   = "Distance to river",
  breaks = seq(0, 15000, by = 2500),
  labels = paste0(seq(0, 15000, by = 2500), " m"),
  guide  = colourbar()),
  title = "Distance to river")

# distance to fault plot
fault_colors <- c("#D500F9",                         # high structural interest
                  "#651FFF", "#311B92",
                  "#424242", "#212121")              # far

p_fault <- plot_covariate("fault_dist", scale_fill_gradientn(
  colors = fault_colors,
  name   = "Distance to fault",
  breaks = seq(0, 73000, by = 10000),
  labels = paste0(seq(0, 73000, by = 10000), " m"),
  guide  = colourbar()),
  title = "Distance to fault")

# distance to gypsum plot
gypsum_colors <- c("#BF360C",                        # the source
                   "#E65100", "#FFB74D",
                   "#FFF9C4", "#EEEEEE")             # far

p_gypsum <- plot_covariate("dist_to_gypse", scale_fill_gradientn(
  colors = gypsum_colors,
  name   = "Distance to gypsum",
  breaks = seq(0, 110000, by = 20000),
  labels = paste0(seq(0, 110000, by = 20000), " m"),
  guide  = colourbar()),
  title = "Distance to gypsum")

# number of depth layers plot
pred_grid$order_discrete <- factor(pmin(11, pmax(1, pred_grid$order_min)), levels = 1:11)

# named palette: layer k always gets colour k, even if some layers are missing
depth_palette <- setNames(c("#FFF9C4", "#FFF176", "#FFD54F", "#FFB300", "#FF8F00",
                            "#EF6C00", "#D84315", "#A1887F", "#6D4C41", "#3E2723", "black"),
                          1:11)

p_depth <- plot_covariate("order_discrete", scale_fill_manual(
  values = depth_palette,
  name   = "Depth layer",
  guide  = guide_legend(ncol = 2, title.position = "top")),
  title = "Number of depth layers")

p_altitude+p_idpr
p_river+p_depth
p_fault/p_gypsum


## All six on one page
# (p_altitude | p_idpr | p_depth) / (p_river | p_fault | p_gypsum)
