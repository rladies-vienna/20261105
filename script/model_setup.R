library(gstat)
library(patchwork)
library(INLA)
library(Matrix)
library(INLAutils)
library(ggplot2)
library(sf)
library(sp) #for coordinates
library(fmesher)

source("functions.R")
load("113AK_tutorial_data.RData")   # stations, pred_grid, boundary, covar_names

## All coordinates are normalised: one shift and one factor for both axes.
## The longer side of the stations' bounding box runs from 0 to 1, so every
## distance (mesh edges, range, CV radii) is a fraction of the domain length.
coords_stations   <- cbind(stations$x,  stations$y)
coords_prediction <- cbind(pred_grid$x, pred_grid$y)
y_log             <- log10(stations$concentration)

################ STEP 0: INVESTIGATE DATA ###################
#investigate data distribution
p1 <- ggplot(stations, aes(x = concentration)) +
  geom_histogram(bins = 30, fill = "#4C78A8", colour = "white", alpha = 0.9) +
  labs(title = "Original scale", x = "concentration (mg/L)", y = "Count") +
  theme_minimal(base_size = 13)

# maybe log-transform data
p2 <- ggplot(stations, aes(x = y_log)) +
  geom_histogram(bins = 30, fill = "#F58518", colour = "white", alpha = 0.9) +
  labs(title = expression(log[10](concentration)),
       x = expression(log[10](concentration)), y = "Count") +
  theme_minimal(base_size = 13)

p1 + p2 +
  plot_annotation(title = "Distribution of concentration before and after log10 transformation")


################ STEP 1: THE MESH ###################
bnd_segm <- fm_as_segm(boundary)      # boundary as mesh segment (no CRS)

est <- estimate_spatial_range(coords_stations, y_log)
plot_variogram_fit(est, unit = "units")

#build mesh in normalised units
mesh_out <- build_mesh_ensemble(coords_stations, bnd_segm, y_log, n_levels = 6)

mesh_list <- mesh_out$mesh_list

cat("The finest mesh has", mesh_list[[1]]$n, "vertices.")
cat("The coarsest mesh has", mesh_list[[6]]$n, "vertices.")

m1 <- ggplot() + geom_fm(data = mesh_list[[1]]) +
  geom_point(aes(x, y), data = stations)

m6 <- ggplot() + geom_fm(data = mesh_list[[6]]) +
  geom_point(aes(x, y), data = stations)

m1 + m6

my_range <- mesh_out$emp_range        # practical range from the variogram, normalised units

#mesh assessment with inla functions
out1 <- inla.mesh.assessment(mesh_list[[1]], spatial.range = my_range, alpha = 2)
out6 <- inla.mesh.assessment(mesh_list[[6]], spatial.range = my_range, alpha = 2)
summary(out1$sd.dev)
summary(out6$sd.dev)

p1 <- ggplot() + geom_sf(data = sf::st_as_sf(out1), aes(colour = sd.dev)) +
  scale_colour_viridis_c()

p6 <- ggplot() + geom_sf(data = sf::st_as_sf(out6), aes(colour = sd.dev)) +
  scale_colour_viridis_c()

p1+p6

#choose final mesh
mesh_113AK <- mesh_list[[6]]

################# STEP 2: BUILD A MATRICES  #################

# alternative: A_matrix <- inla.spde.make.A(mesh_113AK, loc = coords_stations)
A_matrix <- fm_basis(mesh_113AK, loc = coords_stations)
A_pred   <- fm_basis(mesh_113AK, loc = coords_prediction)

cat("The dimension of A_matrix is", dim(A_matrix))
cat("Rows: nr of measurements", dim(A_matrix)[1])
cat("Columns: nr of mesh points", dim(A_matrix)[2])

plot_A_row(mesh_113AK, A_matrix, coords_stations, i = 10)


################ STEP 3: BUILD SPDE #################

#no priors
spde_no_priors <- inla.spde2.matern(mesh_113AK, alpha = 2)

#highly recommend to use Penalised Complexity (PC) priors
spde <- inla.spde2.pcmatern(
  mesh        = mesh_113AK,
  #estimates from estimate_spatial_range, in normalised units
  prior.range = c(my_range, 0.5),  #P(range < my_range) = 50%
  prior.sigma = c(0.3, 0.5),       #P(sigma > 0.3) = 50%, log10 scale, no unit change
  alpha       = 2)

#index: names the latent field in the formula
#spde$n.spde = mesh_113AK$n
index_space <- inla.spde.make.index("spatial_field", n.spde = spde$n.spde)
str(index_space) # spatial_field = 1..n.spde, plus .group and .repl (used for space-time models)

#define the formula
formula <- y ~ -1 + Intercept +  #-1 removes the default intercept
  altitude + order_min + IDPR + river_dist + fault_dist + dist_to_gypse +  #the covariates
  f(spatial_field, model = spde)  #the spatial field, called via the index name

################ STEP 4: BUILD STACK #################

#observation stack
stk_est <- inla.stack(
  data = list(y = y_log),
  A = list(A_matrix, 1), # one projection matrix per effect block:
  # A_matrix for the field, 1 (identity) for the covariates
  effects = list(index_space,  # block 1: spatial field, indexed by mesh node
                 data.frame(Intercept = 1, stations[, covar_names])),
  tag = "est")  #through this tag you find the observation data in the stack

#for prediction: YOU NEED COVARIATE INFO AT EACH GRID POINT
stk_pred <- inla.stack(
  data = list(y = NA),   #prediction data NA (to be estimated)
  A = list(A_pred, 1),
  effects = list(index_space,
                 data.frame(Intercept = 1, pred_grid[, covar_names])),
  tag = "pred")  #through this tag you find the prediction data in the stack

#join stack
stk <- inla.stack(stk_est, stk_pred)
dim(inla.stack.A(stk))   # stations + grid points rows

## A must have one row per observation
stopifnot(nrow(A_matrix) == nrow(stations),
          nrow(A_pred)   == nrow(pred_grid))

################ STEP 5: INLA CALL and MODEL COMPARISON #################
###### MODEL 1 #####
# 5.1 INLA call
res <- inla(formula, family = "gaussian",
            data = inla.stack.data(stk, spde = spde),
            #control.family = list(link='log'),   #for different distribution, different link functions
            #quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),  #uncomment but heavier computation
            control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE, #cpo, config, marginals feed CV, sampling and INLAutils
                                   #dic and waic information criteria
                                   return.marginals.predictor = FALSE),
            control.predictor = list(A = inla.stack.A(stk), compute = TRUE, link = 1), #link = 1 uses the family link for NA rows
            # control.inla = list(int.strategy = "eb",  # Empirical Bayes (faster)
            #                    strategy = "adaptive"),
            verbose = FALSE)   #TRUE and you can see detailed computations

#summary statistics show no effect of distance covariates
summary(res)

res$cpu.used          # seconds: pre-processing, running, post-processing

###### MODEL 2 #####
#scale covariates
cont_vars <- c("altitude", "IDPR", "river_dist", "fault_dist", "dist_to_gypse")

# centre and sd from the observations only, applied to both data sets
ctr <- sapply(stations[, cont_vars], mean, na.rm = TRUE)
scl <- sapply(stations[, cont_vars], sd,   na.rm = TRUE)

covars_s <- stations
pred_s   <- pred_grid
covars_s[, cont_vars] <- scale(stations[, cont_vars],  center = ctr, scale = scl)
pred_s[, cont_vars]   <- scale(pred_grid[, cont_vars], center = ctr, scale = scl)


stk_est_scaled <- inla.stack(
  data    = list(y = y_log),
  A       = list(A_matrix, 1),       # one projection matrix per effect block:
  # A_matrix for the field, 1 (identity) for the covariates
  effects = list(index_space,        # block 1: spatial field, indexed by mesh node
                 data.frame(Intercept = 1, covars_s[, covar_names])),   # block 2: scaled covariates
  tag     = "est")                   # find these rows with inla.stack.index(stk, "est")$data

stk_pred_scaled <- inla.stack(
  data    = list(y = NA),            # unknown, INLA returns the posterior
  A       = list(A_pred, 1),
  effects = list(index_space,
                 data.frame(Intercept = 1, pred_s[, covar_names])),     # same scaling as the observations
  tag     = "pred")                  # find these rows with inla.stack.index(stk, "pred")$data

stk_scaled <- inla.stack(stk_est_scaled, stk_pred_scaled)

# 5.2 INLA call with scaled distance covariates
res_scaled <- inla(formula, family = "gaussian",
                   data = inla.stack.data(stk_scaled, spde = spde),
                   control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE,
                                          return.marginals.predictor = TRUE),
                   control.predictor = list(A = inla.stack.A(stk_scaled), compute = TRUE, link = 1),
                   verbose = FALSE)

#summary statistics
summary(res_scaled)

res_scaled$cpu.used          # seconds: pre-processing, running, post-processing

###### MODEL 3 #####
# remove distance covariates
formula_no_dist <- y ~ -1 + Intercept + altitude + order_min + IDPR + f(spatial_field, model = spde)

# 5.3 INLA call with no distance covariates
res_scaled_no_dist <- inla(formula_no_dist, family = "gaussian",
                           data = inla.stack.data(stk_scaled, spde = spde),
                           control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE,
                                                  return.marginals.predictor = FALSE),
                           control.predictor = list(A = inla.stack.A(stk_scaled), compute = TRUE, link = 1),
                           verbose = FALSE)

summary(res_scaled_no_dist)

# test which model is the best according to DIC and WAIC
fits <- list(
  original                  = res,
  scaled                    = res_scaled,
  scaled_no_distance_covars = res_scaled_no_dist
)

model_scores <- function(fits) {
  do.call(rbind, lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    data.frame(
      model        = nm,
      DIC          = f$dic$dic,      # the lower the better
      WAIC         = f$waic$waic,    # the lower the better
      LOO_logscore = -mean(log(f$cpo$cpo), na.rm = TRUE), # NA rows (predictions) are skipped, the lower the better
      CPO_failures = sum(f$cpo$failure > 0, na.rm = TRUE)  #0 failures: you can trust the log score
    )
  }))
}

scores <- model_scores(fits)
scores[order(scores$DIC), ]          # best model first


################ STEP 6: MODEL OUTPUT OF BEST PERFORMING ONE #################
# tags let us find rows again after fitting
idx_est  <- inla.stack.index(stk_scaled, tag = "est")$data
idx_pred <- inla.stack.index(stk_scaled, tag = "pred")$data

###### 6a. Fixed effects ######
head(res_scaled$summary.fixed)     #per 1 sd of each covariate

######  6b. Hyperparameters ######
res_scaled$summary.hyperpar                   # noise precision, range (normalised units), field sd

#noise sd instead of precision
noise_sd <- inla.tmarginal(function(x) 1 / sqrt(x),
                           res_scaled$marginals.hyperpar[["Precision for the Gaussian observations"]])
inla.zmarginal(noise_sd)

# nugget share: how much residual variance is measurement noise or is the spatial field really worth it?
s <- inla.hyperpar.sample(4000, res_scaled)                                    # columns in the order of summary.hyperpar
noise_var <- 1 / s[, "Precision for the Gaussian observations"]                # precision -> variance
field_var <- s[, "Stdev for spatial_field"]^2                                  # Stdev for spatial_field -> variance
ratio <- noise_var / (noise_var + field_var)
quantile(ratio, c(0.025, 0.5, 0.975))

######  6c. Fit at the stations  ######
fitted_est <- res_scaled$summary.fitted.values[idx_est, c("mean", "sd")]

## Observed vs fitted at the data locations
plot(y_log, fitted_est$mean, xlab = "observed", ylab = "fitted", asp = 1)
abline(0, 1, col = "red")

#alternative with intervals
ggplot_inla_residuals(res_scaled, observed = y_log, CI = TRUE)


################ STEP 7: PLOTS #################

#shows posterior distributions of effects
#autoplot(res_scaled)   #takes a bit of time

#plot posterior distributions of effects
#which=1 fixed effects, which=2 random effects, which=3 spatial field
autoplot(res_scaled, which = 1, CI = TRUE)

#plot spatial mean to the mesh
proj <- fm_evaluator(mesh_113AK, dims = c(200, 200))
mean_lattice <- fm_evaluate(proj, field = res_scaled$summary.random$spatial_field$mean)

mean_df <- data.frame(
  x    = rep(proj$x, times = length(proj$y)),
  y    = rep(proj$y, each  = length(proj$x)),
  mean = as.vector(mean_lattice)
)

#The field has mean zero
#It shows where the log10-data is above (orange) or below (blue) what covariates predict
#A diverging scale centred at 0 makes that visible.
ggplot() +
  geom_raster(data = mean_df, aes(x, y, fill = mean)) +
  geom_fm(data = mesh_113AK, alpha = 0, linewidth = 0.05, colour = "grey70") +  # sets coord_sf
  geom_point(data = stations, aes(x, y), colour = "black", size = 0.8) +
  scale_fill_gradient2(low = "#2F6F8F", mid = "#F7F5F0", high = "#A54E0E",
                       midpoint = 0, na.value = "transparent",
                       name = "mean of u(s)") +
  labs(x = "x", y = "y",
       title = "Posterior mean of the spatial field") +
  theme_minimal()

pred_grid$pred_mean <- res_scaled$summary.fitted.values$mean[idx_pred]
pred_grid$pred_sd   <- res_scaled$summary.fitted.values$sd[idx_pred]

plot_mean_log10 <- ggplot(pred_grid, aes(x, y, fill = pred_mean)) +
  geom_raster() +
  scale_fill_viridis_c(name = "mean\nlog10") +
  coord_equal()

plot_sd_log10 <- ggplot(pred_grid, aes(x, y, fill = pred_sd)) +
  geom_raster() +
  scale_fill_viridis_c(name = "sd\nlog10") +
  coord_equal()

plot_mean_log10 + plot_sd_log10

#alternative: project SD on mesh
sd_lattice <- fm_evaluate(proj, field = res_scaled$summary.random$spatial_field$sd)

sd_df <- data.frame(
  x  = rep(proj$x, times = length(proj$y)),
  y  = rep(proj$y, each  = length(proj$x)),
  sd = as.vector(sd_lattice)
)

ggplot() +
  geom_raster(data = sd_df, aes(x, y, fill = sd)) +
  geom_fm(data = mesh_113AK, alpha = 0, linewidth = 0.05, colour = "grey70") +  # sets coord_sf
  geom_point(data = stations, aes(x, y), colour = "white", size = 0.8) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent",
                       name = "sd of u(s)") +
  labs(x = "x (normalised)", y = "y (normalised)",
       title = "Posterior sd of the spatial field") +
  theme_minimal()


#transform data back from log10-scale
s <- res_scaled$summary.fitted.values[idx_pred, ]

#10^x is monotone, so quantiles carry over exactly:
pred_grid$conc_median <- 10^s$`0.5quant`     # posterior median, original units
pred_grid$conc_lo     <- 10^s$`0.025quant`   # 95% interval, lower
pred_grid$conc_hi     <- 10^s$`0.975quant`   # 95% interval, upper

ggplot(pred_grid, aes(x, y, fill = conc_median)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(trans = "log10", name = "median")

#for the mean and sd, the posterior of η is close to Gaussian, so the lognormal formulas apply:
k  <- log(10)
mu <- s$mean
sg <- s$sd

pred_grid$conc_mean <- 10^mu * exp((k * sg)^2 / 2)
pred_grid$conc_sd   <- pred_grid$conc_mean * sqrt(exp((k * sg)^2) - 1)

plot_mean <- ggplot(pred_grid, aes(x, y, fill = conc_mean)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(trans = "log10", name = "posterior mean")

plot_sd <- ggplot(pred_grid, aes(x, y, fill = conc_sd)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(trans = "log10", name = "posterior sd")

plot_mean + plot_sd


################ STEP 8: CROSS-VALIDATION ################

#8a. LOO and LGOCV (no refit, all stations)
# Leave-one-out checks
cpo <- res_scaled$cpo$cpo[idx_est]        # check the CPO at observed stations

i <- which.min(cpo)    # the station the model predicts worst
cat("The worst predicted station is", i)
cat("It is located at (normalised) coords", coords_stations[i, ])
cat("The quantity on that station is", stations$concentration[i], "mg/L.")

hist(res_scaled$cpo$pit[idx_est])      # PIT: roughly flat (uniform dist) = well calibrated

#alternative: use inla.group.cv with num.level.sets = -1 for LOO
loo <- inla.group.cv(res_scaled, num.level.sets = -1)  # -1 = LOO

#leave group out: for each station i, INLA builds a group I_i of rows that are strongly
#correlated with eta_i in the posterior. The groups do not depend on distance.
#It then removes y_i together with the whole group and predicts y_i from the rest.
lgocv <- inla.group.cv(res_scaled, num.level.sets = 3)    # takes about 2 min

-mean(log(loo$cv[idx_est]))          # LOO log score
-mean(log(lgocv$cv[idx_est]))        # LGOCV log score
lgocv$groups[[1]]$idx                # rows left out together with station 1


# 8b. Disc CV: group i = all stations within rad of station i (no refit)
# radii as multiples of the practical range from the variogram (normalised units)
rads <- c(0, 0.25, 0.5, 1, 1.5, 2) * my_range

# Row = station number, because the est rows come first in stk_scaled
D <- as.matrix(dist(coords_stations))       # station-to-station distances (normalised units)

disc_groups <- function(rad) lapply(seq_along(idx_est), function(i) idx_est[D[i, ] <= rad])
disc_cv <- lapply(rads, function(r) inla.group.cv(res_scaled, groups = disc_groups(r)))

k <- length(rads)                              # any radius from rads, here the largest
i <- which.min(res_scaled$cpo$cpo[idx_est])    # any station, here the worst predicted

plot_cv_group(i, disc_cv[[k]]$groups[[i]]$idx, rad = rads[k],
              title = sprintf("Disc CV, rad = %.2f (%.1f x range)", rads[k], rads[k] / my_range))


# Columns of cv_all:
# method     : LOO         = leave one station out
#              LGOCV m = 3 = leave out a group of correlated rows
#              Disc CV     = leave out all stations within rad
#              All three use res_scaled without refit: priors and hyperparameters
#              stay as in the full model.
# rad_m      : radius of the removed disc in normalised units. NA = not distance based.
# n          : number of stations scored (all stations).
# group_size : mean number of stations left out per test station, including
#              the test station itself. 1 = ordinary LOO.
# logscore   : -mean(log predictive density at the observed value). Lower is better.
#              Rewards an accurate mean and a correct spread. Can be negative.
# RMSE       : root mean squared error of the predictive mean, log10 units. Lower is better.
# MAE        : mean absolute error of the predictive mean, log10 units.
#              Less sensitive to single large errors than RMSE.
# CRPS       : continuous ranked probability score, log10 units. Lower is better.
#              Compares the whole predictive distribution with the observed value.
#              Equals the MAE when the predictive sd shrinks to 0.
# cover95    : share of stations whose observed value falls inside the 95 %
#              predictive interval.
# The predictive sd is the sd of the linear predictor only (no measurement noise).
# The scores refer to the mean surface, so cover95 falls below 0.95 even for a
# correct model.
cv_all <- rbind(
  as_row(cv_scores(loo),   "LOO"),
  as_row(cv_scores(lgocv), "LGOCV m = 3"),
  do.call(rbind, Map(function(g, r) as_row(cv_scores(g), "Disc CV", r), disc_cv, rads)))
rownames(cv_all) <- NULL

cv_all


