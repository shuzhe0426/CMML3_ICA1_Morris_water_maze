set.seed(1)

pool_diameter <- 1.4
platform_radius <- 0.06

platform_x <- cos(-pi/4) * pool_diameter/4
platform_y <- sin(-pi/4) * pool_diameter/4

strad <- pool_diameter/2 * 0.85
starting_xs <- strad * c(cos(pi/6), cos(pi/3), cos(7*pi/6), cos(4*pi/3))
starting_ys <- strad * c(sin(pi/6), sin(pi/3), sin(7*pi/6), sin(4*pi/3))
starting_x <- starting_xs[1]
starting_y <- starting_ys[1]

th <- (0:200)/100*pi

load_run_trial <- function(file_path, stop_at_pattern) {
  code <- readLines(file_path, warn = FALSE)
  stop_idx <- match(TRUE, grepl(stop_at_pattern, code))
  if (is.na(stop_idx)) stop(paste0("Cannot find stop pattern in ", file_path, ": ", stop_at_pattern))
  code <- code[1:(stop_idx - 1)]
  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  if (!exists("run_trial", envir = env)) stop(paste0("run_trial not defined in truncated: ", file_path))
  get("run_trial", envir = env)
}

run_trial_place <- load_run_trial("place-cell model.R", "^N_pc <-")
run_trial_distance <- load_run_trial("distance-cell model.R", "^N_dc <-")
run_trial_combined <- load_run_trial("place distance cells combined model.R", "^# =========================== main")

N_ac <- 36

place_params <- list(
  N_pc = 211,
  sigma_pc = 0.1,
  sigma_ac = 2,
  etdecay = 0.83,
  beta = 6,
  alpha = 0.01,
  gamma = 0.85,
  Vdecay = 0.82,
  ac_const = 0.02,
  Wnoise = 0.0004,
  Wmult = 0.1,
  hitwall = 0.5,
  speed = 0.175
)

distance_params <- list(
  N_dc = 10,
  sigma_dc = 0.1,
  sigma_ac = 2,
  etdecay = 0.83,
  beta = 6,
  alpha = 0.01,
  gamma = 0.85,
  Vdecay = 0.82,
  ac_const = 0.02,
  Wnoise = 0.0004,
  Wmult = 0.1,
  hitwall = 0.5,
  speed = 0.175
)

combined_params <- list(
  N_pc = 211,
  N_dc = 10,
  sigma_pc = 0.1,
  sigma_dc = 0.1,
  sigma_ac = 2,
  etdecay = 0.83,
  beta = 6,
  alpha = 0.006,
  gamma = 0.85,
  Vdecay = 0.82,
  ac_const = 0.02,
  Wnoise = 0.0004,
  Wmult = 0.1,
  hitwall = 0,
  speed = 0.175,
  weight_wall = 0.2
)

make_place_cells <- function(N_pc) {
  PC_x <- rep(0, N_pc)
  PC_y <- rep(0, N_pc)
  for (i in 1:N_pc) {
    PC_x[i] <- (stats::runif(1) - 0.5) * pool_diameter
    PC_y[i] <- (stats::runif(1) - 0.5) * pool_diameter
    while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)) {
      PC_x[i] <- (stats::runif(1) - 0.5) * pool_diameter
      PC_y[i] <- (stats::runif(1) - 0.5) * pool_diameter
    }
  }
  list(PC_x = PC_x, PC_y = PC_y)
}

make_distance_cells <- function(N_dc) {
  DC <- rep(0, N_dc)
  for (i in 1:N_dc) DC[i] <- stats::runif(1) * (pool_diameter/2)
  DC
}

pcs <- make_place_cells(place_params$N_pc)
dcs <- make_distance_cells(distance_params$N_dc)

weights_place <- matrix(stats::runif(place_params$N_pc * N_ac), nrow = place_params$N_pc) * place_params$Wmult
weights_distance <- matrix(stats::runif(distance_params$N_dc * N_ac), nrow = distance_params$N_dc) * distance_params$Wmult

weights_pc <- matrix(stats::runif(combined_params$N_pc * N_ac), nrow = combined_params$N_pc) * combined_params$Wmult
weights_dc <- matrix(stats::runif(combined_params$N_dc * N_ac), nrow = combined_params$N_dc) * combined_params$Wmult

res_place <- run_trial_place(
  weights_place,
  place_params$Wmult,
  place_params$sigma_pc,
  place_params$sigma_ac,
  pcs$PC_x,
  pcs$PC_y,
  place_params$Vdecay,
  place_params$ac_const,
  place_params$beta,
  place_params$etdecay,
  place_params$alpha,
  place_params$gamma,
  place_params$Wnoise,
  platform_x,
  platform_y,
  starting_x,
  starting_y,
  place_params$speed,
  place_params$hitwall
)

res_distance <- run_trial_distance(
  weights_distance,
  distance_params$Wmult,
  distance_params$sigma_dc,
  distance_params$sigma_ac,
  dcs,
  distance_params$Vdecay,
  distance_params$ac_const,
  distance_params$beta,
  distance_params$etdecay,
  distance_params$alpha,
  distance_params$gamma,
  distance_params$Wnoise,
  platform_x,
  platform_y,
  starting_x,
  starting_y,
  distance_params$speed,
  distance_params$hitwall
)

res_combined <- run_trial_combined(
  weights_pc,
  weights_dc,
  combined_params$Wmult,
  combined_params$sigma_pc,
  combined_params$sigma_dc,
  combined_params$sigma_ac,
  pcs$PC_x,
  pcs$PC_y,
  dcs,
  combined_params$Vdecay,
  combined_params$ac_const,
  combined_params$beta,
  combined_params$etdecay,
  combined_params$alpha,
  combined_params$gamma,
  combined_params$Wnoise,
  platform_x,
  platform_y,
  starting_x,
  starting_y,
  combined_params$speed,
  combined_params$hitwall,
  combined_params$weight_wall,
  1 - combined_params$weight_wall
)

pdf("figS2_example_trajectories_fixed.pdf", width = 10.5, height = 3.8)
par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))

plot_one <- function(track_x, track_y, title) {
  plot(pool_diameter/2*cos(th), pool_diameter/2*sin(th), type = "l", asp = 1, xlab = "x", ylab = "y", main = title)
  lines(track_x, track_y, col = "black")
  lines(platform_x + platform_radius*cos(th), platform_y + platform_radius*sin(th), col = "red")
}

plot_one(res_place[[2]], res_place[[3]], "Place (fixed)")
plot_one(res_distance[[2]], res_distance[[3]], "Distance (fixed)")
plot_one(res_combined[[3]], res_combined[[4]], "Combined w=0.2 (fixed)")

dev.off()
