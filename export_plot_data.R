set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[1], fixed = TRUE)
}

Nruns <- as.integer(get_arg("nruns", "50"))
Ndays <- as.integer(get_arg("ndays", "8"))
Ntrials <- as.integer(get_arg("ntrials", "4"))
out_path <- get_arg("out", "plot_data_runs.csv")
combined_weights <- strsplit(get_arg("weights", "0.2,0.5,0.8"), ",", fixed = TRUE)[[1]]
combined_weights <- as.numeric(combined_weights)

update_scalar <- function(code, name, value) {
  idx <- grep(paste0("^", name, " <- "), code)
  if (length(idx) != 1) stop(paste0("Expected exactly one assignment for: ", name))
  code[idx] <- sub(paste0("^", name, " <- [^#]*"), paste0(name, " <- ", value, " "), code[idx])
  code
}

update_equals <- function(code, name, value) {
  idx <- grep(paste0("^", name, " = "), code)
  if (length(idx) != 1) stop(paste0("Expected exactly one assignment for: ", name))
  code[idx] <- sub(paste0("^", name, " = [^#]*"), paste0(name, " = ", value, " "), code[idx])
  code
}

disable_write_csv <- function(code) {
  idx <- grep("write\\.csv\\(", code)
  if (length(idx) > 0) code[idx] <- "invisible(NULL)"
  code
}

insert_after_all <- function(code, needle, new_lines) {
  idxs <- which(grepl(needle, code, fixed = TRUE))
  if (length(idxs) < 1) stop(paste0("Cannot find insertion point: ", needle))
  for (idx in rev(idxs)) {
    code <- c(code[1:idx], new_lines, code[(idx + 1):length(code)])
  }
  code
}

make_action_selection_safe_combined <- function(code) {
  code <- insert_after_all(
    code,
    "ACsel_pc <- ACsel_pc / sum(ACsel_pc)",
    c("    if (any(!is.finite(ACsel_pc))) { ACsel_pc <- rep(1/N_ac, N_ac) }")
  )
  code <- insert_after_all(
    code,
    "ACsel_dc <- ACsel_dc / sum(ACsel_dc)",
    c("    if (any(!is.finite(ACsel_dc))) { ACsel_dc <- rep(1/N_ac, N_ac) }")
  )
  code <- insert_after_all(
    code,
    "ACsel <- ACsel / sum(ACsel)",
    c("    if (any(!is.finite(ACsel))) { ACsel <- rep(1/N_ac, N_ac) }")
  )
  code
}

extract_run_level <- function(PMs, metric_indices) {
  d <- dim(PMs)
  if (length(d) == 3) PMs <- array(PMs, dim = c(d[1], d[2], d[3], 1))
  d <- dim(PMs)
  if (length(d) == 5) PMs <- PMs[, , , , 1, drop = TRUE]
  d <- dim(PMs)
  if (length(d) != 4) stop("Unsupported PMs dimensions")

  Ndays_local <- d[2]
  Nruns_local <- d[4]
  out <- data.frame()

  for (r in seq_len(Nruns_local)) {
    for (day in seq_len(Ndays_local)) {
      vals <- vapply(metric_indices, function(mi) mean(PMs[mi, day, , r], na.rm = TRUE), numeric(1))
      row <- data.frame(day = day, run = r)
      for (k in names(metric_indices)) row[[k]] <- vals[k]
      out <- rbind(out, row)
    }
  }
  out
}

run_place <- function(task) {
  code <- readLines("place-cell model.R", warn = FALSE)
  code <- update_scalar(code, "plot_trajectories", 0)
  code <- update_scalar(code, "plot_cognitive_maps", 0)
  code <- update_scalar(code, "Nruns", Nruns)
  code <- update_scalar(code, "Ntrials", Ntrials)
  code <- update_scalar(code, "Ndays", Ndays)

  if (task == "variable") {
    code <- insert_after_all(
      code,
      "platform_y <-",
      c("platform_x0 <- platform_x", "platform_y0 <- platform_y")
    )
    code <- insert_after_all(
      code,
      "for (trial in 1:Ntrials){",
      c(
        "      whichplatform=sample(c(1,-1),1)",
        "      px = platform_x0*whichplatform",
        "      py = platform_y0*whichplatform"
      )
    )
    mod_idxs <- grep("modresults <- run_trial", code, fixed = TRUE)
    if (length(mod_idxs) < 1) stop("Cannot find run_trial call to patch")
    code[mod_idxs] <- gsub("platform_x, platform_y", "px, py", code[mod_idxs], fixed = TRUE)
  }

  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  if (!exists("PMs", envir = env)) stop("PMs not found after running place script")

  PMs <- get("PMs", envir = env)
  extract_run_level(PMs, c(latency = 1, target_quadrant = 3, wall_zone = 5))
}

run_distance <- function(task) {
  code <- readLines("distance-cell model.R", warn = FALSE)
  code <- update_scalar(code, "plot_trajectories", 0)
  code <- update_scalar(code, "plot_cognitive_maps", 0)
  code <- update_scalar(code, "Nruns", Nruns)
  code <- update_scalar(code, "Ntrials", Ntrials)
  code <- update_scalar(code, "Ndays", Ndays)
  code <- update_scalar(code, "variable_platform", if (task == "variable") 1 else 0)
  code <- disable_write_csv(code)

  if (task == "variable") {
    code <- insert_after_all(
      code,
      "platform_y <-",
      c("platform_x0 <- platform_x", "platform_y0 <- platform_y")
    )
    code <- gsub("platform_x = platform_x*whichplatform", "platform_x = platform_x0*whichplatform", code, fixed = TRUE)
    code <- gsub("platform_y = platform_y*whichplatform", "platform_y = platform_y0*whichplatform", code, fixed = TRUE)
  }

  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  if (!exists("PMs", envir = env)) stop("PMs not found after running distance script")

  PMs <- get("PMs", envir = env)
  extract_run_level(PMs, c(latency = 1, target_quadrant = 3, wall_zone = 5))
}

run_combined <- function(task, weight_wall) {
  code <- readLines("place distance cells combined model.R", warn = FALSE)

  code <- update_scalar(code, "plot_trajectories", 0)
  code <- update_scalar(code, "plot_cognitive_maps", 0)
  code <- update_scalar(code, "plot_integrated_cognitive_map", 0)
  code <- update_scalar(code, "Nruns", Nruns)
  code <- update_scalar(code, "Ntrials", Ntrials)
  code <- update_scalar(code, "Ndays", Ndays)

  code <- update_equals(code, "weight_wall", weight_wall)
  code <- update_equals(code, "weight_place", "1 - weight_wall")
  code <- disable_write_csv(code)
  code <- make_action_selection_safe_combined(code)

  if (task == "fixed") {
    idx <- grep("^\\s*whichplatform=sample\\(c\\(1,-1\\),1\\)", code)
    if (length(idx) < 1) stop("Cannot find whichplatform sampling line to disable")
    code[idx] <- sub("whichplatform=sample\\(c\\(1,-1\\),1\\)", "whichplatform = 1", code[idx])
  }

  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  if (!exists("PMs", envir = env)) stop("PMs not found after running combined script")

  PMs <- get("PMs", envir = env)
  extract_run_level(PMs, c(latency = 1, target_quadrant = 3, wall_zone = 5))
}

all_rows <- list()

df <- run_place("fixed")
df$task <- "fixed"
df$model <- "place"
df$weight_wall <- NA_real_
all_rows[[length(all_rows) + 1]] <- df

df <- run_place("variable")
df$task <- "variable"
df$model <- "place"
df$weight_wall <- NA_real_
all_rows[[length(all_rows) + 1]] <- df

df <- run_distance("fixed")
df$task <- "fixed"
df$model <- "distance"
df$weight_wall <- NA_real_
all_rows[[length(all_rows) + 1]] <- df

df <- run_distance("variable")
df$task <- "variable"
df$model <- "distance"
df$weight_wall <- NA_real_
all_rows[[length(all_rows) + 1]] <- df

for (w in combined_weights) {
  df <- run_combined("fixed", w)
  df$task <- "fixed"
  df$model <- "combined"
  df$weight_wall <- w
  all_rows[[length(all_rows) + 1]] <- df

  df <- run_combined("variable", w)
  df$task <- "variable"
  df$model <- "combined"
  df$weight_wall <- w
  all_rows[[length(all_rows) + 1]] <- df
}

all <- do.call(rbind, all_rows)
all <- all[, c("task", "model", "weight_wall", "day", "run", "latency", "wall_zone", "target_quadrant")]

write.csv(all, out_path, row.names = FALSE)
