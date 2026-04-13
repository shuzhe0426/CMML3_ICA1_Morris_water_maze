set.seed(1)

source_path <- "place distance cells combined model.R"
code <- readLines(source_path, warn = FALSE)

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

summarize_pms <- function(PMs) {
  d <- dim(PMs)
  if (length(d) == 3) PMs <- array(PMs, dim = c(d[1], d[2], d[3], 1))
  d <- dim(PMs)
  if (length(d) == 5) PMs <- PMs[, , , , 1, drop = TRUE]
  d <- dim(PMs)
  if (length(d) != 4) stop("Unsupported PMs dimensions")
  Ndays <- d[2]
  data.frame(
    day = seq_len(Ndays),
    latency = apply(PMs[1, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    dist = apply(PMs[2, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    target_quadrant = apply(PMs[3, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    opposite_quadrant = apply(PMs[4, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    wall_zone = apply(PMs[5, , , , drop = FALSE], 2, mean, na.rm = TRUE)
  )
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

make_action_selection_safe <- function(code) {
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

disable_variable_platform <- function(code) {
  idx <- grep("^\\s*whichplatform=sample\\(c\\(1,-1\\),1\\)", code)
  if (length(idx) < 1) stop("Cannot find whichplatform sampling line to disable")
  code[idx] <- sub("whichplatform=sample\\(c\\(1,-1\\),1\\)", "whichplatform = 1", code[idx])
  code
}

code <- update_scalar(code, "plot_trajectories", 0)
code <- update_scalar(code, "plot_cognitive_maps", 0)
code <- update_scalar(code, "plot_integrated_cognitive_map", 0)
code <- update_scalar(code, "Nruns", 50)
code <- update_scalar(code, "Ntrials", 4)
code <- update_scalar(code, "Ndays", 8)

code <- update_equals(code, "weight_wall", 0.8)
code <- update_equals(code, "weight_place", "1 - weight_wall")

code <- disable_write_csv(code)
code <- make_action_selection_safe(code)
code <- disable_variable_platform(code)

env <- new.env(parent = globalenv())
eval(parse(text = code), envir = env)

if (!exists("PMs", envir = env)) stop("PMs not found after running script")
batch <- summarize_pms(get("PMs", envir = env))
write.csv(batch, "combined_fixed_w08_batch.csv", row.names = FALSE)
