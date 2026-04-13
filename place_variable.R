set.seed(1)

source_path <- "place-cell model.R"
code <- readLines(source_path, warn = FALSE)

update_scalar <- function(code, name, value) {
  idx <- grep(paste0("^", name, " <- "), code)
  if (length(idx) != 1) stop(paste0("Expected exactly one assignment for: ", name))
  code[idx] <- sub(paste0("^", name, " <- [^#]*"), paste0(name, " <- ", value, " "), code[idx])
  code
}

summarize_pms <- function(PMs) {
  d <- dim(PMs)
  if (length(d) == 3) PMs <- array(PMs, dim = c(d[1], d[2], d[3], 1))
  d <- dim(PMs)
  if (length(d) == 5) PMs <- PMs[, , , , 1, drop = FALSE]
  d <- dim(PMs)
  if (length(d) != 4) stop("Unsupported PMs dimensions")
  Ndays <- d[2]
  out <- data.frame(
    day = seq_len(Ndays),
    latency = apply(PMs[1, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    dist = apply(PMs[2, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    target_quadrant = apply(PMs[3, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    opposite_quadrant = apply(PMs[4, , , , drop = FALSE], 2, mean, na.rm = TRUE),
    wall_zone = apply(PMs[5, , , , drop = FALSE], 2, mean, na.rm = TRUE)
  )
  out
}

insert_after_first <- function(code, needle, new_lines) {
  idx <- match(TRUE, grepl(needle, code, fixed = TRUE))
  if (is.na(idx)) stop(paste0("Cannot find insertion point: ", needle))
  c(code[1:idx], new_lines, code[(idx + 1):length(code)])
}

insert_after_all <- function(code, needle, new_lines) {
  idxs <- which(grepl(needle, code, fixed = TRUE))
  if (length(idxs) < 1) stop(paste0("Cannot find insertion point: ", needle))
  for (idx in rev(idxs)) {
    code <- c(code[1:idx], new_lines, code[(idx + 1):length(code)])
  }
  code
}

code <- update_scalar(code, "plot_trajectories", 0)
code <- update_scalar(code, "plot_cognitive_maps", 0)
code <- update_scalar(code, "Nruns", 50)
code <- update_scalar(code, "Ntrials", 4)
code <- update_scalar(code, "Ndays", 8)

code <- insert_after_first(
  code,
  "platform_y <-",
  c("platform_x0 <- platform_x", "platform_y0 <- platform_y", "variable_platform <- 1")
)

code <- insert_after_all(
  code,
  "for (trial in 1:Ntrials){",
  c(
    "    whichplatform=sample(c(1,-1),1)",
    "    px = platform_x0*whichplatform",
    "    py = platform_y0*whichplatform"
  )
)

mod_idxs <- grep("modresults <- run_trial", code, fixed = TRUE)
if (length(mod_idxs) < 1) stop("Cannot find run_trial call to patch")
code[mod_idxs] <- gsub("platform_x, platform_y", "px, py", code[mod_idxs], fixed = TRUE)

env <- new.env(parent = globalenv())
eval(parse(text = code), envir = env)

if (!exists("PMs", envir = env)) stop("PMs not found after running script")
batch <- summarize_pms(get("PMs", envir = env))
write.csv(batch, "place_variable_batch.csv", row.names = FALSE)
