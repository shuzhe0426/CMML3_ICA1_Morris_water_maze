set.seed(1)

source_path <- "distance-cell model.R"
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

disable_write_csv <- function(code) {
  idx <- grep("write\\.csv\\(", code)
  if (length(idx) > 0) code[idx] <- "invisible(NULL)"
  code
}

insert_after_first <- function(code, needle, new_lines) {
  idx <- match(TRUE, grepl(needle, code, fixed = TRUE))
  if (is.na(idx)) stop(paste0("Cannot find insertion point: ", needle))
  c(code[1:idx], new_lines, code[(idx + 1):length(code)])
}

code <- update_scalar(code, "plot_trajectories", 0)
code <- update_scalar(code, "plot_cognitive_maps", 0)
code <- update_scalar(code, "Nruns", 50)
code <- update_scalar(code, "variable_platform", 1)
code <- update_scalar(code, "Ntrials", 4)
code <- update_scalar(code, "Ndays", 8)
code <- disable_write_csv(code)

code <- insert_after_first(
  code,
  "platform_y <-",
  c("platform_x0 <- platform_x", "platform_y0 <- platform_y")
)
code <- gsub("platform_x = platform_x\\*whichplatform", "platform_x = platform_x0*whichplatform", code, fixed = TRUE)
code <- gsub("platform_y = platform_y\\*whichplatform", "platform_y = platform_y0*whichplatform", code, fixed = TRUE)

env <- new.env(parent = globalenv())
eval(parse(text = code), envir = env)

if (!exists("PMs", envir = env)) stop("PMs not found after running script")
batch <- summarize_pms(get("PMs", envir = env))
write.csv(batch, "distance_variable_batch.csv", row.names = FALSE)
