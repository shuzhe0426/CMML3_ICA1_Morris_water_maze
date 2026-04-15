in_path <- if (file.exists("plot_data_runs.csv")) {
  "plot_data_runs.csv"
} else if (file.exists("plot_data_runs_small.csv")) {
  "plot_data_runs_small.csv"
} else {
  stop("plot_data_runs.csv not found; run export_plot_data.R first")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
library(ggplot2)

d <- read.csv(in_path)

ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(mean = mean(x), lo = NA_real_, hi = NA_real_))
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  t <- qt(0.975, df = n - 1)
  c(mean = m, lo = m - t * se, hi = m + t * se)
}

extract_ci_from_agg <- function(agg, value_col) {
  v <- agg[[value_col]]
  if (is.matrix(v) && all(c("mean", "lo", "hi") %in% colnames(v))) {
    list(mean = v[, "mean"], lo = v[, "lo"], hi = v[, "hi"])
  } else if (is.data.frame(v) && all(c("mean", "lo", "hi") %in% names(v))) {
    list(mean = v$mean, lo = v$lo, hi = v$hi)
  } else if (is.list(v)) {
    m <- do.call(rbind, v)
    list(mean = m[, "mean"], lo = m[, "lo"], hi = m[, "hi"])
  } else {
    stop("Unexpected aggregate column structure")
  }
}

plot_day_ci <- function(metric_col, y_label, out_pdf, include_combined_weight) {
  sub <- d[d$model != "combined" | (is.finite(d$weight_wall) & abs(d$weight_wall - include_combined_weight) < 1e-9), ]
  sub$model_label <- ifelse(sub$model == "combined", paste0("combined (w=", include_combined_weight, ")"), sub$model)
  sub$model_label <- factor(sub$model_label, levels = c("place", "distance", paste0("combined (w=", include_combined_weight, ")")))
  sub$task <- factor(sub$task, levels = c("fixed", "variable"))

  agg <- aggregate(sub[[metric_col]] ~ task + model_label + day, sub, function(x) ci95(x))
  value_col <- setdiff(names(agg), c("task", "model_label", "day"))[1]
  ci <- extract_ci_from_agg(agg, value_col)
  agg$mean <- ci$mean
  agg$lo <- ci$lo
  agg$hi <- ci$hi

  p <- ggplot(agg, aes(x = day, y = mean, color = model_label, fill = model_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, linewidth = 0) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ task, nrow = 1) +
    labs(x = "Day", y = y_label, color = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(out_pdf, p, width = 9, height = 4.2)
}

plot_day_ci("latency", "Latency (s)", "fig1_latency_vs_day_ci.pdf", 0.2)
plot_day_ci("wall_zone", "Wall zone (%)", "figS1_wall_zone_vs_day_ci.pdf", 0.2)

dc <- d[d$model == "combined" & is.finite(d$weight_wall), ]
dc$task <- factor(dc$task, levels = c("fixed", "variable"))

per_run <- aggregate(cbind(latency, wall_zone, target_quadrant) ~ task + weight_wall + run, dc, mean)

ratio_ci <- function(x) ci95(x)

agg_latency <- aggregate(latency ~ task + weight_wall, per_run, ratio_ci)
lat_value_col <- setdiff(names(agg_latency), c("task", "weight_wall"))[1]
lat_ci <- extract_ci_from_agg(agg_latency, lat_value_col)
agg_latency$mean <- lat_ci$mean
agg_latency$lo <- lat_ci$lo
agg_latency$hi <- lat_ci$hi
agg_latency$metric <- "Latency (s)"
agg_latency$value <- agg_latency$mean

agg_wall <- aggregate(wall_zone ~ task + weight_wall, per_run, ratio_ci)
wall_value_col <- setdiff(names(agg_wall), c("task", "weight_wall"))[1]
wall_ci <- extract_ci_from_agg(agg_wall, wall_value_col)
agg_wall$mean <- wall_ci$mean
agg_wall$lo <- wall_ci$lo
agg_wall$hi <- wall_ci$hi
agg_wall$metric <- "Wall zone (%)"
agg_wall$value <- agg_wall$mean

agg_target <- aggregate(target_quadrant ~ task + weight_wall, per_run, ratio_ci)
target_value_col <- setdiff(names(agg_target), c("task", "weight_wall"))[1]
target_ci <- extract_ci_from_agg(agg_target, target_value_col)
agg_target$mean <- target_ci$mean
agg_target$lo <- target_ci$lo
agg_target$hi <- target_ci$hi
agg_target$metric <- "Target quadrant (%)"
agg_target$value <- agg_target$mean

ratio <- rbind(
  agg_latency[, c("task", "weight_wall", "metric", "value", "lo", "hi")],
  agg_wall[, c("task", "weight_wall", "metric", "value", "lo", "hi")],
  agg_target[, c("task", "weight_wall", "metric", "value", "lo", "hi")]
)

ratio$metric <- factor(ratio$metric, levels = c("Latency (s)", "Wall zone (%)", "Target quadrant (%)"))

p_ratio <- ggplot(ratio, aes(x = weight_wall, y = value)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03, linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_grid(metric ~ task, scales = "free_y") +
  labs(x = "weight_wall", y = NULL) +
  theme_bw(base_size = 12)

ggsave("fig2_combined_ratio_tradeoff_ci.pdf", p_ratio, width = 9, height = 7.2)
