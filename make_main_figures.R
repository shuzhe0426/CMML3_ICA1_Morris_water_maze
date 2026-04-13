files <- data.frame(
  group = c(
    "place_fixed",
    "place_variable",
    "distance_fixed",
    "distance_variable",
    "combined_fixed_w02",
    "combined_variable_w02"
  ),
  path = c(
    "place_fixed_batch.csv",
    "place_variable_batch.csv",
    "distance_fixed_batch.csv",
    "distance_variable_batch.csv",
    "combined_fixed_w02_batch.csv",
    "combined_variable_w02_batch.csv"
  ),
  stringsAsFactors = FALSE
)

read_one <- function(group, path) {
  d <- read.csv(path)
  d$group <- group
  d
}

all <- do.call(rbind, Map(read_one, files$group, files$path))

all$model <- ifelse(grepl("^place_", all$group), "place",
             ifelse(grepl("^distance_", all$group), "distance", "combined"))
all$cond <- ifelse(grepl("_fixed$", all$group), "fixed", "variable")

all$model_label <- ifelse(all$model == "combined", "combined (w=0.2)", all$model)
all$model_label <- factor(all$model_label, levels = c("place", "distance", "combined (w=0.2)"))
all$cond <- factor(all$cond, levels = c("fixed", "variable"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
library(ggplot2)

p_latency <- ggplot(all, aes(x = day, y = latency, color = model_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ cond, nrow = 1) +
  labs(x = "Day", y = "Latency (s)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("fig1_latency_vs_day.pdf", p_latency, width = 9, height = 4.2)

p_wall <- ggplot(all, aes(x = day, y = wall_zone, color = model_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ cond, nrow = 1) +
  labs(x = "Day", y = "Wall zone (%)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("fig2_wall_zone_vs_day.pdf", p_wall, width = 9, height = 4.2)

ratio_files <- data.frame(
  weight_wall = c(0.2, 0.5, 0.8, 0.2, 0.5, 0.8),
  cond = c("fixed", "fixed", "fixed", "variable", "variable", "variable"),
  path = c(
    "combined_fixed_w02_batch.csv",
    "combined_fixed_w05_batch.csv",
    "combined_fixed_w08_batch.csv",
    "combined_variable_w02_batch.csv",
    "combined_variable_w05_batch.csv",
    "combined_variable_w08_batch.csv"
  ),
  stringsAsFactors = FALSE
)

if (all(file.exists(ratio_files$path))) {
  read_ratio <- function(weight_wall, cond, path) {
    d <- read.csv(path)
    data.frame(weight_wall = weight_wall, cond = cond, d)
  }
  ratio_all <- do.call(rbind, Map(read_ratio, ratio_files$weight_wall, ratio_files$cond, ratio_files$path))
  ratio_mean <- aggregate(cbind(latency, wall_zone, target_quadrant) ~ weight_wall + cond, ratio_all, mean)

  ratio_long <- rbind(
    data.frame(weight_wall = ratio_mean$weight_wall, cond = ratio_mean$cond, metric = "Latency (s)", value = ratio_mean$latency),
    data.frame(weight_wall = ratio_mean$weight_wall, cond = ratio_mean$cond, metric = "Wall zone (%)", value = ratio_mean$wall_zone),
    data.frame(weight_wall = ratio_mean$weight_wall, cond = ratio_mean$cond, metric = "Target quadrant (%)", value = ratio_mean$target_quadrant)
  )
  ratio_long$cond <- factor(ratio_long$cond, levels = c("fixed", "variable"))
  ratio_long$metric <- factor(ratio_long$metric, levels = c("Latency (s)", "Wall zone (%)", "Target quadrant (%)"))

  p_ratio <- ggplot(ratio_long, aes(x = weight_wall, y = value, group = cond)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(metric ~ cond, scales = "free_y") +
    labs(x = "weight_wall", y = NULL) +
    theme_bw(base_size = 12)

  ggsave("fig3_combined_ratio_tradeoff.pdf", p_ratio, width = 9, height = 7.2)
}
