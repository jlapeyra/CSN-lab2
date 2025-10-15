#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(scales)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_cdf.R <input_file> <output_file> [title]")

input_file  <- args[1]
output_file <- args[2]
plot_title  <- ifelse(length(args) >= 3, args[3], "CDF Comparison (log-log scale)")

# Read data
data <- read.table(input_file, header = TRUE)

# Convert to long format
data_long <- pivot_longer(data, cols = c(EmpiricalCDF, TheoreticalCDF),
                          names_to = "Type", values_to = "CDF")

# Line type: solid for empirical, dotted for theoretical
data_long$LineType <- ifelse(data_long$Type == "EmpiricalCDF", "solid", "dotted")

# Colors: blue for empirical, gray for theoretical
line_colors <- c("EmpiricalCDF" = "blue", "TheoreticalCDF" = "gray40")

# Plot
p <- ggplot(data_long, aes(x = K, y = CDF, color = Type, linetype = LineType)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = Type), size = 1.5) +
  scale_linetype_identity() +
  scale_color_manual(values = line_colors) +
  labs(title = plot_title, x = "Degree (k)", y = "CDF") +
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

# Save plot
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
