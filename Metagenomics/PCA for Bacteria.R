# Load required libraries
library(tidyverse)
library(ggplot2)
library(scales)

# Clear environment
rm(list=ls())

# Read the data
data <- read.csv("Abundance_Results.csv")
metadata <- read.csv("NEWmetadata.csv")

# Process metadata 
metadata <- subset(metadata, Gruops %in% c("Healthy","CD","UC"))

# Prepare bacterial abundance data
abundance_data <- data[, -1]  # Remove Name column
bacteria_names <- data$Name
rownames(abundance_data) <- bacteria_names

# Transpose data (samples as rows, bacteria as columns)
abundance_data_t <- t(abundance_data)

# Fix sample names - remove 'X' prefix from abundance data
rownames(abundance_data_t) <- sub("^X", "", rownames(abundance_data_t))

# Print first few sample names to verify matching
cat("Sample names in abundance data (first 5):", head(rownames(abundance_data_t), 5), "\n")
cat("Sample names in metadata (first 5):", head(metadata$filename, 5), "\n")

# Log transform the data without scaling
transformed_data <- log1p(abundance_data_t)

# Perform PCA
pca_result <- prcomp(transformed_data, scale. = FALSE)

# Match groups to samples
sample_groups <- metadata$Gruops[match(rownames(abundance_data_t), metadata$filename)]

# Create data frame with PCA results
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Gruops = sample_groups  # Use the matched groups
)

# Define custom colors
custom_colors <- c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")

# Create PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Gruops)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_bw() +
  labs(
    title = "PCA of Bacterial Community Data",
    x = paste0("PC1 (", round(100 * pca_result$sdev[1]^2/sum(pca_result$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_result$sdev[2]^2/sum(pca_result$sdev^2), 1), "%)"),
    color = "Groups"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = custom_colors)

# Print group summary
cat("\nGroup assignment summary:\n")
print(table(sample_groups, useNA = "always"))
pca_plot
# Save plot
ggsave("pca_plot_for_Bacterial_Community.png", pca_plot, width = 10, height = 8, dpi = 300)
