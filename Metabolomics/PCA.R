# Load required libraries
library(tidyverse)
library(ggplot2)
library(scales)
#----------------Cleaning the environment
rm(list=ls())
# Read the data
data <- read.csv("Data__after_Cleanup.csv")
metadata <- read.csv("Metadata.csv")

# Fix column names
colnames(data) <- sub("X","", colnames(data))
metadata <- subset(metadata,Gruops %in% c("Healthy","CD","UC"))
# Remove the first three columns (ID, m/z, RT) and convert to matrix
peak_data <- as.matrix(data[, -c(1:3)])

# Transpose so samples are rows and features are columns
peak_data_t <- t(peak_data)

# Create sample names
sample_names <- rownames(peak_data_t)
metadata$sample_clean <- gsub(".mzML$", "", metadata$filename)

Gruops <- metadata$Gruops

# Log transform and scale the data
scaled_data <- scale(log1p(peak_data_t))

# Perform PCA
pca_result <- prcomp(scaled_data)

# Create data frame with PCA results
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2]
)

# Add sample names as a column
pca_df$sample_clean <- sample_names
pca_df$Gruops <- Gruops
# Merge with metadata
pca_df <- merge(pca_df, metadata, by = "sample_clean", all.x = TRUE)

# Print diagnostic information
cat("Dimensions of PCA dataframe:", dim(pca_df), "\n")
cat("Number of NA values in Groups:", sum(is.na(pca_df$Gruops)), "\n")

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)
total_var <- round(sum(var_explained[1:2]) * 100, 1)

# Define custom colors for 'Gruops'
custom_colors <- c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")

# Create PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Gruops)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "PCA of Metabolomics Data",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Groups",
    shape = "Groups"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = custom_colors)

pca_plot

# Add confidence ellipses (95%)
pca_plot <- pca_plot + 
  stat_ellipse(level = 0.95, type = "t")

pca_plot

# Save the plot
ggsave("pca_plot_for_Metabolites.png", pca_plot, width = 10, height = 8, dpi = 300)

# Print variance information
cat("\nExplained variance:\n")
cat(paste("PC1:", pc1_var, "%\n"))
cat(paste("PC2:", pc2_var, "%\n"))
cat(paste("\nTotal variance explained by PC1 and PC2:", total_var, "%\n"))

# Print first few rows of final data to check
print(head(pca_df))
