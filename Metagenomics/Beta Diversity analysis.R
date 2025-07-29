# Beta Diversity Analysis for IBD Microbiome Data
# This script analyzes metagenomic data to examine differences in bacterial community
# composition between healthy controls, Crohn's disease, and ulcerative colitis patients

# Load required libraries
library(tidyverse)   # For data manipulation and visualization
library(vegan)       # For diversity analysis and ordination methods
library(ggplot2)     # For creating publication-quality plots
library(ggrepel)     # For better label placement in plots
library(RColorBrewer) # For color palettes
library(gridExtra)   # For arranging multiple plots
library(patchwork)   # For combining plots
library(ggsci)       # For scientific journal color palettes

# Set seed for reproducibility
set.seed(123)

# ======= DATA IMPORT AND PREPARATION =======

# Import the bacterial abundance data
microbial_data <- read.csv("Abundance_Results.csv", row.names = 1)
# Fix column names
colnames(microbial_data) <- sub("X","", colnames(microbial_data))
colnames(microbial_data) <- sub(".mzML$","", colnames(microbial_data))
# Import metadata containing group information
metadata <- read.csv("Metadata.csv")

# Clean up sample names in metadata and filter for relevant groups
metadata <- metadata %>%
  mutate(sample_clean = gsub(".mzML$", "", filename)) %>%
  filter(Gruops %in% c("Healthy", "CD", "UC"))

# Transpose the abundance data to have samples as rows and bacterial species as columns
# This is the required format for most diversity analyses
abundance_matrix <- t(as.matrix(microbial_data[, -1]))  # Remove the first column which contains species names
rownames(abundance_matrix) <- colnames(microbial_data)[-1]  # Set sample names as row names

# Match abundance data with metadata
# Only keep samples that have metadata and vice versa
common_samples <- intersect(rownames(abundance_matrix), metadata$sample_clean)
abundance_matrix <- abundance_matrix[common_samples, ]
metadata <- metadata[match(common_samples, metadata$sample_clean), ]

# Check if data dimensions match
cat("Dimensions of abundance matrix:", dim(abundance_matrix), "\n")
cat("Number of samples in metadata:", nrow(metadata), "\n")

# Add a check to ensure we have sufficient data for analysis
if(nrow(abundance_matrix) < 10 || ncol(abundance_matrix) < 10) {
  stop("Not enough data for meaningful analysis (< 10 samples or < 10 species)")
}

# ======= DATA TRANSFORMATION AND NORMALIZATION =======

# Apply a log transformation to reduce the effect of highly abundant species
# Add a small constant to avoid log(0)
epsilon <- min(abundance_matrix[abundance_matrix > 0]) / 2
log_abundance <- log1p(abundance_matrix)

# Normalize data (relative abundance - converting to percentages)
relative_abundance <- sweep(abundance_matrix, 1, rowSums(abundance_matrix), "/") * 100

# Check for samples with very low total counts that might skew results
low_count_samples <- which(rowSums(abundance_matrix) < quantile(rowSums(abundance_matrix), 0.05))
if(length(low_count_samples) > 0) {
  cat("Warning: Samples with very low counts detected:", rownames(abundance_matrix)[low_count_samples], "\n")
}

# ======= BETA DIVERSITY CALCULATION =======

# Calculate Bray-Curtis dissimilarity matrix
bc_dist <- vegdist(relative_abundance, method = "bray")

# Calculate Jaccard dissimilarity (presence/absence data)
jaccard_dist <- vegdist(relative_abundance > 0, method = "jaccard")

# Calculate Euclidean distance for comparison
euclidean_dist <- vegdist(log_abundance, method = "euclidean")

# Export distance matrices
write.csv(as.matrix(bc_dist), "bray_curtis_distance_matrix.csv")
write.csv(as.matrix(jaccard_dist), "jaccard_distance_matrix.csv")

# ======= ORDINATION METHODS =======

# Principal Coordinates Analysis (PCoA) using Bray-Curtis distances
pcoa_bc <- cmdscale(bc_dist, k = 5, eig = TRUE)
pcoa_bc_var <- round(100 * pcoa_bc$eig / sum(pcoa_bc$eig), 1)

# Create a dataframe for PCoA results
pcoa_bc_df <- as.data.frame(pcoa_bc$points)
colnames(pcoa_bc_df) <- paste0("PCo", 1:ncol(pcoa_bc_df))
pcoa_bc_df$sample_id <- rownames(pcoa_bc_df)
pcoa_bc_df <- merge(pcoa_bc_df, metadata, by.x = "sample_id", by.y = "sample_clean")

# Calculate centroid positions for each group
centroids <- pcoa_bc_df %>%
  group_by(Gruops) %>%
  summarize(PCo1 = mean(PCo1),
            PCo2 = mean(PCo2),
            .groups = 'drop')

# Non-metric Multidimensional Scaling (NMDS) using Bray-Curtis distances
# Try several dimensions to find the best stress value
best_nmds <- NULL
best_stress <- Inf
best_k <- 0

for(k in 2:5) {
  try({
    nmds_attempt <- metaMDS(bc_dist, k = k, trymax = 50, autotransform = FALSE)
    if(nmds_attempt$stress < best_stress) {
      best_nmds <- nmds_attempt
      best_stress <- nmds_attempt$stress
      best_k <- k
    }
  }, silent = TRUE)
}

if(is.null(best_nmds)) {
  warning("NMDS did not converge with any of the attempted dimensions. Using 2D NMDS.")
  nmds_result <- metaMDS(bc_dist, k = 2, trymax = 100, autotransform = FALSE)
} else {
  nmds_result <- best_nmds
  cat("Best NMDS results with k =", best_k, "dimensions (stress =", best_stress, ")\n")
}

# Create a dataframe for NMDS results
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- paste0("NMDS", 1:ncol(nmds_df))
nmds_df$sample_id <- rownames(nmds_df)
nmds_df <- merge(nmds_df, metadata, by.x = "sample_id", by.y = "sample_clean")

# Export ordination coordinates
write.csv(pcoa_bc_df, "pcoa_coordinates.csv")
write.csv(nmds_df, "nmds_coordinates.csv")

# ======= STATISTICAL TESTS =======

# PERMANOVA test (Permutational Multivariate Analysis of Variance)
# Tests if the centroids of groups differ
permanova_result <- adonis2(bc_dist ~ Gruops, data = metadata, permutations = 999)
print("PERMANOVA results for Bray-Curtis distances:")
print(permanova_result)

# ANOSIM (Analysis of Similarities)
# Tests if the between-group distances are greater than within-group distances
anosim_result <- anosim(bc_dist, metadata$Gruops, permutations = 999)
print("ANOSIM results for Bray-Curtis distances:")
print(anosim_result)

# ======= VISUALIZATION =======

# Define a consistent color scheme for the groups
custom_colors <- c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")

# PCoA plot with confidence ellipses (95%)
pcoa_plot <- ggplot(pcoa_bc_df, aes(x = PCo1, y = PCo2, color = Gruops)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, type = "t") +
  geom_text_repel(
    data = centroids,
    aes(x = PCo1, y = PCo2, label = Gruops),
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.5,
    force = 1,
    segment.color = "grey50"
  ) +
  labs(
    title = "Principal Coordinates Analysis (PCoA)",
    x = paste0("PCo1 (", pcoa_bc_var[1], "% of variance)"),
    y = paste0("PCo2 (", pcoa_bc_var[2], "% of variance)"),
    color = "Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = custom_colors)

# NMDS plot with confidence ellipses (95%)
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Gruops)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    title = paste0("NMDS Ordination (stress = ", round(nmds_result$stress, 3), ")"),
    x = "NMDS1",
    y = "NMDS2",
    color = "Group"
  ) +
  annotate(
    "text",
    x = min(nmds_df$NMDS1), y = max(nmds_df$NMDS2),
    label = paste("Stress =", round(nmds_result$stress, 3)),
    hjust = 0, vjust = 1
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = custom_colors)

# Add PERMANOVA and ANOSIM statistics to plots
permanova_pval <- permanova_result$`Pr(>F)`[1]
anosim_pval <- anosim_result$signif

pcoa_plot <- pcoa_plot +
  annotate(
    "text",
    x = min(pcoa_bc_df$PCo1),
    y = min(pcoa_bc_df$PCo2),
    label = sprintf("PERMANOVA: p = %.4f\nANOSIM: p = %.4f", permanova_pval, anosim_pval),
    hjust = 0.1, vjust = 0.9
  )

# Combine both plots into a single figure
combined_plot <- pcoa_plot + nmds_plot + plot_layout(ncol = 2)

# Combine plots with statistics on the side
# First convert the ggplot objects to grobs
pcoa_grob <- ggplotGrob(pcoa_plot)
nmds_grob <- ggplotGrob(nmds_plot)

# Arrange the plots with the stats text
combined_layout <- grid.arrange(
  pcoa_grob, nmds_grob, stats_grob,
  ncol = 3, 
  widths = c(4, 4, 1.5),
  top = "Beta Diversity Analysis of IBD Microbiome"
)
# Save plots
ggsave("pcoa_plot.pdf", pcoa_plot, width = 8, height = 6, dpi = 300)
ggsave("nmds_plot.pdf", nmds_plot, width = 8, height = 6, dpi = 300)
ggsave("combined_ordination_plots.pdf", combined_plot, width = 16, height = 6, dpi = 300)
ggsave("pcoa_plot.png", pcoa_plot, width = 8, height = 6, dpi = 300)
ggsave("nmds_plot.png", nmds_plot, width = 8, height = 6, dpi = 300)
ggsave("combined_ordination_plots.png", combined_plot, width = 16, height = 6, dpi = 300)

# ======= ADDITIONAL ANALYSES =======

# Calculate distance to group centroid (beta dispersion)
# This tests if the within-group variation differs between groups
beta_dispersion <- betadisper(bc_dist, metadata$Gruops)
disp_anova <- anova(beta_dispersion)
print("Beta dispersion ANOVA results:")
print(disp_anova)

# Tukey's post-hoc test for beta dispersion
tukey_result <- TukeyHSD(beta_dispersion)
print("Tukey's HSD test for beta dispersion:")
print(tukey_result)

# Plot beta dispersion
dispersion_plot <- plot(beta_dispersion, hull = FALSE, ellipse = TRUE, main = "Beta Dispersion")

# Save dispersion plot
pdf("beta_dispersion_plot.pdf", width = 8, height = 6)
plot(beta_dispersion, hull = FALSE, ellipse = TRUE, main = "Beta Dispersion")
dev.off()

# ======= REPORTING FINAL RESULTS =======

# Create a summary dataframe for the most important results
analysis_summary <- data.frame(
  Analysis = c("PERMANOVA (Bray-Curtis)", "ANOSIM (Bray-Curtis)", "Beta Dispersion"),
  Statistic = c(
    paste("F =", round(permanova_result$F[1], 3)),
    paste("R =", round(anosim_result$statistic, 3)),
    paste("F =", round(disp_anova$F[1], 3))
  ),
  P_value = c(
    permanova_pval,
    anosim_pval,
    disp_anova$`Pr(>F)`[1]
  ),
  Interpretation = c(
    ifelse(permanova_pval < 0.05, "Group centroids differ significantly", "No significant difference between group centroids"),
    ifelse(anosim_pval < 0.05, "Between-group distances significantly greater than within-group", "No significant difference in between vs within group distances"),
    ifelse(disp_anova$`Pr(>F)`[1] < 0.05, "Groups have significantly different dispersions", "No significant difference in group dispersions")
  )
)

write.csv(analysis_summary, "beta_diversity_analysis_summary.csv", row.names = FALSE)

cat("\n===== Summary of Beta Diversity Analysis =====\n")
print(analysis_summary)
cat("\nAnalysis complete. Results and plots have been saved to the working directory.\n")

# ======= SESSION INFO (for reproducibility) =======

session_info <- sessionInfo()
write_lines(capture.output(session_info), "beta_diversity_session_info.txt")
print(session_info)