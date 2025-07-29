# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

# Clear environment
rm(list=ls())

# Read and prepare data
data <- read.csv("Abundance_Results.csv")
metadata <- read.csv("Metadata.csv")
# Fix column names
colnames(data) <- sub("X","", colnames(data))
# Clean metadata
metadata <- subset(metadata, Gruops %in% c("Healthy", "CD", "UC"))

# Prepare abundance matrix
abundance_matrix <- as.matrix(data[, -1])  # Remove Name column

rownames(abundance_matrix) <- data$Name

# Calculate relative abundances
rel_abundance <- t(apply(abundance_matrix, 1, function(x) x/sum(x)))

# Get top 10 most abundant bacteria
mean_abundance <- rowMeans(rel_abundance)
top_bacteria <- names(sort(mean_abundance, decreasing = TRUE))[1:10]
top_abundance <- rel_abundance[top_bacteria, ]

# Transpose for diversity calculations
abundance_for_diversity <- t(abundance_matrix)

# Calculate diversity indices
shannon_div <- diversity(abundance_for_diversity, index = "shannon")
simpson_div <- diversity(abundance_for_diversity, index = "simpson")

# Create diversity dataframe
div_df <- data.frame(
  Sample = names(shannon_div),
  Shannon = shannon_div,
  Simpson = simpson_div
)

# Clean sample names
div_df$Sample <- gsub(".mzML$", "", div_df$Sample)
metadata$filename <- gsub(".mzML$", "", metadata$filename)
# Merge with metadata
div_df <- merge(div_df, metadata, by.x = "Sample", by.y = "filename")

# Define custom colors
custom_colors <- c(
  "Healthy" = "#6d8abf", 
  "CD" = "#ccbb9c", 
  "UC" = "#ebe8e1"
)

#--------------Shannon Diversity 

# Perform Kruskal-Wallis test
kw_test <- kruskal.test(Shannon ~ Gruops, data = div_df)

# Perform pairwise Wilcoxon tests
pw_test <- pairwise.wilcox.test(div_df$Shannon, div_df$Gruops, p.adjust.method = "bonferroni")

# Function to convert p-values to significance stars
p_to_stars <- function(p_value) {
  if (is.na(p_value)) return("ns")
  if (p_value <= 0.0001) return("****")
  if (p_value <= 0.001) return("***")
  if (p_value <= 0.01) return("**")
  if (p_value <= 0.05) return("*")
  return()
}

shannon_plot <- ggplot(div_df, aes(x = Gruops, y = Shannon, fill = Gruops)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  labs(
    title = "Shannon Diversity Index by Group",
    x = "Groups",
    y = "Shannon Diversity Index"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  # Add significance bracket
  stat_compare_means(
    comparisons = list(
      c("Healthy", "CD"),
      c("Healthy", "UC"),
      c("CD", "UC")
    ),
    method = "wilcox.test",
    p.adjust.method = "bonferroni",
    label = "p.signif" ,# This will convert p-values to stars
    hide.ns = TRUE,  # This will remove non-significant annotations 
    #bracket = FALSE  # This removes the brackets
  )

shannon_plot

#----------Abunance Bacteria

# Prepare abundance data for plotting
top_abundance_df <- as.data.frame(t(top_abundance))
top_abundance_df$Sample <- rownames(top_abundance_df)
top_abundance_long <- pivot_longer(top_abundance_df, 
                                   cols = -Sample,
                                   names_to = "Bacteria", 
                                   values_to = "Abundance")

# Add group information
top_abundance_long$Sample <- gsub(".mzML$", "", top_abundance_long$Sample)
top_abundance_long <- merge(top_abundance_long, metadata, 
                            by.x = "Sample", by.y = "filename")
# Function to perform Kruskal-Wallis test for each bacteria
perform_kruskal_test <- function(data) {
  kruskal.test(Abundance ~ Gruops, data = data)
}

# Perform tests for each bacteria
significance_tests <- top_abundance_long %>%
  group_by(Bacteria) %>%
  summarise(
    p_value = perform_kruskal_test(cur_data())$p.value
  )

# Add significance to abundance plot
abundance_plot <- ggplot(top_abundance_long, 
                         aes(x = Gruops, y = Abundance, fill = Bacteria)) +
  geom_boxplot(position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  labs(
    title = "Top 10 Bacteria Relative Abundance",
    x = "Groups",
    y = "Relative Abundance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  # Add significance brackets for each bacteria
  stat_compare_means(
    aes(group = Bacteria),
    label = "p.signif",
    method = "kruskal.test",
    bracket.size = 10
  )

abundance_plot

# Perform PERMANOVA test
dist_matrix <- vegdist(abundance_for_diversity, method = "bray")
permanova_result <- adonis2(dist_matrix ~ metadata$Gruops)

# Print results
print("PERMANOVA Results:")
print(permanova_result)

# Perform pairwise comparisons for diversity
print("Pairwise Wilcoxon test results for Shannon diversity:")
print(pairwise.wilcox.test(div_df$Shannon, div_df$Gruops, 
                           p.adjust.method = "bonferroni"))

# Print summary statistics for diversity indices
print("Summary statistics for Shannon diversity by group:")
print(aggregate(Shannon ~ Gruops, data = div_df, 
                FUN = function(x) c(mean = mean(x), sd = sd(x))))


#---------------Save plots
ggsave("shannon_diversity.png", shannon_plot, width = 8, height = 6)
ggsave("relative_abundance.png", abundance_plot, width = 15, height = 10)
