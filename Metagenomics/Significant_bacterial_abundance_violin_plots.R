# Load required libraries
{library(tidyverse)
library(ggplot2)
library(reshape2)
library(stats)
library(ggpubr)
library(ggforce)  # For jittered points
library(patchwork)
}
# Clear environment
rm(list=ls())


# Read the correlation data to get unique bacteria
correlations <- read.csv("sterco-correllation.csv")
unique_bacteria <- unique(correlations$bacteria)

# Read the RNA data and metadata
rna_data <- read.csv("16RNA_results_renamed.csv")
metadata <- read.csv("Metadata.csv")

# Fix column names
colnames(rna_data) <- sub("X","", colnames(rna_data))
# Function to perform ANOVA and get p-value
get_anova_pvalue <- function(bacteria_name) {
  # Prepare data
  bacteria_row <- rna_data[rna_data$Name == bacteria_name, ]
  
  if(nrow(bacteria_row) == 0) {
    warning(paste("Bacteria", bacteria_name, "not found in RNA data"))
    return(NULL)
  }
  
  # Extract intensity columns (ending with _aq.mzML)
  intensity_cols <- grep("_aq.mzML$", colnames(bacteria_row), value = TRUE)
  
  # Prepare data for ANOVA
  abundance_data <- data.frame(
    abundance = as.numeric(bacteria_row[1, intensity_cols]),
    group = metadata$Gruops[match(gsub(".mzML", "", intensity_cols), 
                                  gsub(".mzML", "", metadata$filename))]
  )
  
  # Remove any NA groups
  abundance_data <- abundance_data[!is.na(abundance_data$group), ]
  
  # Perform one-way ANOVA
  anova_result <- aov(abundance ~ group, data = abundance_data)
  p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
  
  return(list(
    data = abundance_data,
    p_value = p_value
  ))
}

# Create violin plots for bacteria with significant differences
plot_bacteria_abundance <- function(bacteria_name) {
  # Get ANOVA results
  anova_results <- get_anova_pvalue(bacteria_name)
  
  if(is.null(anova_results) || anova_results$p_value >= 0.05) {
    return(NULL)
  }
  
  # Prepare data
  bacteria_data <- anova_results$data
  p_value <- anova_results$p_value
  
  # Create violin plot with individual points
  p <- ggplot(bacteria_data, aes(x = group, y = abundance, fill = group, color = group)) +
    geom_violin(trim = FALSE, alpha = 1) +  # Change alpha to 1 for full opacity
    geom_jitter(width = 0.2, size = 2, alpha = 1) +  # Change point alpha to 1
    labs(
      title = paste("Abundance of", bacteria_name),
      x = "Group",
      y = "Abundance %"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),  # Add white background
      plot.background = element_rect(fill = "white"),   # Add white plot background
      #panel.grid = element_blank()  # Remove grid lines for cleaner look
    ) +
    scale_fill_manual(values = c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")) +
    scale_color_manual(values = c("Healthy" = "black", "CD" = "black", "UC" = "black"))
  
  # Add p-value annotation
  p <- p + stat_compare_means(
    comparisons = list(
      c("Healthy", "CD"),
      c("Healthy", "UC"),
      c("CD", "UC")
    ),
    method = "t.test",
    label = "p.signif",  # This will add * notation
    hide.ns = TRUE,  # This will remove non-significant annotations 
  )
  
  # Save the plot as PNG
  ggsave(
    filename = paste0(gsub(" ", "_", bacteria_name), "_abundance_violin.png"), 
    plot = p, 
    width = 10,  # Increased width to accommodate annotations
    height = 7,
    dpi = 300
  )
  
  return(p)
}
setwd("P:/Omer/Bilirubin & IDB/drafts for paper/5.2.2025/last/Paper/ster_plots")
# Generate plots for bacteria with significant differences
plots <- lapply(unique_bacteria, plot_bacteria_abundance)

# Print names of successfully generated plots
successful_plots <- unique_bacteria[!sapply(plots, is.null)]
cat("Bacteria with significant abundance differences:\n")
print(successful_plots)





#==================saving the all the graphs together =============

# Load required libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stats)
library(ggpubr)
library(ggforce)  # For jittered points
library(patchwork)  # For combining plots

# Clear environment
rm(list=ls())

# Set working directory
setwd("P:/Omer/Bilirubin & IDB/drafts for paper/5.2.2025/last/Paper")

# Read the correlation data to get unique bacteria
correlations <- read.csv("corr with ster.csv")
unique_bacteria <- unique(correlations$bacteria)

# Read the RNA data and metadata
rna_data <- read.csv("16RNA_results_renamed%.csv")
metadata <- read.csv("NEWmetadata.csv")

# Fix column names
colnames(rna_data) <- sub("X","", colnames(rna_data))

# Function to perform ANOVA and get p-value
get_anova_pvalue <- function(bacteria_name) {
  # Prepare data
  bacteria_row <- rna_data[rna_data$Name == bacteria_name, ]
  
  if(nrow(bacteria_row) == 0) {
    warning(paste("Bacteria", bacteria_name, "not found in RNA data"))
    return(NULL)
  }
  
  # Extract intensity columns (ending with _aq.mzML)
  intensity_cols <- grep("_aq.mzML$", colnames(bacteria_row), value = TRUE)
  
  # Prepare data for ANOVA
  abundance_data <- data.frame(
    abundance = as.numeric(bacteria_row[1, intensity_cols]),
    group = metadata$Gruops[match(gsub(".mzML", "", intensity_cols), 
                                  gsub(".mzML", "", metadata$filename))]
  )
  
  # Remove any NA groups
  abundance_data <- abundance_data[!is.na(abundance_data$group), ]
  
  # Perform one-way ANOVA
  anova_result <- aov(abundance ~ group, data = abundance_data)
  p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
  
  return(list(
    data = abundance_data,
    p_value = p_value
  ))
}

# Function to create violin plot for a single bacteria
plot_bacteria_abundance <- function(bacteria_name) {
  # Get ANOVA results
  anova_results <- get_anova_pvalue(bacteria_name)
  
  if(is.null(anova_results) || anova_results$p_value >= 0.05) {
    return(NULL)
  }
  
  # Prepare data
  bacteria_data <- anova_results$data
  p_value <- anova_results$p_value
  
  # Create violin plot with individual points
  p <- ggplot(bacteria_data, aes(x = group, y = abundance, fill = group, color = group)) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_jitter(width = 0.2, size = 1, alpha = 1) +
    labs(
      title = bacteria_name,
      x = "Group",
      y = "Abundance %"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 7),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    scale_fill_manual(values = c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")) +
    scale_color_manual(values = c("Healthy" = "black", "CD" = "black", "UC" = "black"))
  
  # Add p-value annotation
  p <- p + stat_compare_means(
    comparisons = list(
      c("Healthy", "CD"),
      c("Healthy", "UC"),
      c("CD", "UC")
    ),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE,
    size = 2.5
  )
  
  return(p)
}

# Generate plots for bacteria with significant differences
plots <- lapply(unique_bacteria, plot_bacteria_abundance)

# Remove NULL values (non-significant plots)
plots <- plots[!sapply(plots, is.null)]

# Calculate number of rows and columns for the grid
n_plots <- length(plots)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots/n_cols)

# Combine all plots using patchwork
combined_plot <- wrap_plots(plots, ncol = n_cols)
setwd("P:/Omer/Bilirubin & IDB/drafts for paper/5.2.2025/last/Paper/ster_plots/ster unique")
# Save the combined plot
ggsave(
  filename = "combined_abundance_violin_plots.png",
  plot = combined_plot,
  width = 3 * n_cols,
  height = 3 * n_rows,
  dpi = 300
)

# Print names of successfully plotted bacteria
successful_plots <- unique_bacteria[!sapply(plots, is.null)]
cat("Bacteria with significant abundance differences:\n")
print(successful_plots)

# Print information about the combined plot
cat(sprintf("\nCreated combined plot with %d significant bacteria\n", n_plots))
cat(sprintf("Grid dimensions: %d rows x %d columns\n", n_rows, n_cols))


#======================# Load required libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stats)
library(ggpubr)
library(ggforce)
library(patchwork)

# Clear environment
rm(list=ls())

# Define the bacteria of interest
bacteria_list <- c(
  "Blautia_sp._SG-772",
  "Lachnospiraceae_bacterium_UBA4711",
  "Lachnospiraceae_bacterium_UBA5902",
  "Lachnospiraceae_bacterium_UBA6452",
  "Oscillibacter_sp."
)

# Read the RNA data and metadata
rna_data <- read.csv("16RNA_results_renamed%.csv")
metadata <- read.csv("NEWmetadata.csv")

# Fix column names
colnames(rna_data) <- sub("X","", colnames(rna_data))

# Same functions as before, just modified to handle these specific bacteria
get_anova_pvalue <- function(bacteria_name) {
  # [Previous function code remains the same]
}

plot_bacteria_abundance <- function(bacteria_name) {
  # [Previous function code remains the same]
}

# Generate plots only for our specific bacteria
plots <- lapply(bacteria_list, plot_bacteria_abundance)

# Remove NULL values (non-significant plots)
plots <- plots[!sapply(plots, is.null)]

# Arrange plots in a 2x3 grid (which will accommodate 5 plots nicely)
combined_plot <- wrap_plots(plots, ncol = 2)

# Save the combined plot
ggsave(
  filename = "combined_abundance_violin_plots.png",
  plot = combined_plot,
  width = 8,  # Adjusted for 2 columns
  height = 12, # Adjusted to accommodate all plots
  dpi = 300
)

# Print names of successfully plotted bacteria
successful_plots <- bacteria_list[!sapply(plots, is.null)]
cat("Bacteria with significant abundance differences:\n")
print(successful_plots)

#=====================


# Load required libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stats)
library(ggpubr)
library(ggforce)
library(patchwork)

# Clear environment
rm(list=ls())

# Define the bacteria of interest
bacteria_list <- c(
  "Blautia_sp._SG-772",
  "Lachnospiraceae_bacterium_UBA4711",
  "Lachnospiraceae_bacterium_UBA5902",
  "Lachnospiraceae_bacterium_UBA6452",
  "Oscillibacter_sp."
)

setwd("P:/Omer/Bilirubin & IDB/drafts for paper/5.2.2025/last/Paper")
# Read the RNA data and metadata
rna_data <- read.csv("16RNA_results_renamed%.csv")
metadata <- read.csv("NEWmetadata.csv")

# Fix column names
colnames(rna_data) <- sub("X","", colnames(rna_data))

# Same functions as before, just modified to handle these specific bacteria
get_anova_pvalue <- function(bacteria_name) {
  # [Previous function code remains the same]
}

plot_bacteria_abundance <- function(bacteria_name) {
  # [Previous function code remains the same]
}

# Generate plots only for our specific bacteria
plots <- lapply(bacteria_list, plot_bacteria_abundance)

# Remove NULL values (non-significant plots)
plots <- plots[!sapply(plots, is.null)]

# Arrange plots in a 2x3 grid (which will accommodate 5 plots nicely)
combined_plot <- wrap_plots(plots, ncol = 2)

# Save the combined plot
ggsave(
  filename = "combined_abundance_violin_plots.png",
  plot = combined_plot,
  width = 8,  # Adjusted for 2 columns
  height = 12, # Adjusted to accommodate all plots
  dpi = 300
)

# Print names of successfully plotted bacteria
successful_plots <- bacteria_list[!sapply(plots, is.null)]
cat("Bacteria with significant abundance differences:\n")
print(successful_plots)