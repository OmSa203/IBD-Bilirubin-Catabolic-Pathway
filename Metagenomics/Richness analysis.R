#==== Libraries and Functions====

#---- Load the libraries -----
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(dunn.test)
library(ggpubr)
#----Functions----

# --- Check for matching sample names ---
# This is CRITICAL. Sample IDs (column names in otu_data, row names in metadata) MUST match.
# Sort to ensure consistent order if needed, but the names themselves are most important.
check_sample_names <- function(otu_data, metadata) {
  otu_cols <- colnames(otu_data)
  meta_rows <- rownames(metadata)
  
  missing_in_metadata <- setdiff(otu_cols, meta_rows)
  missing_in_otu_data <- setdiff(meta_rows, otu_cols)
  
  if (length(missing_in_metadata) > 0 || length(missing_in_otu_data) > 0) {
    stop_message <- "Sample names in abundance table columns and metadata row names do not match!\n"
    if (length(missing_in_metadata) > 0) {
      stop_message <- paste0(stop_message, "Missing in metadata: ", paste(missing_in_metadata, collapse = ", "), "\n")
    }
    if (length(missing_in_otu_data) > 0) {
      stop_message <- paste0(stop_message, "Missing in otu_data: ", paste(missing_in_otu_data, collapse = ", "), "\n")
    }
    stop(stop_message)
  } else {
    message("Sample names match. Proceeding...")
  }
}

# --- cmopers between two groups and save an md file---
richness_stats_and_report <- function(richness_metrics_data,grouping_var, output_filename = "richness_statistical_results.md") {
  
  # Input validation
  if (!inherits(richness_metrics_data, "data.frame")) {
    stop("richness_metrics_data must be a data frame.")
  }
  if (!grouping_var %in% colnames(richness_metrics_data)) {
    stop(paste0("Grouping variable '", grouping_var, "' not found in richness_metrics_data."))
  }
  
  # Ensure richness metrics columns exist
  required_metrics <- c("Observed", "Chao1", "ACE")
  if (!all(required_metrics %in% colnames(richness_metrics_data))) {
    missing_metrics <- setdiff(required_metrics, colnames(richness_metrics_data))
    stop(paste0("Missing required richness metrics columns: ", paste(missing_metrics, collapse = ", ")))
  }
  
  # Extract the grouping variable vector
  groups <- richness_metrics_data[[grouping_var]]
  num_unique_groups <- length(unique(groups[!is.na(groups)])) # Exclude NA for count
  print(num_unique_groups)
  # Start capturing output to the Markdown file and console
  # 'append = TRUE' adds to the file, 'split = TRUE' sends to console AND file
  sink(output_filename, append = TRUE, split = TRUE)
  
  # Add a Markdown heading and timestamp to the file
  cat(paste0("# Statistical Analysis Results for '", grouping_var, "'\n\n"))
  cat(paste0("Analysis performed on: ", Sys.time(), "\n\n"))
  cat("```R\n") # Start a Markdown code block for R syntax highlighting
  
  if (num_unique_groups == 2) {
    message(paste("Performing Wilcoxon Rank-Sum tests for", grouping_var))
    cat(paste("Performing Wilcoxon Rank-Sum tests for", grouping_var, "\n"))
    
    # Observed Richness
    wilcox_obs <- tryCatch({
      wilcox.test(Observed ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Wilcoxon test for Observed Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(wilcox_obs)) {
      print(paste("Observed Richness (Wilcoxon): p-value =", wilcox_obs$p.value))
      cat(paste("Observed Richness (Wilcoxon): p-value =", wilcox_obs$p.value, "\n"))
    } else {
      cat("Observed Richness (Wilcoxon): Test could not be performed.\n")
    }
    
    
    # Chao1 Richness
    wilcox_chao1 <- tryCatch({
      wilcox.test(Chao1 ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Wilcoxon test for Chao1 Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(wilcox_chao1)) {
      print(paste("Chao1 Richness (Wilcoxon): p-value =", wilcox_chao1$p.value))
      cat(paste("Chao1 Richness (Wilcoxon): p-value =", wilcox_chao1$p.value, "\n"))
    } else {
      cat("Chao1 Richness (Wilcoxon): Test could not be performed.\n")
    }
    
    # ACE Richness
    wilcox_ace <- tryCatch({
      wilcox.test(ACE ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Wilcoxon test for ACE Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(wilcox_ace)) {
      print(paste("ACE Richness (Wilcoxon): p-value =", wilcox_ace$p.value))
      cat(paste("ACE Richness (Wilcoxon): p-value =", wilcox_ace$p.value, "\n"))
    } else {
      cat("ACE Richness (Wilcoxon): Test could not be performed.\n")
    }
    
  } else if (num_unique_groups > 2) {
    # --- B. If you have MORE THAN TWO groups ---
    message(paste("Performing Kruskal-Wallis tests for", grouping_var))
    cat(paste("Performing Kruskal-Wallis tests for", grouping_var, "\n"))
    
    # Observed Richness
    kruskal_obs <- tryCatch({
      kruskal.test(Observed ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Kruskal-Wallis test for Observed Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(kruskal_obs)) {
      print(paste("Observed Richness (Kruskal-Wallis): p-value =", kruskal_obs$p.value))
      cat(paste("Observed Richness (Kruskal-Wallis): p-value =", kruskal_obs$p.value, "\n"))
      if (kruskal_obs$p.value < 0.05) {
        message("  Performing Dunn's post-hoc test for Observed Richness:")
        cat("  Performing Dunn's post-hoc test for Observed Richness:\n")
        dunn_obs_output <- capture.output(dunn.test(richness_metrics_data$Observed, groups, method = "bonferroni"))
        cat(paste(dunn_obs_output, collapse = "\n"), "\n")
      }
    } else {
      cat("Observed Richness (Kruskal-Wallis): Test could not be performed.\n")
    }
    
    # Chao1 Richness
    kruskal_chao1 <- tryCatch({
      kruskal.test(Chao1 ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Kruskal-Wallis test for Chao1 Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(kruskal_chao1)) {
      print(paste("Chao1 Richness (Kruskal-Wallis): p-value =", kruskal_chao1$p.value))
      cat(paste("Chao1 Richness (Kruskal-Wallis): p-value =", kruskal_chao1$p.value, "\n"))
      if (kruskal_chao1$p.value < 0.05) {
        message("  Performing Dunn's post-hoc test for Chao1 Richness:")
        cat("  Performing Dunn's post-hoc test for Chao1 Richness:\n")
        dunn_chao1_output <- capture.output(dunn.test(richness_metrics_data$Chao1, groups, method = "bonferroni"))
        cat(paste(dunn_chao1_output, collapse = "\n"), "\n")
      }
    } else {
      cat("Chao1 Richness (Kruskal-Wallis): Test could not be performed.\n")
    }
    
    
    # ACE Richness
    kruskal_ace <- tryCatch({
      kruskal.test(ACE ~ groups, data = richness_metrics_data)
    }, error = function(e) {
      message("Error with Kruskal-Wallis test for ACE Richness: ", e$message)
      return(NULL)
    })
    if (!is.null(kruskal_ace)) {
      print(paste("ACE Richness (Kruskal-Wallis): p-value =", kruskal_ace$p.value))
      cat(paste("ACE Richness (Kruskal-Wallis): p-value =", kruskal_ace$p.value, "\n"))
      if (kruskal_ace$p.value < 0.05) {
        message("  Performing Dunn's post-hoc test for ACE Richness:")
        cat("  Performing Dunn's post-hoc test for ACE Richness:\n")
        dunn_ace_output <- capture.output(dunn.test(richness_metrics_data$ACE, groups, method = "bonferroni"))
        cat(paste(dunn_ace_output, collapse = "\n"), "\n")
      }
    } else {
      cat("ACE Richness (Kruskal-Wallis): Test could not be performed.\n")
    }
    
  } else {
    message("Warning: Grouping variable '", grouping_var, "' has fewer than 2 unique levels or contains only NAs after filtering.")
    cat("Warning: Grouping variable '", grouping_var, "' has fewer than 2 unique levels or contains only NAs after filtering.\n")
  }
  
  cat("```\n") # End the code block in Markdown
  cat("\n--- End of Analysis ---\n")
  
  # --- Stop capturing output ---
  sink()
  message(paste("Output also written to", output_filename))
}

plot_alpha_richness <- function(richness_metrics_data, grouping_var, save_plots = FALSE, save_dir = NULL) {
  
  # Plot for Observed Richness
  p_observed <- ggplot(alpha_df, aes_string(x = grouping_variable, y = "Observed", fill = grouping_variable)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) + # Show individual data points
    labs(
      title = paste("OTUs"),
      y = "Observed Bacteria (ASVs/OTUs)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none", # Adjust legend as needed
      plot.title = element_text(hjust = 0.5) # Center the plot title
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1"))+
    # Add significance bracket
    stat_compare_means(
      comparisons = list(
        c("Healthy", "CD"),
        c("Healthy", "UC")
        #c("CD", "UC")
      ),
      method = "wilcox.test",
      p.adjust.method = "bonferroni",
      label = "p.signif" ,# This will convert p-values to stars
      hide.ns = FALSE  # This will remove non-significant annotations 
      #bracket = TRUE  # This removes the brackets
    )
  
  # Plot for Chao1 Richness
  p_chao1 <- ggplot(alpha_df, aes_string(x = grouping_variable, y = "Chao1", fill = grouping_variable)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    labs(
      title = paste("Chao1 Richness"),
      y = "Estimated Chao1 Richness"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none", # Adjust legend as needed
      plot.title = element_text(hjust = 0.5) # Center the plot title
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1"))+
    # Add significance bracket
    stat_compare_means(
      comparisons = list(
        c("Healthy", "CD"),
        c("Healthy", "UC")
        #c("CD", "UC")
      ),
      method = "wilcox.test",
      p.adjust.method = "bonferroni",
      label = "p.signif" ,# This will convert p-values to stars
      hide.ns = FALSE  # This will remove non-significant annotations 
      #bracket = FALSE  # This removes the brackets
    )
  
  # Plot for ACE Richness
  p_ace <- ggplot(alpha_df, aes_string(x = grouping_variable, y = "ACE", fill = grouping_variable)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    labs(
      title = paste("ACE Richness"),
      y = "Estimated ACE Richness"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none", # Adjust legend as needed
      plot.title = element_text(hjust = 0.5) # Center the plot title
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#6d8abf", "CD" = "#ccbb9c", "UC" = "#ebe8e1")) +
    # Add significance bracket
    stat_compare_means(
      comparisons = list(
        c("Healthy", "CD"),
        c("Healthy", "UC")
        #c("CD", "UC")
      ),
      method = "wilcox.test",
      p.adjust.method = "bonferroni",
      label = "p.signif" ,# This will convert p-values to stars
      hide.ns = FALSE  # This will remove non-significant annotations 
      #bracket = TRUE  # This removes the brackets
    )
  
  
  # --- Automatic Saving Logic ---
  if (save_plots) {
    if (!is.null(save_dir)) {
      # Create the directory if it doesn't exist
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
      # Set the directory for saving
      current_save_path <- save_dir
    } else {
      # Use current working directory if save_dir is NULL
      current_save_path <- getwd()
      message(paste("Saving plots to current working directory:", current_save_path))
    }
    
    # Save each plot
    ggsave(filename = paste0(current_save_path, "/Observed_Richness.png"), plot = p_observed, width = 7, height = 5, units = "in", dpi = 300)
    ggsave(filename = paste0(current_save_path, "/Chao1_Richness.png"), plot = p_chao1, width = 7, height = 5, units = "in", dpi = 300)
    ggsave(filename = paste0(current_save_path, "/ACE_Richness.png"), plot = p_ace, width = 7, height = 5, units = "in", dpi = 300)
    
    message("Plots saved successfully.")
  }
  
  # Return the plots as a list
  return(list(
    p_observed = p_observed,
    p_chao1 = p_chao1,
    p_ace = p_ace
  ))
}


#==== Preparing Data =====

# ---- Import Bacteria Abundance Data & Metadata ----
# IMPORTANT: Ensure your bacteria are in ROWS and samples in COLUMNS.
# phyloseq expects taxa (bacteria) in rows for otu_table(taxa_are_rows = TRUE).
# If your samples are rows and bacteria are columns, you'll need to transpose:
# otu_data <- t(otu_data)

otu_data <- read.csv("Abundance_Results%.csv", row.names = 1, header = TRUE)
colnames(otu_data) <- sub("X","", colnames(otu_data)) #If needed 

# --- Import Metadata ---
# Read the metadata table. Assuming first column is Sample ID.
metadata <- read.csv("Metadata.csv", row.names = 1, header = TRUE)

# Find differences
check_sample_names(otu_data, metadata)

#If missing names are found and it for example: Solvent.mzML found in the metadata but not in the samples table, I deleted this row. 
metadata <- metadata %>%
  filter(row.names(metadata) != "Solvent.mzML")

#Second check after deletion
check_sample_names(otu_data, metadata)


#----- Create a phyloseq Object -----
# This step integrates your data into a single, convenient object for downstream analysis.


# Convert data frames to phyloseq components
OTU <- otu_table(otu_data, taxa_are_rows = TRUE)
SAMP <- sample_data(metadata)

# Create the phyloseq object
# If you also have a taxonomy table for your bacteria, you'd add it here:
# TAX <- tax_table(as.matrix(your_taxonomy_df))
# ps_obj <- phyloseq(OTU, SAMP, TAX)

ps_obj <- phyloseq(OTU, SAMP)

# View a summary of your phyloseq object
ps_obj


#===== Richness Analysis ====


#----Rarefaction (Normalization for Richness)----
#Rarefaction is essential for comparing richness across samples with different sequencing depths. 
#It subsamples each sample to the same number of reads (the minimum observed library size).


# Check library sizes (total counts per sample)
sample_sums(ps_obj)

# Determine the minimum library size (rarefaction depth)
min_lib_size <- min(sample_sums(ps_obj))
message(paste("Minimum library size:", min_lib_size))

# Rarefy your phyloseq object
# Set a random seed for reproducibility
ps_obj_rare <- rarefy_even_depth(ps_obj, sample.size = min_lib_size, rngseed = 123)

# Verify that all samples now have the same library size
sample_sums(ps_obj_rare)


#----Calculate Richness Metrics----

alpha_diversity_metrics <- estimate_richness(ps_obj_rare, measures = c("Observed", "Chao1", "ACE"))

# View the first few rows of the calculated metrics
head(alpha_diversity_metrics)
rownames(alpha_diversity_metrics) <- sub("X","", rownames(alpha_diversity_metrics)) #If needed
head(alpha_diversity_metrics) #check if the "X" deleted 

# Add these metrics back to your phyloseq object's sample data for easy access
sample_data(ps_obj_rare) <- cbind(sample_data(ps_obj_rare), alpha_diversity_metrics)

# You can also create a separate data frame for easier manipulation and plotting
alpha_df <- data.frame(sample_data(ps_obj_rare))

# View the combined data frame
head(alpha_df)
write.csv(alpha_df,"alpha_df.csv", row.names = FALSE)

#----Statistical Analysis (Comparing Richness Between Groups)----

# Identify your grouping variable (replace 'Group' with your actual column name)
grouping_variable <- "Groups" # e.g., "Treatment", "PatientType", etc.

# Check the unique levels of your grouping variable
unique(alpha_df[[grouping_variable]])

#the analysis and the report 
richness_stats_and_report(alpha_df,grouping_variable)

#----Visualision----
plot_alpha_richness(alpha_df,grouping_variable,TRUE)

