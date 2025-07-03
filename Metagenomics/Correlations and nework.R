library(tidyverse)
library(tidygraph)
library(ggraph)
library(reshape2)

#----------------Cleaning the environment
rm(list=ls())

# Read the data
rna_data <- read.csv("16RNA_results_renamed.csv")
metab_data <- read.csv("Data__after_Cleanup.csv")
bilirubin_analogs <- read.csv("Bilirubin's analogs.csv")

# Define bacteria of interest
target_bacteria <- rna_data$Name

# Filter RNA data for target bacteria
bacteria_indices <- grep(paste(target_bacteria, collapse="|"), rna_data$Name)
bacteria_data <- rna_data[bacteria_indices,]

# Function to group and merge features
merge_features <- function(metab_data, target_mz, ppm_tolerance, rt_tolerance = 0.5) {
  # Calculate mass error for this metabolite
  mass_error <- (ppm_tolerance * target_mz) / 1000000
  
  # Find features matching this mass
  matching_indices <- which(abs(metab_data$row.m.z - target_mz) < mass_error)
  
  if(length(matching_indices) == 0) return(NULL)
  
  matching_features <- metab_data[matching_indices,]
  
  # Sort by retention time
  matching_features <- matching_features[order(matching_features$row.retention.time),]
  
  # Group by retention time windows
  rt_groups <- list()
  current_group <- 1
  rt_groups[[1]] <- matching_features[1,]
  
  if(nrow(matching_features) > 1) {
    for(i in 2:nrow(matching_features)) {
      if(abs(matching_features$row.retention.time[i] - matching_features$row.retention.time[i-1]) <= rt_tolerance) {
        # Add to current group
        rt_groups[[current_group]] <- rbind(rt_groups[[current_group]], matching_features[i,])
      } else {
        # Start new group
        current_group <- current_group + 1
        rt_groups[[current_group]] <- matching_features[i,]
      }
    }
  }
  
  # Merge features within each group
  merged_features <- lapply(rt_groups, function(group) {
    # Get peak area columns
    intensity_cols <- grep("_aq.mzML$", colnames(group))
    
    # Create merged feature
    merged <- group[1,]
    # Sum peak areas
    merged[,intensity_cols] <- colSums(as.matrix(group[,intensity_cols]))
    # Average retention time
    merged$row.retention.time <- mean(group$row.retention.time)
    # Store row IDs of merged features
    merged$feature_ids <- paste(group$row.ID, collapse=", ")
    merged$feature_count <- nrow(group)
    
    return(merged)
  })
  
  return(do.call(rbind, merged_features))
}

# Process each bilirubin analog
ppm <- 10  # 10 ppm tolerance
merged_metabolites <- list()

for(i in 1:nrow(bilirubin_analogs)) {
  merged <- merge_features(metab_data, bilirubin_analogs$m.z[i], ppm)
  if(!is.null(merged)) {
    merged$bilirubin_analog <- bilirubin_analogs$Name[i]
    merged_metabolites[[i]] <- merged
  }
}

bilirubin_data <- do.call(rbind, merged_metabolites)

# Get common samples
common_samples <- intersect(
  grep("_aq.mzML$", colnames(bacteria_data), value = TRUE),
  grep("_aq.mzML$", colnames(bilirubin_data), value = TRUE)
)

# Convert to matrices
bacteria_mat <- as.matrix(bacteria_data[, common_samples])
bilirubin_mat <- as.matrix(bilirubin_data[, common_samples])

# Calculate correlations
correlations <- list()
for(i in 1:nrow(bacteria_mat)) {
  for(j in 1:nrow(bilirubin_mat)) {
    cor_result <- cor.test(as.numeric(bacteria_mat[i,]), 
                           as.numeric(bilirubin_mat[j,]), 
                           method = "pearson")
    
    if(!is.na(cor_result$estimate)) {
      correlations[[length(correlations) + 1]] <- data.frame(
        bacteria = bacteria_data$Name[i],
        metabolite = sprintf("m/z %.4f @ RT %.2f (features: %s)", 
                             bilirubin_data$row.m.z[j],
                             bilirubin_data$row.retention.time[j],
                             bilirubin_data$feature_ids[j]),
        bilirubin_analog = bilirubin_data$bilirubin_analog[j],
        correlation = cor_result$estimate,
        p_value = cor_result$p.value
      )
    }
  }
}

# Combine and filter results
results <- do.call(rbind, correlations)
sig_correlations <- results %>%
  filter(p_value < 0.05) %>%
  arrange(bacteria, desc(abs(correlation)))

# Save results
write.csv(sig_correlations, "bacteria_bilirubin_correlations6.csv", row.names = FALSE)

# Print correlations for each bacteria and bilirubin analog
print("Significant correlations by bacteria:")
for(bact in unique(sig_correlations$bacteria)) {
  cat("\n", bact, ":\n")
  print(sig_correlations %>% 
          filter(bacteria == bact) %>% 
          select(bilirubin_analog, metabolite, correlation, p_value) %>%
          arrange(desc(abs(correlation))) %>%
          head(5))
}


# Prepare data for network
nodes_bacteria <- data.frame(
  name = unique(sig_correlations$bacteria),
  type = "Bacteria"
)
nodes_metabolites <- data.frame(
  name = unique(sig_correlations$bilirubin_analog),
  type = "Metabolite"
)
nodes <- rbind(nodes_bacteria, nodes_metabolites)

edges <- sig_correlations %>%
  filter(correlation > 0.5 | correlation < 0) %>%  # This is the only new line
  select(from = bacteria, 
         to = bilirubin_analog, 
         correlation,
         p_value)

# Get top 10 bacteria by absolute correlation
top_10_bacteria <- edges %>%
  arrange(desc(abs(correlation))) %>%
  select(from) %>%
  distinct() %>%
  head(10) %>%
  pull(from)

# Get only the connected nodes
connected_nodes <- unique(c(edges$from, edges$to))

# Filter nodes to keep only connected ones
nodes <- nodes %>%
  filter(name %in% connected_nodes) %>%
  mutate(show_label = case_when(
    type == "Metabolite" ~ TRUE,
    name %in% top_10_bacteria ~ TRUE,
    TRUE ~ FALSE
  ))

# Create graph
graph <- tbl_graph(nodes = nodes,
                   edges = edges,
                   directed = FALSE)

# Create network plot
ggraph(graph, layout = "kk") +
  geom_edge_link(aes(edge_alpha = abs(correlation),
                     edge_width = abs(correlation),
                     color = correlation)) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(data = . %>% filter(show_label),  # Only show labels for selected nodes
                 aes(label = name), 
                 repel = TRUE, 
                 size = 4) +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_manual(values = c(Bacteria = 4, Metabolite = 5)) +
  theme_graph() +
  labs(title = "Bacteria-Bilirubin Correlation Network")
# Save plots
ggsave("correlation_network.png", width = 12, height = 8)
pdf("correlation_heatmap.png", width = 10, height = 8)

dev.off()


# Calculate adjusted p-values (q-values) using Benjamini-Hochberg
adjusted_p_values <- p.adjust(sig_correlations$p_value , method = "BH")
final_result <- sig_correlations
final_result$padj <- adjusted_p_values

edges <- final_result %>%
  filter(correlation > 0.5 | correlation < 0) %>%  # This is the only new line
  filter(padj < 0.05) %>%
  select(from = bacteria, 
         to = bilirubin_analog, 
         correlation,
         padj)

# Get top 10 bacteria by absolute correlation
top_10_bacteria <- edges %>%
  arrange(desc(abs(correlation))) %>%
  select(from) %>%
  distinct() %>%
  head(10) %>%
  pull(from)

# Get only the connected nodes
connected_nodes <- unique(c(edges$from, edges$to))

# Filter nodes to keep only connected ones
nodes <- nodes %>%
  filter(name %in% connected_nodes) %>%
  mutate(show_label = case_when(
    type == "Metabolite" ~ TRUE,
    name %in% top_10_bacteria ~ TRUE,
    TRUE ~ FALSE
  ))

# Create graph
graph <- tbl_graph(nodes = nodes,
                   edges = edges,
                   directed = FALSE)

# Create network plot
ggraph(graph, layout = "kk") +
  geom_edge_link(aes(edge_alpha = abs(correlation),
                     edge_width = abs(correlation),
                     color = correlation)) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(data = . %>% filter(show_label),  # Only show labels for selected nodes
                 aes(label = name), 
                 repel = TRUE, 
                 size = 4) +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_manual(values = c(Bacteria = 4, Metabolite = 5)) +
  theme_graph() +
  labs(title = "Bacteria-Bilirubin Correlation Network")

edges
# Get only the connected nodes
unique_correlation_bacteria <- unique(edges)

#Save results
write.csv(unique_correlation_bacteria, "unique_bacteria_bilirubin_correlations_padj.csv", row.names = FALSE)
# Save results
write.csv(final_result, "bacteria_bilirubin_correlations_padj.csv", row.names = FALSE)
# Save plots
ggsave("correlation_network.png", width = 12, height = 8)
