#### filter_high_conf_tfbs ####
filter_high_conf_tfbs <- function(df,
                                  absScore_cutoff = 10, 
                                  min_hits = 5, 
                                  cluster = TRUE, 
                                  max_avg_dist = 1000,
                                  focus_tfs = NULL,
                                  Plot = TRUE,
                                  top_n_plot = 30,
                                  label_top_n = 5, 
                                  sort_by="weighted_score") {
  
  library(dplyr)
  library(ggplot2)
  
  # Step 1: Score-based filtering
  df_filtered <- df %>%
    filter(absScore >= absScore_cutoff)
  
  # Step 2: TFs with sufficient hits
  tf_counts <- df_filtered %>%
    group_by(TF_name) %>%
    filter(n() >= min_hits) %>%
    ungroup()
  
  # Step 3: Clustering analysis
  tf_density <- tf_counts %>%
    group_by(TF_name) %>%
    summarise(
      avg_dist = ifelse(n() > 1, mean(diff(sort(adjusted_start))), NA),
      n_hits = n(),
      mean_absScore = mean(absScore),
      max_absScore = max(absScore),
      .groups = "drop"
    )
  
  # Step 4: Priority scoring
  tf_density <- tf_density %>%
    mutate(
      cluster_penalty = ifelse(
        is.na(avg_dist), 
        0.8,  # Still penalize unknown clusters
        pmin(1, max(0.5, 1 - (avg_dist / (2 * max_avg_dist))))
      ),
      weighted_score = (mean_absScore * cluster_penalty) * log2(n_hits + 1)
    )
  
  # Step 5: Add dynamic cluster classes
  tf_density <- tf_density %>%
    mutate(
      cluster_class = case_when(
        is.na(avg_dist) ~ "Unknown",
        avg_dist <= 500 ~ "Tight Cluster",
        avg_dist <= 1500 ~ "Medium Cluster",
        TRUE ~ "Loose Cluster"
      )
    ) %>% arrange(desc(weighted_score))
  
  # Step 6: Focus on specific TFs if requested
  if (!is.null(focus_tfs)) {
    tf_density <- tf_density %>% 
      filter(TF_name %in% focus_tfs)
  }
  
  # Step 7: Optional plotting
  if (Plot) {
    if (sort_by == "weighted_score"){
      tf_density_top <- tf_density %>%
        arrange(desc(weighted_score)) %>%
        slice_head(n = top_n_plot)
    }else if (sort_by == "max_absScore"){
      tf_density_top <- tf_density %>%
        arrange(desc(max_absScore)) %>%
        slice_head(n = top_n_plot)
    }else if (sort_by == "mean_absScore"){
      tf_density_top <- tf_density %>%
        arrange(desc(mean_absScore)) %>%
        slice_head(n = top_n_plot)
    }
    
    
    # Format to paste the top Tfs into GXD (https://www.informatics.jax.org/expression.shtml)
    # Collapse into a single comma-separated string
    tf_list_string <- tf_density[1:top_n_plot,]$TF_name %>% paste(collapse = ",")
    cat(tf_list_string)
    
    # Identify top N TFs for labeling
    top_labels <- tf_density_top %>%
      slice_max(weighted_score, n = label_top_n) %>%
      pull(TF_name)
    
    p <- ggplot(tf_density_top, aes(x = reorder(TF_name, weighted_score),
                                    y = weighted_score,
                                    fill = cluster_class)) +
      geom_col() +
      coord_flip() +
      geom_text(aes(label = ifelse(TF_name %in% top_labels, TF_name, "")),
                hjust = -0.1, size = 4, color = "black", fontface = "bold") +
      scale_fill_manual(values = c(
        "Tight Cluster" = "#FF5C5C",     # Red
        "Medium Cluster" = "#FFB84D",    # Orange
        "Loose Cluster" = "#66B2FF",     # Blue
        "Unknown" = "grey"
      )) +
      labs(title = "Top Prioritized TFs (High Confidence)",
           x = "TF Name", y = "Weighted Score",
           fill = "Cluster Type") +
      theme_minimal() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold"),
            legend.position = "bottom") +
      ylim(0, max(tf_density_top$weighted_score) * 1.15)  # Add space for text
    
    print(p)
  }
  
  return(tf_density)
}


library(ggplot2)
library(dplyr)

plot_tf_counts_comparison <- function(filtered_df, clustered_df) {
  # Count TF hits
  counts_filtered <- filtered_df %>%
    count(TF_name, name = "n_hits") %>%
    mutate(type = "Filtered")
  print(counts_filtered %>% head())
  counts_clustered <- clustered_df %>%
    count(TF_name, name = "n_hits") %>%
    mutate(type = "Clustered")
  print(counts_clustered %>% head())
  combined <- bind_rows(counts_filtered, counts_clustered)
  print(combined %>% head())
  # Keep top TFs to reduce clutter
  top_tfs <- combined %>%
    group_by(TF_name) %>%
    summarise(max_hits = max(n_hits)) %>%
    top_n(20, max_hits) %>%
    pull(TF_name)
  print(top_tfs %>% head())
  combined <- combined %>% filter(TF_name %in% top_tfs)
  
  ggplot(combined, aes(x = reorder(TF_name, n_hits), y = n_hits, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = "TF Motif Hit Counts: Filtered vs Clustered",
         x = "Transcription Factor", y = "Number of Hits") +
    scale_fill_manual(values = c("Filtered" = "#4575b4", "Clustered" = "#d73027")) +
    theme_minimal()
}

plot_tf_heatmap <- function(df, bin_size = 500, return_summary = TRUE) {
  library(dplyr)
  library(ggplot2)
  library(viridis)
  
  # Step 1: Bin the data
  df_binned <- df %>%
    mutate(bin = floor(adjusted_start / bin_size) * bin_size) %>%
    count(TF_name, bin)
  
  # Step 2: Create the heatmap
  p <- ggplot(df_binned, aes(x = bin, y = TF_name, fill = n)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", name = "Hit Count") +
    labs(title = "TF Binding Site Distribution Across Region",
         x = "Genomic Bin (bp)", y = "TF Name") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  if (return_summary) {
    # Step 3: Summarize TF density per bin
    tf_hot_bins <- df_binned %>%
      group_by(bin) %>%
      summarise(
        n_TFs = n_distinct(TF_name),
        total_hits = sum(n),
        .groups = "drop"
      ) %>%
      arrange(desc(n_TFs))
    
    return(list(plot = p, summary = tf_hot_bins))
  } else {
    return(p)
  }
}

