#### Functions-TFBS Analysis####

library(ggplot2)
library(dplyr)


filter_by_conservation <- function(df, min_conservation = 0.8, phast_track = NULL) {
  if (is.null(phast_track)) {
    stop("phast_track (GScores object) must be provided.")
  }
  
  # Create GRanges for motif hits
  hit_gr <- GRanges(
    seqnames = "chr11",  # adjust if variable
    ranges = IRanges(start = df$adjusted_start, end = df$adjusted_end),
    strand = df$strand
  )
  
  # Get conservation scores
  conservation_scores <- gscores(phast_track, hit_gr)
  
  # Add mean conservation score to df
  df$conservation_score <- mcols(conservation_scores)$default
  
  # Filter
  df_filtered <- df %>% filter(conservation_score >= min_conservation)
  
  return(df_filtered)
}

#' Filter High Confidence Transcription Factor Binding Sites
#'
#' This function filters transcription factor (TF) binding sites (TFBS) based on score thresholds,
#' minimum number of hits per TF, and optional motif clustering. It also prioritizes TFs using a weighted
#' scoring system that considers both motif strength and clustering density. Optionally, conservation
#' filtering based on phastCons scores can be applied to refine hits to evolutionarily conserved regions.
#'
#' @param df A data frame containing motif search results. Must include columns: TF_name, adjusted_start, absScore.
#' @param absScore_cutoff Numeric. Minimum absolute score threshold for filtering motifs (default = 10).
#' @param min_hits Integer. Minimum number of binding sites required per TF to be considered (default = 5).
#' @param cluster Logical. Whether to consider clustering of TF motifs based on their genomic spacing (default = TRUE).
#' @param max_avg_dist Numeric. Maximum average distance between TF binding sites to consider a cluster tight (default = 1000).
#' @param focus_tfs Character vector. Optional. A list of specific TFs to retain after filtering.
#' @param Plot Logical. Whether to generate a summary barplot of top TFs (default = TRUE).
#' @param top_n_plot Integer. Number of top TFs to plot in the barplot (default = 30).
#' @param label_top_n Integer. Number of top TFs to annotate with labels in the plot (default = 5).
#' @param sort_by Character. How to rank TFs for plotting. Options: "weighted_score", "max_absScore", "mean_absScore" (default = "weighted_score").
#' @param use_conservation Logical. Whether to filter hits based on evolutionary conservation scores (default = FALSE).
#' @param phast_track A GScores object (e.g., from phastCons35way.UCSC.mm39). Required if use_conservation = TRUE.
#' @param min_conservation Numeric. Minimum average phastCons conservation score required to retain a motif hit (default = 0.8).
#'
#' @return A filtered and prioritized data frame of TFs with calculated metrics:
#' \itemize{
#'   \item \code{avg_dist}: Average distance between binding sites
#'   \item \code{n_hits}: Number of hits per TF
#'   \item \code{mean_absScore}: Mean absolute score per TF
#'   \item \code{max_absScore}: Maximum absolute score per TF
#'   \item \code{cluster_class}: Cluster type based on spacing (Tight, Medium, Loose, Unknown)
#'   \item \code{weighted_score}: Priority score combining binding strength, clustering, and hit count
#' }
#'
#' @details
#' The function applies the following workflow:
#' \enumerate{
#'   \item Optionally filters motif hits based on conservation (if \code{use_conservation = TRUE}).
#'   \item Filters motifs above a minimum \code{absScore}.
#'   \item Retains TFs with at least \code{min_hits}.
#'   \item Calculates motif clustering (average inter-site distance).
#'   \item Assigns a weighted priority score considering motif strength, cluster tightness, and number of hits.
#'   \item Optionally visualizes the top TFs.
#' }
#'
#' This is useful for identifying likely regulatory TFs in specific genomic regions (e.g., deleted or inserted segments).
#'
#' @examples
#' \dontrun{
#' library(phastCons35way.UCSC.mm39)
#' phast <- getGScores("phastCons35way.UCSC.mm39")
#' result <- filter_high_conf_tfbs(all_hits_df,
#'                                 absScore_cutoff = 10,
#'                                 min_hits = 5,
#'                                 use_conservation = TRUE,
#'                                 phast_track = phast,
#'                                 min_conservation = 0.85)
#' }
#'
#' @export
filter_high_conf_tfbs <- function(df,
                                  absScore_cutoff = 10,
                                  min_hits = 5,
                                  cluster = TRUE,
                                  max_avg_dist = 1000,
                                  focus_tfs = NULL,
                                  Plot = TRUE,
                                  top_n_plot = 30,
                                  label_top_n = 5,
                                  sort_by = "weighted_score",
                                  fontsize = 8,
                                  use_conservation = FALSE,
                                  phast_track = NULL,
                                  min_conservation = 0.8 ) {
  
  # Optionally filter by conservation
  if (use_conservation) {
    if (is.null(phast_track)) {
      stop("To use conservation filtering, supply phast_track (e.g., phastCons35way.UCSC.mm39).")
    }
    df <- filter_by_conservation(df, min_conservation = min_conservation, phast_track = phast_track)
  }
  
  # Compute motif center if not already present
  if (!"motif_center" %in% colnames(df)) {
    df$motif_center <- floor((df$adjusted_start + df$adjusted_end) / 2)
  }
  
  
  # Score-based filtering
  df_filtered <- df %>%
    filter(absScore >= absScore_cutoff)
  
  # TFs with sufficient hits
  tf_counts <- df_filtered %>%
    group_by(TF_name) %>%
    filter(n() >= min_hits) %>%
    ungroup()
  
  # Clustering analysis
  tf_density <- tf_counts %>%
    group_by(TF_name) %>%
    summarise(
      avg_dist = ifelse(n() > 1, mean(diff(sort(motif_center))), NA),
      n_hits = n(),
      mean_absScore = mean(absScore),
      max_absScore = max(absScore),
      .groups = "drop"
    )
  
  tf_density <- tf_density %>%
    mutate(
      cluster_penalty = ifelse(
        is.na(avg_dist),
        0.8,
        pmin(1, max(0.5, 1 - (avg_dist / (2 * max_avg_dist))))
      ),
      weighted_score = (mean_absScore * cluster_penalty) * log2(n_hits + 1)
    )
  
  # Cluster classes
  tf_density <- tf_density %>%
    mutate(
      cluster_class = case_when(
        is.na(avg_dist) ~ "Unknown",
        avg_dist <= 500 ~ "Tight Cluster",
        avg_dist <= 1500 ~ "Medium Cluster",
        TRUE ~ "Loose Cluster"
      )
    ) %>% arrange(desc(weighted_score))
  
  # Optional focus TFs
  if (!is.null(focus_tfs)) {
    tf_density <- tf_density %>%
      filter(TF_name %in% focus_tfs)
  }
  
  # Plotting
  if (Plot) {
    tf_density_top <- tf_density %>%
      arrange(desc(get(sort_by))) %>%
      slice_head(n = top_n_plot)
    
    tf_list_string <- paste(tf_density_top$TF_name, collapse = ",")
    print(tf_list_string)
    
    top_labels <- tf_density_top %>%
      slice_max(weighted_score, n = label_top_n) %>%
      pull(TF_name)
    
    p <- ggplot(tf_density_top, aes(x = reorder(TF_name, weighted_score),
                                    y = weighted_score,
                                    fill = cluster_class)) +
      geom_col() +
      coord_flip() +
      geom_text(aes(label = ifelse(TF_name %in% top_labels, TF_name, "")),
                hjust = -0.1, size = fontsize/2, color = "black", fontface = "bold") +
      scale_fill_manual(values = c(
        "Tight Cluster" = "#FF5C5C",
        "Medium Cluster" = "#FFB84D",
        "Loose Cluster" = "#66B2FF",
        "Unknown" = "grey"
      )) +
      labs(title = paste0("Top Prioritized TFs (Filtered by ", sort_by, ")"),
           x = "TF Name", y = "Weighted Score", fill = "Cluster Type") +
      theme_minimal() +
      theme(axis.text = element_text(size = fontsize),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold"),
            legend.position = "bottom") +
      ylim(0, max(tf_density_top$weighted_score) * 1.15)
    
    print(p)
  }
  
  return(tf_density)
}


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

bar_plot_motif_count <- function(motif_df){
  
  ## Count motifs per TF and plot as a stacked barplot
  
  motif_count_plot <- motif_df %>%
    group_by(TF_name) %>%
    summarise(n_hits = n()) %>%
    arrange(desc(n_hits))
  
  ggplot(motif_count_plot, aes(x = reorder(TF_name, n_hits), y = n_hits, fill = TF_name)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    coord_flip() +
    labs(title = "Motif Count per TF",
         x = "Transcription Factor", y = "Number of Motif Hits") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10))
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

