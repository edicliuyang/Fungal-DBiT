################################################################################
# Spatial fungi (TOP7 prefixes only) + host analysis — single combined pipeline
# - EVERYTHING fungi-related is restricted to genes from cfg$prefixes (Top 7)
################################################################################

############################
# CONFIG (edit these)
############################
cfg <- list(
  workdir = "~/MYH20401",
  image_path = "MYH20401_small.jpg",      # optional background image
  fungi_matrix_path = "fungi_expression_matrix.tsv",
  host_matrix_path  = "non_prefix_expression_matrix.tsv",  # optional; set NULL if none
  
  # Coordinate grid bounds
  xlim = c(0, 51),
  ylim = c(51, 1),
  
  # Top 7 Prefix list (ONLY these are used anywhere below)
  prefixes <- c("AFUA_", "ANI_1_", "CAALFM_", "CAGL0", "C5L36_", "CPAR2_", "CTRG_"),
  
  # Plot tuning
  point_size = 3.5,
  scale_mult = 4,
  prefix_upper_mult = 10,
  min_total_umi_for_pie = 5,
  pie_scaling_factor = 0.04,
  
  # Violin y-limits
  violin_ylim1 = 500,
  violin_ylim2 = 300,
  
  # Correlation thresholds
  prefix_cor_threshold = 0.3,
  fungi_host_cor_threshold = 0.3,
  
  # Host gene selection (if host matrix provided)
  host_top_variable_n = 50,
  host_remove_regex = "^(Gm|LOC|mt\\.|mt-|mt_|MT-|Rn18s|Rn28s|Rn5s)",
  
  # K-means / niche clustering
  niche_k = 3,
  niche_seed = 123,
  
  # Output prefix
  out_prefix = "OUT"
)

############################
# Packages
############################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(scales)
})

# Optional packages (used if installed)
has_pkg <- function(x) requireNamespace(x, quietly = TRUE)

if (has_pkg("OpenImageR")) library(OpenImageR)
if (has_pkg("ggforce"))   library(ggforce)
if (has_pkg("pheatmap"))  library(pheatmap)
if (has_pkg("vegan"))     library(vegan)
if (has_pkg("igraph"))    library(igraph)
if (has_pkg("ggraph"))    library(ggraph)
if (has_pkg("patchwork")) library(patchwork)
if (has_pkg("ggrepel"))   library(ggrepel)
if (has_pkg("broom"))     library(broom)

############################
# Helpers
############################
msg <- function(...) cat(sprintf(...), "\n")
safe_mkdir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

read_expr_matrix <- function(path) {
  stopifnot(file.exists(path))
  read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
}

add_coords_from_rownames <- function(df) {
  df <- data.frame(Spot = rownames(df), df, check.names = FALSE)
  df %>% separate(Spot, into = c("A", "B"), sep = "x", convert = TRUE)
}

# ----------- TOP7 restriction core -----------
# Keep ONLY fungi genes whose names start with one of cfg$prefixes
filter_fungi_by_prefixes <- function(df_fungi, prefixes) {
  gene_cols <- setdiff(colnames(df_fungi), c("A", "B"))
  patt <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  keep_genes <- gene_cols[grepl(patt, gene_cols)]
  df_fungi %>% dplyr::select(A, B, all_of(keep_genes))
}

compute_total_umi <- function(df, exclude_cols = c("A", "B")) {
  gene_cols <- setdiff(colnames(df), exclude_cols)
  rowSums(df[, gene_cols, drop = FALSE], na.rm = TRUE)
}

# Safer detected genes: count >0
compute_detected_genes <- function(df, exclude_cols = c("A", "B")) {
  gene_cols <- setdiff(colnames(df), exclude_cols)
  rowSums(df[, gene_cols, drop = FALSE] > 0, na.rm = TRUE)
}

get_bg_grob <- function(image_path) {
  if (is.null(image_path) || !file.exists(image_path) || !has_pkg("OpenImageR")) return(NULL)
  img <- OpenImageR::readImage(image_path)
  rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
}

plot_spatial_value <- function(df, value, title, outfile,
                               xlim = c(0, 51), ylim = c(51, 1),
                               point_size = 3.5,
                               viridis_opt = "inferno",
                               lower_limit = 0, upper_limit = NULL,
                               show_axes = TRUE) {
  stopifnot(value %in% colnames(df))
  if (is.null(upper_limit)) {
    med_val <- median(df[[value]], na.rm = TRUE)
    upper_limit <- med_val * cfg$scale_mult
  }
  
  p <- ggplot(df, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[value]])) +
    geom_point(shape = 16, size = point_size) +
    scale_color_viridis_c(option = viridis_opt,
                          limits = c(lower_limit, upper_limit),
                          oob = scales::squish) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    ggtitle(title) +
    coord_equal(xlim = xlim, ylim = ylim) +
    scale_y_reverse() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  if (show_axes) {
    p <- p + theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16, face = "bold")
    )
  } else {
    p <- p + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  }
  
  pdf(outfile, width = 8.6, height = 8.6)
  print(p)
  dev.off()
}

sum_prefix_umi_per_spot <- function(df_fungi, prefixes) {
  out <- df_fungi %>% dplyr::select(A, B)
  for (prefix in prefixes) {
    genes <- grep(paste0("^", prefix), colnames(df_fungi), value = TRUE)
    out[[prefix]] <- if (length(genes) > 0) rowSums(df_fungi[, genes, drop = FALSE], na.rm = TRUE) else 0
  }
  out
}

build_umi_long <- function(spot_prefix_umi, prefixes) {
  spot_prefix_umi %>%
    pivot_longer(cols = all_of(prefixes), names_to = "Prefix", values_to = "UMI") %>%
    group_by(A, B) %>%
    mutate(
      Total_UMI = sum(UMI, na.rm = TRUE),
      Fraction = ifelse(Total_UMI > 0, UMI / Total_UMI, 0)
    ) %>%
    arrange(A, B, Prefix) %>%
    mutate(
      start = 2 * pi * cumsum(lag(Fraction, default = 0)),
      end   = 2 * pi * cumsum(Fraction)
    ) %>%
    ungroup()
}

default_prefix_colors <- function(n) {
  cols <- c(
    '#F0CE58', '#B487B7', '#289E92', '#EB545C', '#5084C2', '#DBA091', '#D7EF9B', '#EF7512',
    '#6A5ACD', '#20B2AA', '#FF69B4', '#87CEEB', '#FFD700', '#A52A2A', '#40E0D0', '#D2691E',
    '#9ACD32', '#00CED1', '#DC143C', '#8A2BE2', '#556B2F', '#FF8C00', '#7FFFD4', '#FF1493',
    '#00FA9A', '#CD5C5C', '#4169E1', '#F08080', '#ADFF2F', '#BA55D3', '#708090', '#8B4513'
  )
  rep(cols, length.out = n)
}

plot_prefix_heatmaps <- function(df_fungi, prefixes, outdir) {
  safe_mkdir(outdir)
  for (prefix in prefixes) {
    genes <- grep(paste0("^", prefix), colnames(df_fungi), value = TRUE)
    if (length(genes) == 0) next
    
    umi_sum <- rowSums(df_fungi[, genes, drop = FALSE], na.rm = TRUE)
    if (all(umi_sum == 0)) next
    
    dfp <- df_fungi %>% dplyr::select(A, B) %>% mutate(UMI_sum = umi_sum)
    
    nonzero <- dfp$UMI_sum[dfp$UMI_sum > 0]
    upper <- if (length(nonzero) >= 10) as.numeric(quantile(nonzero, 0.99, na.rm = TRUE)) else max(dfp$UMI_sum, na.rm = TRUE)
    upper <- max(upper * 1.5, 1)
    
    outfile <- file.path(outdir, paste0(prefix, "_UMI_heatmap_linear_lowOK.pdf"))
    
    p <- ggplot(dfp, aes(x = as.numeric(A), y = as.numeric(B), color = UMI_sum)) +
      geom_point(shape = 16, size = cfg$point_size) +
      scale_color_gradient(low = "white", high = "red",
                           limits = c(0, upper), oob = scales::squish) +
      ggtitle(paste0(prefix, " UMI")) +
      guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
      coord_equal(xlim = cfg$xlim, ylim = cfg$ylim) +
      scale_y_reverse() +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    pdf(outfile, width = 8.6, height = 8.6)
    print(p)
    dev.off()
  }
}

plot_spatial_pie <- function(umi_long, prefixes, prefix_colors, outfile) {
  if (!has_pkg("ggforce")) {
    msg("Package ggforce not installed; skip spatial pie chart (%s).", outfile)
    return(invisible(NULL))
  }
  
  umi_long_f <- umi_long %>%
    filter(Total_UMI > cfg$min_total_umi_for_pie) %>%
    mutate(radius = sqrt(Total_UMI) * cfg$pie_scaling_factor)
  
  pdf(outfile, width = 8.6, height = 8.6)
  p <- ggplot(umi_long_f) +
    ggforce::geom_arc_bar(aes(
      x0 = as.numeric(A), y0 = as.numeric(B),
      r0 = 0, r = radius,
      start = start, end = end,
      fill = Prefix
    ), color = NA, alpha = 0.6) +
    scale_fill_manual(values = prefix_colors) +
    coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Spatial Pie Chart of Fungal Prefix UMI Fractions (Top7)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  print(p)
  dev.off()
}

plot_violin <- function(umi_long, prefix_colors, outfile, ylim_max) {
  pdf(outfile, width = 8, height = 6)
  p <- ggplot(umi_long, aes(x = Prefix, y = UMI, fill = Prefix)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    scale_fill_manual(values = prefix_colors) +
    scale_y_continuous(limits = c(0, ylim_max),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "UMI Counts per Spot",
         title = "Distribution of UMI Counts per Fungal Prefix (Top7)") +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.position = "none"
    )
  print(p)
  dev.off()
}

plot_top_genes_spatial <- function(df_fungi, prefixes, outdir, top_n = 10) {
  safe_mkdir(outdir)
  for (prefix in prefixes) {
    genes <- grep(paste0("^", prefix), colnames(df_fungi), value = TRUE)
    if (length(genes) == 0) next
    
    gene_totals <- colSums(df_fungi[, genes, drop = FALSE], na.rm = TRUE)
    top_genes <- names(sort(gene_totals, decreasing = TRUE))[1:min(top_n, length(gene_totals))]
    
    df_subset <- df_fungi %>% dplyr::select(A, B, all_of(top_genes))
    
    pdf(file.path(outdir, paste0(prefix, "_top", top_n, "_genes_spatial.pdf")),
        width = 8, height = 7)
    for (gene in top_genes) {
      p <- ggplot(df_subset, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[gene]])) +
        geom_point(size = 2.2) +
        scale_color_gradientn(colors = c("white", "red"), oob = scales::squish) +
        ggtitle(gene) +
        coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
        scale_y_reverse() +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)
        )
      print(p)
    }
    dev.off()
    
    if (top_n >= 9 && has_pkg("patchwork")) {
      top9 <- top_genes[1:9]
      plot_list <- lapply(top9, function(gene) {
        ggplot(df_subset, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[gene]])) +
          geom_point(size = 1) +
          scale_color_gradientn(colors = c("white", "red"), oob = scales::squish) +
          ggtitle(gene) +
          coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
          scale_y_reverse() +
          theme_void() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "none",
            panel.background = element_rect(fill = "white", color = NA)
          )
      })
      combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
        (plot_list[[4]] | plot_list[[5]] | plot_list[[6]]) /
        (plot_list[[7]] | plot_list[[8]] | plot_list[[9]])
      
      pdf(file.path(outdir, paste0(prefix, "_top9_genes_spatial_3x3.pdf")),
          width = 12, height = 12)
      print(combined)
      dev.off()
    }
  }
}

plot_prefix_totals_pie <- function(spot_prefix_umi, prefixes, prefix_colors, outfile) {
  totals <- spot_prefix_umi %>%
    summarise(across(all_of(prefixes), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Prefix", values_to = "Total_UMI") %>%
    mutate(Fraction = Total_UMI / sum(Total_UMI, na.rm = TRUE)) %>%
    arrange(desc(Fraction)) %>%
    mutate(Prefix = factor(Prefix, levels = Prefix))
  
  pdf(outfile, width = 6, height = 6)
  p <- ggplot(totals, aes(x = "", y = Fraction, fill = Prefix)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = prefix_colors) +
    theme_void() +
    ggtitle("Prefix UMI Proportion Across All Spots (Top7)") +
    geom_text(aes(label = paste0(round(Fraction * 100, 1), "%")),
              position = position_stack(vjust = 0.5), size = 4) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  print(p)
  dev.off()
}

plot_prefix_correlation_heatmap <- function(spot_prefix_umi, prefixes, outfile) {
  if (!has_pkg("pheatmap")) {
    msg("Package pheatmap not installed; skip %s", outfile)
    return(invisible(NULL))
  }
  
  mat <- spot_prefix_umi %>% dplyr::select(all_of(prefixes)) %>% as.matrix()
  cor_m <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  cor_m[is.na(cor_m)] <- 0
  
  pdf(outfile, width = 10, height = 10)
  pheatmap::pheatmap(
    cor_m,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize = 12,
    fontsize_number = 9,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    border_color = NA,
    main = "Prefix Co-localization (Pearson r) (Top7)"
  )
  dev.off()
}

plot_prefix_network <- function(spot_prefix_umi, prefixes, outfile, cor_threshold = 0.3) {
  if (!has_pkg("igraph") || !has_pkg("ggraph")) {
    msg("Need igraph + ggraph; skip %s", outfile)
    return(invisible(NULL))
  }
  
  mat <- spot_prefix_umi %>% dplyr::select(all_of(prefixes)) %>% as.matrix()
  cor_m <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  cor_m[is.na(cor_m)] <- 0
  
  edges <- as.data.frame(as.table(cor_m)) %>%
    filter(Var1 != Var2) %>%
    filter(abs(Freq) >= cor_threshold) %>%
    rename(from = Var1, to = Var2, weight = Freq)
  
  if (nrow(edges) == 0) {
    msg("No edges above threshold %.2f; skip %s", cor_threshold, outfile)
    return(invisible(NULL))
  }
  
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  total_umi <- spot_prefix_umi %>%
    summarise(across(all_of(prefixes), ~ sum(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Prefix", values_to = "TotalUMI")
  
  igraph::V(g)$TotalUMI <- total_umi$TotalUMI[match(igraph::V(g)$name, total_umi$Prefix)]
  igraph::V(g)$Degree <- as.numeric(igraph::degree(g))
  
  pdf(outfile, width = 8, height = 8)
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8) +
    ggraph::geom_node_point(aes(size = Degree, color = TotalUMI)) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_color_viridis_c(option = "plasma", name = "Total UMI") +
    scale_edge_color_gradient2(low = "blue", mid = "white", high = "red",
                               midpoint = 0, name = "Correlation") +
    theme_void() +
    ggtitle("Prefix Co-occurrence Network (Top7)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right"
    )
  print(p)
  dev.off()
}

plot_dominant_prefix <- function(spot_prefix_umi, prefixes, prefix_colors, outfile) {
  mat <- spot_prefix_umi %>% dplyr::select(all_of(prefixes)) %>% as.matrix()
  dom_idx <- apply(mat, 1, which.max)
  dom <- prefixes[dom_idx]
  
  dfp <- spot_prefix_umi %>% mutate(DominantPrefix = dom)
  
  pdf(outfile, width = 8, height = 8)
  p <- ggplot(dfp, aes(x = as.numeric(A), y = as.numeric(B), color = DominantPrefix)) +
    geom_point(size = 2) +
    scale_color_manual(values = prefix_colors) +
    coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Dominant Prefix per Spot (Top7)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  print(p)
  dev.off()
}

plot_entropy <- function(spot_prefix_umi, prefixes, outfile) {
  if (!has_pkg("vegan")) {
    msg("Package vegan not installed; skip %s", outfile)
    return(invisible(NULL))
  }
  mat <- spot_prefix_umi %>% dplyr::select(all_of(prefixes))
  ent <- vegan::diversity(mat, index = "shannon")
  
  dfp <- spot_prefix_umi %>% mutate(Entropy = ent)
  
  pdf(outfile, width = 8, height = 8)
  p <- ggplot(dfp, aes(x = as.numeric(A), y = as.numeric(B), color = Entropy)) +
    geom_point(size = 2.5) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Spatial Shannon Entropy (Top7)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.text = element_text(size = 10)
    )
  print(p)
  dev.off()
}

run_niche_clustering <- function(spot_prefix_umi, prefixes, k, seed, outfile) {
  mat <- spot_prefix_umi %>% dplyr::select(all_of(prefixes)) %>% as.matrix()
  total <- rowSums(mat, na.rm = TRUE)
  frac <- ifelse(total > 0, mat / total, 0)
  
  pca <- prcomp(frac, scale. = TRUE)
  set.seed(seed)
  km <- kmeans(pca$x[, 1:min(3, ncol(pca$x)), drop = FALSE], centers = k)
  
  dfc <- spot_prefix_umi %>% mutate(NicheCluster = factor(km$cluster))
  
  pdf(outfile, width = 8, height = 8)
  p <- ggplot(dfc, aes(x = as.numeric(A), y = as.numeric(B), color = NicheCluster)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("red", "blue", "green", "purple", "orange")[1:k]) +
    coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Spatial Niches Based on Prefix Composition (Top7)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
  print(p)
  dev.off()
  
  dfc
}

marker_heatmap_by_cluster <- function(df_fungi, cluster_df, prefixes, outfile) {
  if (!has_pkg("pheatmap")) {
    msg("Package pheatmap not installed; skip %s", outfile)
    return(invisible(NULL))
  }
  
  df_full <- df_fungi %>%
    left_join(cluster_df %>% dplyr::select(A, B, NicheCluster), by = c("A", "B"))
  
  gene_cols <- setdiff(colnames(df_full), c("A", "B", "NicheCluster"))
  clusters <- levels(df_full$NicheCluster)
  marker_list <- list()
  
  for (cl in clusters) {
    in_cl  <- df_full$NicheCluster == cl
    out_cl <- df_full$NicheCluster != cl
    
    pvals <- sapply(gene_cols, function(g) suppressWarnings(wilcox.test(df_full[in_cl, g], df_full[out_cl, g])$p.value))
    p_adj <- p.adjust(pvals, method = "fdr")
    
    logFC <- sapply(gene_cols, function(g) {
      log2(mean(df_full[in_cl, g] + 1, na.rm = TRUE) / mean(df_full[out_cl, g] + 1, na.rm = TRUE))
    })
    
    tab <- data.frame(Gene = gene_cols, adj_p = p_adj, logFC = logFC) %>%
      filter(is.finite(logFC), adj_p < 1, logFC > 0) %>%
      arrange(desc(logFC))
    
    marker_list[[cl]] <- tab
  }
  
  top_markers <- unique(unlist(lapply(marker_list, function(x) head(x$Gene, 20))))
  if (length(top_markers) == 0) {
    msg("No marker genes found; skip %s", outfile)
    return(invisible(NULL))
  }
  
  avg <- df_full %>%
    group_by(NicheCluster) %>%
    summarise(across(all_of(top_markers), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()
  
  rownames(avg) <- avg$NicheCluster
  avg$NicheCluster <- NULL
  
  mat <- t(avg)
  mat_z <- t(scale(mat))
  
  pdf(outfile, width = 9, height = 11)
  pheatmap::pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red", "darkred"))(200),
    breaks = seq(-3, 3, length.out = 201),
    fontsize_row = 8,
    fontsize_col = 12,
    main = "Top Marker Genes (Top7-only fungi genes; Z-scored across clusters)",
    border_color = NA
  )
  dev.off()
}

# ---- Fungi–host cross correlation block ----
run_fungi_host_module <- function(df_fungi, df_host, prefixes, outdir) {
  safe_mkdir(outdir)
  if (is.null(df_host)) {
    msg("No host matrix provided; skip fungi–host analyses.")
    return(invisible(NULL))
  }
  
  host_var <- df_host %>%
    dplyr::select(-A, -B) %>%
    summarise(across(everything(), ~ var(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Gene", values_to = "Variance") %>%
    arrange(desc(Variance))
  
  top_host <- host_var %>% slice_head(n = cfg$host_top_variable_n) %>% pull(Gene)
  host_clean <- top_host[!grepl(cfg$host_remove_regex, top_host)]
  host_clean <- intersect(host_clean, colnames(df_host))
  
  if (length(host_clean) < 3) {
    msg("Too few host genes after filtering; skip fungi–host module.")
    return(invisible(NULL))
  }
  
  spot_prefix_umi <- sum_prefix_umi_per_spot(df_fungi, prefixes)
  
  combined <- spot_prefix_umi %>%
    inner_join(df_host %>% dplyr::select(A, B, all_of(host_clean)), by = c("A", "B"))
  
  fungal_columns <- prefixes
  host_columns <- host_clean
  
  cross <- matrix(NA_real_, nrow = length(fungal_columns), ncol = length(host_columns),
                  dimnames = list(fungal_columns, host_columns))
  
  for (i in seq_along(fungal_columns)) {
    for (j in seq_along(host_columns)) {
      x <- combined[[fungal_columns[i]]]
      y <- combined[[host_columns[j]]]
      cross[i, j] <- if (all(x == 0) || all(y == 0)) NA_real_
      else suppressWarnings(cor(x, y, method = "pearson", use = "pairwise.complete.obs"))
    }
  }
  cross[is.na(cross)] <- 0
  
  if (has_pkg("pheatmap")) {
    pdf(file.path(outdir, "fungi_host_cross_correlation_heatmap.pdf"), width = 12, height = 10)
    pheatmap::pheatmap(
      cross,
      cluster_rows = TRUE, cluster_cols = TRUE,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      main = "Fungal Prefix (Top7) × Host Gene Cross-Correlation",
      fontsize = 12, fontsize_row = 10, fontsize_col = 9,
      border_color = NA, angle_col = 45
    )
    dev.off()
  }
  
  top_pairs <- as.data.frame(as.table(cross)) %>%
    rename(FungalPrefix = Var1, HostGene = Var2, Correlation = Freq) %>%
    mutate(AbsCorrelation = abs(Correlation)) %>%
    arrange(desc(AbsCorrelation)) %>%
    slice_head(n = 30)
  
  write.table(top_pairs, file.path(outdir, "top30_fungi_host_correlation_pairs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  combined <- combined %>%
    mutate(FungalTotalUMI = rowSums(dplyr::select(., all_of(prefixes)), na.rm = TRUE))
  med <- median(combined$FungalTotalUMI, na.rm = TRUE)
  combined <- combined %>% mutate(FungalGroup = ifelse(FungalTotalUMI > med, "HighFungi", "LowFungi"))
  
  de <- lapply(host_columns, function(g) {
    tt <- suppressWarnings(t.test(combined[[g]] ~ combined$FungalGroup))
    data.frame(
      Gene = g,
      logFC = mean(combined[[g]][combined$FungalGroup == "HighFungi"], na.rm = TRUE) -
        mean(combined[[g]][combined$FungalGroup == "LowFungi"], na.rm = TRUE),
      pval = tt$p.value
    )
  }) %>% bind_rows() %>%
    mutate(adj_pval = p.adjust(pval, method = "fdr"))
  
  if (has_pkg("ggrepel")) {
    pdf(file.path(outdir, "volcano_hostgenes_fungal_burden.pdf"), width = 8, height = 7)
    p <- ggplot(de, aes(x = logFC, y = -log10(adj_pval))) +
      geom_point(aes(color = adj_pval < 0.05), size = 2) +
      scale_color_manual(values = c("grey60", "red"), guide = "none") +
      ggrepel::geom_text_repel(aes(label = ifelse(adj_pval < 0.05, Gene, NA)),
                               size = 3, max.overlaps = 12) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_classic(base_size = 16) +
      labs(title = "Host gene response: High vs Low fungal burden (Top7)",
           x = "Mean difference (High - Low)", y = "-log10(FDR)")
    print(p)
    dev.off()
  }
  
  if (has_pkg("igraph") && has_pkg("ggraph")) {
    cor_long <- as.data.frame(as.table(cross)) %>%
      rename(Fungi = Var1, HostGene = Var2, Correlation = Freq) %>%
      filter(abs(Correlation) >= cfg$fungi_host_cor_threshold)
    
    if (nrow(cor_long) > 0) {
      g <- igraph::graph_from_data_frame(cor_long, directed = FALSE)
      V(g)$Type <- ifelse(V(g)$name %in% fungal_columns, "Fungi", "Host")
      
      pdf(file.path(outdir, "fungi_host_correlation_network.pdf"), width = 10, height = 8)
      p <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(aes(width = abs(Correlation), color = Correlation), alpha = 0.8) +
        ggraph::geom_node_point(aes(color = Type), size = 5) +
        ggraph::geom_node_text(aes(label = name, color = Type), repel = TRUE, size = 3) +
        scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        scale_color_manual(values = c("Fungi" = "orange", "Host" = "steelblue")) +
        theme_void() +
        labs(title = paste0("Fungi–Host correlation network (Top7; |r| ≥ ", cfg$fungi_host_cor_threshold, ")")) +
        theme(legend.position = "right")
      print(p)
      dev.off()
    }
  }
  
  invisible(list(cross = cross, combined = combined, host_genes = host_columns))
}

############################
# RUN
############################
setwd(cfg$workdir)
safe_mkdir(file.path(cfg$workdir, "results"))

# Load fungi matrix (full), then restrict to TOP7 genes only
msg("Loading fungi matrix: %s", cfg$fungi_matrix_path)
df_fungi_raw <- read_expr_matrix(cfg$fungi_matrix_path)
df_fungi_all <- add_coords_from_rownames(df_fungi_raw)

# Restrict fungi genes strictly to cfg$prefixes
df_fungi <- filter_fungi_by_prefixes(df_fungi_all, cfg$prefixes)
msg("Fungi genes kept after Top7 prefix filtering: %d", ncol(df_fungi) - 2)

# Load host matrix (optional; not prefix-restricted)
df_host <- NULL
if (!is.null(cfg$host_matrix_path) && file.exists(cfg$host_matrix_path)) {
  msg("Loading host matrix: %s", cfg$host_matrix_path)
  df_host_raw <- read_expr_matrix(cfg$host_matrix_path)
  df_host <- add_coords_from_rownames(df_host_raw)
} else {
  msg("No host matrix found (or not provided). Host analyses will be skipped.")
}

# Background image (optional)
bg_grob <- get_bg_grob(cfg$image_path)

# Basic UMI & gene counts (TOP7-only!)
df_fungi$UMI_total  <- compute_total_umi(df_fungi, exclude_cols = c("A", "B"))
df_fungi$Gene_count <- compute_detected_genes(df_fungi, exclude_cols = c("A", "B"))

umi_upper  <- median(df_fungi$UMI_total,  na.rm = TRUE) * cfg$scale_mult
gene_upper <- median(df_fungi$Gene_count, na.rm = TRUE) * cfg$scale_mult

plot_spatial_value(
  df_fungi, "UMI_total", "UMI (Top7 prefixes only)",
  outfile = file.path("results", paste0(cfg$out_prefix, "_UMI_heatmap_top7.pdf")),
  xlim = cfg$xlim, ylim = cfg$ylim,
  point_size = cfg$point_size,
  lower_limit = 0, upper_limit = umi_upper,
  show_axes = TRUE
)

plot_spatial_value(
  df_fungi, "Gene_count", "Detected genes (Top7 prefixes only)",
  outfile = file.path("results", paste0(cfg$out_prefix, "_Gene_heatmap_top7.pdf")),
  xlim = cfg$xlim, ylim = cfg$ylim,
  point_size = cfg$point_size + 0.5,
  lower_limit = 0, upper_limit = gene_upper,
  show_axes = FALSE
)

cfg$prefixes = c("AFUA_", "ANI_1_", "CAALFM_", "CAGL0", "C5L36_", "CPAR2_", "CTRG_")

# Prefix sums (TOP7 only by definition)
spot_prefix_umi <- sum_prefix_umi_per_spot(df_fungi, cfg$prefixes)

# Per-prefix UMI heatmaps
plot_prefix_heatmaps(df_fungi, cfg$prefixes, outdir = file.path("results", "prefix_heatmaps_top7"))

# Long format + pie + violin
umi_long <- build_umi_long(spot_prefix_umi, cfg$prefixes)
prefix_colors <- setNames(default_prefix_colors(length(cfg$prefixes)), cfg$prefixes)

plot_spatial_pie(
  umi_long, cfg$prefixes, prefix_colors,
  outfile = file.path("results", paste0(cfg$out_prefix, "_spatial_pie_chart_top7.pdf"))
)

plot_violin(
  umi_long, prefix_colors,
  outfile = file.path("results", paste0(cfg$out_prefix, "_violin_ylim_", cfg$violin_ylim1, "_top7.pdf")),
  ylim_max = cfg$violin_ylim1
)
plot_violin(
  umi_long, prefix_colors,
  outfile = file.path("results", paste0(cfg$out_prefix, "_violin_ylim_", cfg$violin_ylim2, "_top7.pdf")),
  ylim_max = cfg$violin_ylim2
)

# Top genes per prefix (genes are Top7-only)
plot_top_genes_spatial(df_fungi, cfg$prefixes, outdir = file.path("results", "top_genes_per_prefix_top7"), top_n = 3)

# Prefix totals pie
plot_prefix_totals_pie(
  spot_prefix_umi, cfg$prefixes, prefix_colors,
  outfile = file.path("results", paste0(cfg$out_prefix, "_prefix_total_UMI_piechart_top7.pdf"))
)

# Prefix co-localization
plot_prefix_correlation_heatmap(
  spot_prefix_umi, cfg$prefixes,
  outfile = file.path("results", paste0(cfg$out_prefix, "_prefix_colocalization_heatmap_top7.pdf"))
)

plot_prefix_network(
  spot_prefix_umi, cfg$prefixes,
  outfile = file.path("results", paste0(cfg$out_prefix, "_prefix_cooccurrence_network_top7.pdf")),
  cor_threshold = cfg$prefix_cor_threshold
)

# Dominant prefix + entropy
plot_dominant_prefix(
  spot_prefix_umi, cfg$prefixes, prefix_colors,
  outfile = file.path("results", paste0(cfg$out_prefix, "_dominant_prefix_map_top7.pdf"))
)

plot_entropy(
  spot_prefix_umi, cfg$prefixes,
  outfile = file.path("results", paste0(cfg$out_prefix, "_spatial_entropy_map_top7.pdf"))
)

# Niche clustering + marker heatmap (marker genes are Top7-only)
cluster_df <- run_niche_clustering(
  spot_prefix_umi, cfg$prefixes,
  k = cfg$niche_k, seed = cfg$niche_seed,
  outfile = file.path("results", paste0(cfg$out_prefix, "_spatial_niche_map_prefixes_top7.pdf"))
)

marker_heatmap_by_cluster(
  df_fungi, cluster_df, cfg$prefixes,
  outfile = file.path("results", paste0(cfg$out_prefix, "_niche_marker_genes_heatmap_top7.pdf"))
)

# Save cluster assignment
write.csv(cluster_df, file.path("results", paste0(cfg$out_prefix, "_spotwise_prefix_niches_top7.csv")), row.names = FALSE)

# Fungi–host module (fungi side is Top7-only by using df_fungi + cfg$prefixes)
fh <- run_fungi_host_module(
  df_fungi = df_fungi,
  df_host  = df_host,
  prefixes = cfg$prefixes,
  outdir   = file.path("results", "fungi_host_top7")
)

# Prefix totals table (Top7)
prefix_totals_tbl <- spot_prefix_umi %>%
  summarise(across(all_of(cfg$prefixes), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Prefix", values_to = "Total_UMI") %>%
  arrange(desc(Total_UMI)) %>%
  mutate(
    Fraction = Total_UMI / sum(Total_UMI, na.rm = TRUE),
    Percent  = 100 * Fraction
  )

write.table(
  prefix_totals_tbl,
  file = file.path("results", paste0(cfg$out_prefix, "_prefix_total_UMI_table_top7.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  prefix_totals_tbl,
  file = file.path("results", paste0(cfg$out_prefix, "_prefix_total_UMI_table_top7.csv")),
  row.names = FALSE
)

print(prefix_totals_tbl)


prefix_totals_tbl <- prefix_totals_tbl %>%
  mutate(Prefix = factor(Prefix, levels = Prefix[order(Total_UMI, decreasing = TRUE)]))

p_bar <- ggplot(prefix_totals_tbl, aes(x = Prefix, y = Total_UMI)) +
  geom_col() +
  geom_text(aes(label = Total_UMI), vjust = -0.3, size = 4) +
  labs(
    title = "Total UMIs per species (prefix)",
    x = "Species (prefix)",
    y = "Total UMI"
  ) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(file.path("results", paste0(cfg$out_prefix, "_UMI_per_species_barplot.pdf")),
       plot = p_bar, width = 8, height = 8)
ggsave(file.path("results", paste0(cfg$out_prefix, "_UMI_per_species_barplot.png")),
       plot = p_bar, width = 8, height = 8, dpi = 300)

print(p_bar)



plot_top_genes_onepdf_7x3 <- function(df_fungi, prefixes, outdir,
                                      top_n = 3,
                                      outfile = "ALLprefix_topgenes_7x3.pdf") {
  safe_mkdir(outdir)
  
  if (!has_pkg("patchwork")) {
    stop("This function requires the 'patchwork' package. Please install.packages('patchwork')")
  }
  
  # Build a list-of-lists: one row per prefix, each containing 3 gene plots
  row_plots <- list()
  
  for (prefix in prefixes) {
    genes <- grep(paste0("^", prefix), colnames(df_fungi), value = TRUE)
    if (length(genes) == 0) {
      # If no genes, fill row with blank plots
      row_plots[[prefix]] <- list(
        ggplot() + theme_void() + ggtitle(paste0(prefix, ": no genes")),
        ggplot() + theme_void(),
        ggplot() + theme_void()
      )
      next
    }
    
    gene_totals <- colSums(df_fungi[, genes, drop = FALSE], na.rm = TRUE)
    top_genes <- names(sort(gene_totals, decreasing = TRUE))[1:min(top_n, length(gene_totals))]
    
    # Ensure exactly 3 slots per prefix (pad with NA if fewer)
    if (length(top_genes) < top_n) top_genes <- c(top_genes, rep(NA, top_n - length(top_genes)))
    
    plots_this_prefix <- lapply(top_genes, function(gene) {
      if (is.na(gene)) return(ggplot() + theme_void())
      
      ggplot(df_fungi, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[gene]])) +
        geom_point(size = 1.2) +
        scale_color_gradientn(colors = c("white", "red"), oob = scales::squish) +
        coord_fixed(xlim = cfg$xlim, ylim = cfg$ylim) +
        scale_y_reverse() +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA)
        ) +
        ggtitle(paste0(prefix, " | ", gene))
    })
    
    row_plots[[prefix]] <- plots_this_prefix
  }
  
  # Assemble as 7 rows × 3 columns, keeping prefix order
  combined <- NULL
  for (prefix in prefixes) {
    row_combo <- row_plots[[prefix]][[1]] | row_plots[[prefix]][[2]] | row_plots[[prefix]][[3]]
    combined <- if (is.null(combined)) row_combo else (combined / row_combo)
  }
  
  outpath <- file.path(outdir, outfile)
  pdf(outpath, width = 12, height = 24)  # tall page: 7 rows
  print(combined)
  dev.off()
  
  msg("Saved combined 7x3 top-genes PDF: %s", outpath)
}


plot_top_genes_onepdf_7x3(
  df_fungi = df_fungi,
  prefixes = cfg$prefixes,
  outdir   = file.path("results", "top_genes_per_prefix_top7"),
  top_n    = 3,
  outfile  = paste0(cfg$out_prefix, "_top3genes_ALLprefix_7rowsx3cols.pdf")
)

msg("DONE. Outputs are in: %s", file.path(cfg$workdir, "results"))






library(patchwork)

# Build cluster-level summary table
cluster_fungi_tbl <- spot_prefix_umi %>%
  left_join(cluster_df %>% dplyr::select(A, B, NicheCluster), by = c("A", "B")) %>%
  mutate(
    Fungal_UMI = rowSums(dplyr::select(., all_of(cfg$prefixes)), na.rm = TRUE)
  ) %>%
  group_by(NicheCluster) %>%
  summarise(
    Total_Fungal_UMI = sum(Fungal_UMI, na.rm = TRUE),
    Mean_Fungal_UMI_per_spot = mean(Fungal_UMI, na.rm = TRUE),
    N_spots = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Fungal_UMI))

# Ensure cluster order is consistent
cluster_fungi_tbl$NicheCluster <- factor(
  cluster_fungi_tbl$NicheCluster,
  levels = cluster_fungi_tbl$NicheCluster
)

# ---- Plot 1: total fungal UMIs per cluster ----
p_total <- ggplot(cluster_fungi_tbl,
                  aes(x = NicheCluster, y = Total_Fungal_UMI)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = Total_Fungal_UMI),
            vjust = -0.4, size = 4) +
  labs(
    title = "Total fungal UMIs per spatial niche cluster",
    x = "Niche cluster",
    y = "Total fungal UMIs (Top 7 species)"
  ) +
  theme_classic(base_size = 16)

# ---- Plot 2: mean fungal UMIs per spot ----
p_mean <- ggplot(cluster_fungi_tbl,
                 aes(x = NicheCluster, y = Mean_Fungal_UMI_per_spot)) +
  geom_col(fill = "darkorange") +
  geom_text(aes(label = round(Mean_Fungal_UMI_per_spot, 1)),
            vjust = -0.4, size = 4) +
  labs(
    title = "Mean fungal UMIs per spot",
    x = "Niche cluster",
    y = "Mean fungal UMIs per spot"
  ) +
  theme_classic(base_size = 16)

# ---- Combine into ONE figure ----
combined_plot <- p_total | p_mean

# ---- Save to ONE PDF ----
outfile <- file.path(
  "results",
  paste0(cfg$out_prefix, "_fungal_UMI_per_cluster_barplots.pdf")
)

ggsave(
  filename = outfile,
  plot = combined_plot,
  width = 12,
  height = 5
)

print(combined_plot)

msg("Saved combined cluster bar plots to: %s", outfile)
