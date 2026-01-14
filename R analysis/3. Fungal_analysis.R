suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(OpenImageR)
  library(ggforce)
  library(pheatmap)
  library(igraph)
  library(ggraph)
  library(vegan)
  library(patchwork)
})

## =========================
## 0) User settings
## =========================
dir <- "D:/data/xiaohui/Analysis/MYH20401"
setwd(dir)

matrix_file <- "fungi_expression_matrix.tsv"
img_file <- "Heart_small.jpg"

# spatial bounds for your 52x52-like grid
xlim_grid <- c(0, 51)
ylim_grid <- c(51, 1)

# Prefix list (keep exactly as your column naming)
prefixes <- c("AFUA_", "ANI_1_", "CAALFM_", "CAGL0", "C5L36_", "CPAR2_", "CTRG_")

# Named colors (important: names must match prefixes)
prefix_colors <- c(
  "AFUA_"   = "#F0CE58",
  "ANI_1_"  = "#B487B7",
  "CAALFM_" = "#EB545C",
  "CAGL0"   = "#5084C2",
  "C5L36_"  = "#289E92",
  "CPAR2_"  = "#DBA091",
  "CTRG_"   = "#D7EF9B"
)

## =========================
## 1) Helpers
## =========================
read_background <- function(path) {
  if (!file.exists(path)) return(NULL)
  OpenImageR::readImage(path)
}

add_bg <- function(p, bg_img) {
  if (is.null(bg_img)) return(p)
  g <- rasterGrob(bg_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
  p + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

# Robust upper limit:
# - your original: 10 * median(UMI_sum)
# - if median==0 but there are nonzeros -> upper becomes 0 and plot breaks visually
calc_upper_limit <- function(x, mult = 10, min_upper = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(min_upper)
  med <- median(x, na.rm = TRUE)
  upper <- mult * med
  if (!is.finite(upper) || upper <= 0) {
    upper <- max(x, na.rm = TRUE)
  }
  max(upper, min_upper)
}

# Build spot x prefix matrix in one place
build_spot_prefix_umi <- function(df, prefixes) {
  out <- df %>% dplyr::select(A, B)
  for (pref in prefixes) {
    genes <- grep(paste0("^", pref), names(df), value = TRUE)
    if (length(genes) == 0) {
      out[[pref]] <- 0
    } else {
      out[[pref]] <- rowSums(df[, genes, drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

# Long format (A,B,Prefix,UMI)
build_umi_long <- function(spot_prefix_umi, prefixes) {
  spot_prefix_umi %>%
    tidyr::pivot_longer(cols = all_of(prefixes), names_to = "Prefix", values_to = "UMI")
}

## =========================
## 2) Load data (ONLY ONCE)
## =========================
df0 <- read.table(matrix_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
df <- data.frame(X = rownames(df0), df0, check.names = FALSE) %>%
  tidyr::separate(X, into = c("A", "B"), sep = "x", convert = TRUE)

bg_img <- read_background(img_file)

# Safety: if A/B failed parsing, stop early
if (any(is.na(df$A)) || any(is.na(df$B))) {
  stop("Some spot IDs could not be parsed into A/B using 'x'. Please check rownames format like '12x5'.")
}

## =========================
## 3) Per-prefix spatial UMI heatmaps
## =========================
for (pref in prefixes) {
  genes <- grep(paste0("^", pref), names(df), value = TRUE)
  
  if (length(genes) == 0) {
    message("No genes found for prefix: ", pref, " (skip)")
    next
  }
  
  umi_sum_df <- df %>%
    dplyr::select(A, B, all_of(genes)) %>%
    mutate(UMI_sum = rowSums(across(all_of(genes)), na.rm = TRUE)) %>%
    dplyr::select(A, B, UMI_sum)
  
  if (sum(umi_sum_df$UMI_sum, na.rm = TRUE) <= 0) {
    message("Prefix ", pref, " has all-zero UMI_sum (skip)")
    next
  }
  
  lower_limit <- 0
  upper_limit <- calc_upper_limit(umi_sum_df$UMI_sum, mult = 10, min_upper = 1)
  
  p <- ggplot(umi_sum_df, aes(x = as.numeric(A), y = as.numeric(B), color = UMI_sum)) +
    geom_point(shape = 16, size = 3.5) +
    scale_color_gradient(
      low = "white",
      high = "red",
      limits = c(lower_limit, upper_limit),
      oob = scales::squish
    ) +
    ggtitle(paste0(pref, " UMI")) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    coord_equal(xlim = xlim_grid, ylim = ylim_grid) +
    scale_y_reverse(expand = expansion(mult = c(-0.013, 0.008))) +
    scale_x_continuous(expand = expansion(mult = c(-0.013, -0.013))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
  
  p <- add_bg(p, bg_img)
  
  pdf(file = paste0(pref, "_UMI_heatmap.pdf"), width = 8.6, height = 8.6)
  print(p)
  dev.off()
}

## =========================
## 4) Build long format + fractions (for pie chart, violin, etc.)
## =========================
spot_prefix_umi <- build_spot_prefix_umi(df, prefixes)

umi_long <- build_umi_long(spot_prefix_umi, prefixes)

spot_totals <- umi_long %>%
  group_by(A, B) %>%
  summarise(Total_UMI = sum(UMI, na.rm = TRUE), .groups = "drop")

umi_long <- umi_long %>%
  left_join(spot_totals, by = c("A", "B")) %>%
  mutate(Fraction = ifelse(Total_UMI > 0, UMI / Total_UMI, 0))

## =========================
## 5) Spatial pie chart (arc bars)
## =========================
umi_long_arc <- umi_long %>%
  group_by(A, B) %>%
  arrange(Prefix, .by_group = TRUE) %>%
  mutate(
    start = 2 * pi * cumsum(lag(Fraction, default = 0)),
    end   = 2 * pi * cumsum(Fraction)
  ) %>%
  ungroup()

# filter + radius scale
umi_long_arc <- umi_long_arc %>%
  filter(Total_UMI > 5) %>%                  # your threshold
  mutate(radius = sqrt(Total_UMI) * 0.01)    # your scaling

p_pie <- ggplot(umi_long_arc) +
  geom_arc_bar(aes(
    x0 = as.numeric(A),
    y0 = as.numeric(B),
    r0 = 0,
    r = radius,
    start = start,
    end = end,
    fill = Prefix
  ), color = NA, alpha = 0.6) +
  scale_fill_manual(values = prefix_colors) +
  coord_fixed(xlim = xlim_grid, ylim = ylim_grid) +
  scale_y_reverse() +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  ) +
  ggtitle("Spatial Pie Chart of Fungal Prefix UMI Fractions")

p_pie <- add_bg(p_pie, bg_img)

pdf("spatial_pie_chart.pdf", width = 8.6, height = 8.6)
print(p_pie)
dev.off()

## =========================
## 6) Violin plot (ONLY ONCE; pick one limit)
## =========================
ymax <- 300  # set once here; previously you wrote 500 then 300 in duplicate blocks

pdf("violin_plot_prefix_UMI.pdf", width = 8, height = 6)
print(
  ggplot(umi_long, aes(x = Prefix, y = UMI, fill = Prefix)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    scale_fill_manual(values = prefix_colors) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "UMI Counts per Spot",
         title = "Distribution of UMI Counts per Fungal Prefix") +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      legend.position = "none"
    )
)
dev.off()

## =========================
## 7) Top genes per prefix (single-page multi-plot 3x3)
## =========================
for (pref in prefixes) {
  genes <- grep(paste0("^", pref), names(df), value = TRUE)
  if (length(genes) == 0) {
    message("No genes found for prefix: ", pref, " (skip top genes)")
    next
  }
  
  gene_totals <- colSums(df[, genes, drop = FALSE], na.rm = TRUE)
  top9 <- names(sort(gene_totals, decreasing = TRUE))[1:min(9, length(gene_totals))]
  
  df_subset <- df %>% dplyr::select(A, B, all_of(top9))
  
  plot_list <- lapply(top9, function(gene) {
    ggplot(df_subset, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[gene]])) +
      geom_point(size = 1) +
      scale_color_gradientn(colors = c("white", "red"), oob = scales::squish) +
      ggtitle(gene) +
      coord_fixed(xlim = xlim_grid, ylim = ylim_grid) +
      scale_y_reverse() +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA)
      )
  })
  
  # Ensure exactly 9 panels for a 3x3 layout (pad with blank plots if needed)
  if (length(plot_list) < 9) {
    plot_list <- c(plot_list, rep(list(ggplot() + theme_void()), 9 - length(plot_list)))
  }
  
  combined <- wrap_plots(plot_list, ncol = 3)
  
  pdf(paste0(pref, "_top9_genes_spatial_3x3.pdf"), width = 12, height = 12)
  print(combined)
  dev.off()
}

## =========================
## 8) Prefix totals pie chart (across all spots)
## =========================
prefix_totals <- spot_prefix_umi %>%
  summarise(across(all_of(prefixes), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Prefix", values_to = "Total_UMI") %>%
  mutate(Fraction = Total_UMI / sum(Total_UMI)) %>%
  mutate(Prefix = factor(Prefix, levels = prefixes))

pdf("prefix_total_UMI_piechart.pdf", width = 6, height = 6)
print(
  ggplot(prefix_totals, aes(x = "", y = Fraction, fill = Prefix)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = prefix_colors) +
    geom_text(aes(label = paste0(round(Fraction * 100, 1), "%")),
              position = position_stack(vjust = 0.5), size = 5) +
    theme_void() +
    ggtitle("Prefix UMI Proportion Across All Spots") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    )
)
dev.off()

## =========================
## 9) Prefix co-localization correlation heatmap
## =========================
prefix_only <- spot_prefix_umi %>% select(all_of(prefixes))
cor_matrix <- cor(prefix_only, method = "pearson", use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)] <- 0

pdf("prefix_colocalization_heatmap.pdf", width = 7, height = 6)
pheatmap(
  cor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 14,
  fontsize_number = 12,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "grey80",
  main = "Prefix Co-localization Heatmap",
  angle_col = 45
)
dev.off()

## =========================
## 10) Dominant prefix map + Shannon entropy map
## =========================
dominant_df <- spot_prefix_umi %>%
  rowwise() %>%
  mutate(DominantPrefix = prefixes[which.max(c_across(all_of(prefixes)))]) %>%
  ungroup()

pdf("dominant_prefix_map.pdf", width = 8, height = 8)
print(
  ggplot(dominant_df, aes(x = as.numeric(A), y = as.numeric(B), color = DominantPrefix)) +
    geom_point(size = 2) +
    scale_color_manual(values = prefix_colors) +
    coord_fixed(xlim = xlim_grid, ylim = ylim_grid) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Dominant Prefix per Spot") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    )
)
dev.off()

entropy_df <- spot_prefix_umi %>%
  mutate(Entropy = vegan::diversity(select(., all_of(prefixes)), index = "shannon"))

pdf("spatial_entropy_map.pdf", width = 8, height = 8)
print(
  ggplot(entropy_df, aes(x = as.numeric(A), y = as.numeric(B), color = Entropy)) +
    geom_point(size = 2.5) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed(xlim = xlim_grid, ylim = ylim_grid) +
    scale_y_reverse() +
    theme_void() +
    ggtitle("Spatial Shannon Entropy (Diversity)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      legend.text = element_text(size = 14)
    )
)
dev.off()

## =========================
## 11) Prefix co-occurrence network (cor threshold)
## =========================
cor_threshold <- 0.3

cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  filter(Var1 != Var2) %>%
  filter(abs(Freq) >= cor_threshold)

if (nrow(cor_long) > 0) {
  edges <- cor_long %>% rename(from = Var1, to = Var2, weight = Freq)
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Node attributes: total abundance + degree
  prefix_total_umi <- spot_prefix_umi %>%
    summarise(across(all_of(prefixes), sum, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Prefix", values_to = "TotalUMI")
  
  V(g)$TotalUMI <- prefix_total_umi$TotalUMI[match(V(g)$name, prefix_total_umi$Prefix)]
  V(g)$Degree <- degree(g)
  
  pdf("prefix_cooccurrence_network.pdf", width = 8, height = 8)
  print(
    ggraph(g, layout = "fr") +
      geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8, show.legend = TRUE) +
      geom_node_point(aes(size = Degree, color = TotalUMI)) +
      geom_node_text(aes(label = name), repel = TRUE, size = 4) +
      scale_color_viridis_c(option = "plasma", name = "Total UMI") +
      scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Correlation") +
      theme_void() +
      ggtitle(paste0("Prefix Co-occurrence Network (|r| ≥ ", cor_threshold, ")")) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.position = "right"
      )
  )
  dev.off()
} else {
  message("No prefix pairs pass correlation threshold = ", cor_threshold, " (skip network plot).")
}

message("Done. PDFs written to: ", dir)
