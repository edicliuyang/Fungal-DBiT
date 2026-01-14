suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(Matrix)
  library(viridis)
  library(OpenImageR)
  library(grid)
  library(raster)
})

## =========================
## 0) Settings
## =========================
workdir <- "~/MYH20401"
setwd(workdir)

# Start from this file (per your request)
filtered_file <- "Heart_50um.tsv"

# Optional background image for overlays
img_file <- "Heart_small.jpg"

# Grid limits (adjust if your grid differs)
xlim_grid <- c(0, 51)
ylim_grid <- c(51, 1)

# Scaling mode for colorbar limits in QC plots
scale_mode <- "median2x"  # or "quantile10_90"

# Fungal prefixes (edit to match your column naming)
prefixes <- c("AFUA_", "ANI_1_", "CAALFM_", "CAGL0", "CPAR2_", "CTRG_")

prefix_species_map <- c(
  "AFUA_"   = "A. fumigatus",
  "ANI_1_"  = "A. niger",
  "CAALFM_" = "C. albicans",
  "CAGL0"   = "C. glabrata",
  "CPAR2_"  = "C. parapsilosis",
  "CTRG_"   = "C. tropicalis"
)

prefix_colors <- c(
  "AFUA_"   = "#1f78b4",
  "ANI_1_"  = "#33a02c",
  "CAALFM_" = "#ff7f00",
  "CAGL0"   = "#6a3d9a",
  "CPAR2_"  = "#b15928",
  "CTRG_"   = "#e31a1c"
)

cluster_colors <- c(
  '0' = '#F0CE58', '1' = '#B487B7', '2' = '#289E92', '3' = '#EB545C',
  '4' = '#D7EF9B', '5' = '#EF7512', '6' = '#5084C2', '7' = '#DBA091',
  '8' = '#878787', '9' = '#FC6FCF', '10' = '#F52831', '11' = '#80FF08'
)

## =========================
## 1) Helpers
## =========================
safe_read_image <- function(path) {
  if (!file.exists(path)) return(NULL)
  OpenImageR::readImage(path)
}

calc_limits <- function(x, mode = c("median2x", "quantile10_90")) {
  mode <- match.arg(mode)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 1))
  
  if (mode == "quantile10_90") {
    lo <- as.numeric(stats::quantile(x, probs = 0.10, na.rm = TRUE))
    hi <- as.numeric(stats::quantile(x, probs = 0.90, na.rm = TRUE))
  } else {
    med <- stats::median(x, na.rm = TRUE)
    lo <- 0
    hi <- med * 2
  }
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi) || hi <= lo) hi <- max(x, na.rm = TRUE)
  c(lo, hi)
}

plot_spatial_value <- function(df_xy, value_col, title, out_pdf,
                               imported_raster = NULL,
                               limits = NULL,
                               point_size = 3.8,
                               show_axes = TRUE) {
  
  if (is.null(limits)) limits <- calc_limits(df_xy[[value_col]], scale_mode)
  
  g <- NULL
  if (!is.null(imported_raster)) {
    g <- rasterGrob(imported_raster, width = unit(1, "npc"), height = unit(1, "npc"),
                    interpolate = FALSE)
  }
  
  p <- ggplot(df_xy, aes(x = as.numeric(A), y = as.numeric(B), color = .data[[value_col]])) +
    geom_point(shape = 16, size = point_size) +
    scale_color_viridis_c(option = "inferno", limits = limits, oob = scales::squish) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
    ggtitle(title) +
    coord_equal(xlim = xlim_grid, ylim = ylim_grid) +
    scale_y_reverse(expand = expansion(mult = c(-0.013, 0.008))) +
    scale_x_continuous(expand = expansion(mult = c(-0.013, -0.013))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  if (!is.null(g)) {
    p <- p + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  }
  
  if (!show_axes) {
    p <- p + theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
  } else {
    p <- p + theme(
      axis.text  = element_text(size = 18),
      axis.title = element_text(size = 18, face = "bold")
    ) + labs(x = "X", y = "Y")
  }
  
  pdf(out_pdf, width = 8.6, height = 8.6)
  print(p)
  dev.off()
}

build_prefix_long <- function(df_wide, prefixes) {
  stopifnot("X" %in% colnames(df_wide))
  
  xy <- df_wide %>%
    tidyr::separate(X, into = c("A", "B"), sep = "x", remove = FALSE) %>%
    mutate(A = as.numeric(A), B = as.numeric(B))
  
  out <- lapply(prefixes, function(pref) {
    cols <- grep(paste0("^", pref), colnames(df_wide), value = TRUE)
    if (length(cols) == 0) return(NULL)
    tibble(
      A = xy$A,
      B = xy$B,
      Prefix = pref,
      UMI = rowSums(df_wide[, cols, drop = FALSE], na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  out
}

## =========================
## 2) Load Filtered_matrix_correct.tsv
## =========================
if (!file.exists(filtered_file)) stop("Missing file: ", filtered_file)

filtered <- read.table(filtered_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# If the first column is spot IDs but not named "X", rename it.
if (!("X" %in% colnames(filtered))) {
  message("Column 'X' not found — assuming the first column is spot IDs and renaming it to 'X'.")
  colnames(filtered)[1] <- "X"
}

# Ensure spot IDs look like "AxB"
if (!all(grepl("^[0-9]+x[0-9]+$", filtered$X))) {
  warning("Some spot IDs in column X do not match the pattern 'numberxnumber'. Check your X column formatting.")
}

# Identify gene columns (everything except X)
gene_cols <- setdiff(colnames(filtered), "X")

## =========================
## 3) QC maps (UMI + detected genes)
## =========================
imported_raster <- safe_read_image(img_file)

count <- rowSums(filtered[, gene_cols, drop = FALSE], na.rm = TRUE)
gene_count <- rowSums(filtered[, gene_cols, drop = FALSE] != 0)

qc_df <- filtered %>%
  tidyr::separate(X, into = c("A", "B"), sep = "x", remove = FALSE) %>%
  mutate(
    A = as.numeric(A),
    B = as.numeric(B),
    count = count,
    gene_count = gene_count
  )

plot_spatial_value(
  df_xy = qc_df,
  value_col = "count",
  title = "UMI",
  out_pdf = "UMI_heatmap.pdf",
  imported_raster = imported_raster,
  show_axes = TRUE
)

plot_spatial_value(
  df_xy = qc_df,
  value_col = "gene_count",
  title = "Detected genes",
  out_pdf = "Gene_heatmap.pdf",
  imported_raster = imported_raster,
  show_axes = FALSE
)

## =========================
## 4) Build Seurat object from filtered matrix
## =========================
# Seurat expects features x cells. Your file is spots x genes, so transpose.
data1 <- filtered
rownames(data1) <- data1$X
data1$X <- NULL

mat <- Matrix::Matrix(as.matrix(t(data1)), sparse = TRUE)  # genes x spots
pbmc <- CreateSeuratObject(mat, min.cells = 1, project = "fungalDBiT")

# If your mitochondrial genes are like "mt-" (mouse) or "MT-" (human), adjust pattern.
# This step is safe even if no features match.
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.8, verbose = FALSE)

pdf("Dimplot1.pdf", width = 8.6, height = 8.6)
print(DimPlot(pbmc, label = TRUE) + NoLegend())
dev.off()

## =========================
## 5) Marker genes + heatmap
## =========================
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.table(pbmc.markers, file = "markers_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top5, file = "markers_top5.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("Heatmap_top5_markers.pdf", width = 10, height = 8)
print(DoHeatmap(pbmc, features = top5$gene, group.colors = cluster_colors))
dev.off()

## =========================
## 6) Spatial cluster map (overlay)
## =========================
cluster_info <- data.frame(
  cell = colnames(pbmc),
  cluster = as.character(Idents(pbmc)),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(cell, into = c("A", "B"), sep = "x", remove = FALSE) %>%
  mutate(A = as.numeric(A), B = as.numeric(B))

# For clusters, we use discrete colors via scale_color_manual.
g <- NULL
if (!is.null(imported_raster)) {
  g <- rasterGrob(imported_raster, width = unit(1, "npc"), height = unit(1, "npc"),
                  interpolate = FALSE)
}

p_cluster <- ggplot(cluster_info, aes(x = A, y = B, color = cluster)) +
  geom_point(size = 4.2) +
  coord_equal(xlim = xlim_grid, ylim = ylim_grid) +
  scale_y_reverse(expand = expansion(mult = c(-0.013, 0.008))) +
  scale_x_continuous(expand = expansion(mult = c(-0.013, -0.013))) +
  scale_color_manual(values = cluster_colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  ggtitle("Clusters")

if (!is.null(g)) {
  p_cluster <- p_cluster + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

pdf("clusters_spatial.pdf", width = 8.6, height = 8.6)
print(p_cluster)
dev.off()

## =========================
## 7) Fungal UMI violin by cluster
## =========================
umi_long <- build_prefix_long(filtered, prefixes)

umi_long2 <- umi_long %>%
  left_join(cluster_info %>% select(A, B, cluster), by = c("A", "B")) %>%
  filter(!is.na(cluster))

# Only keep prefixes that actually exist in data
present_prefixes <- intersect(unique(umi_long2$Prefix), names(prefix_species_map))
umi_long2 <- umi_long2 %>% filter(Prefix %in% present_prefixes)

pdf("violin_plot_all_clusters_species.pdf", width = 14, height = 8)
print(
  ggplot(umi_long2, aes(x = Prefix, y = UMI, fill = Prefix)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", fill = "white") +
    scale_fill_manual(values = prefix_colors[present_prefixes]) +
    scale_x_discrete(labels = prefix_species_map[present_prefixes]) +
    facet_wrap(~ cluster, ncol = 4, scales = "fixed") +
    labs(x = NULL, y = "UMI Counts per Spot",
         title = "UMI Distribution per Fungal Species Across Clusters") +
    theme_classic(base_size = 16) +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.position = "none"
    )
)
dev.off()

message("Done. Outputs: UMI_heatmap.pdf, Gene_heatmap.pdf, Dimplot1.pdf, Heatmap_top5_markers.pdf, clusters_spatial.pdf, violin_plot_all_clusters_species.pdf")
