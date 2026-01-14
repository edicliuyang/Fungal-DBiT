## =========================
## Fungal-only analysis
## =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(Matrix)
  library(OpenImageR)
  library(grid)
  library(viridis)   # for scale_color_viridis_c
})

dir <- "D:/data/xiaohui/Analysis/MYH20401"
setwd(dir)

img_file <- "Heart_small.jpg"
filtered_file <- "Filtered_matrix_correct.tsv"

# fungal prefixes
prefixes <- c("AFUA_", "ANI_1_", "CAALFM_", "CAGL0", "C5L36_", "CPAR2_", "CTRG_")

# spatial grid
xlim_grid <- c(0, 51)
ylim_grid <- c(51, 1)

# cluster colors (optional; safe even if you have fewer clusters)
cluster_colors <- c(
  '0' = '#F0CE58','1' = '#B487B7','2' = '#289E92','3' = '#EB545C','4' = '#D7EF9B',
  '5' = '#EF7512','6' = '#5084C2','7' = '#DBA091','8' = '#878787','9' = '#FC6FCF',
  '10' = '#F52831','11' = '#80FF08','12' = '#CC66FF','13' = '#FFFF0A','14' = '#FF9900'
)

## -------------------------
## helper: background raster
## -------------------------
imported_raster <- NULL
if (file.exists(img_file)) {
  imported_raster <- OpenImageR::readImage(img_file)
}
add_bg <- function(p) {
  if (is.null(imported_raster)) return(p)
  g <- rasterGrob(imported_raster, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
  p + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

## =========================
## 1) Load filtered matrix
## =========================
expr <- read.table(filtered_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Ensure first column is spot id X (e.g., "12x5")
if (!("X" %in% colnames(expr))) {
  colnames(expr)[1] <- "X"
}

## =========================
## 2) Extract fungal genes by prefixes
## =========================
# collect ALL fungal gene columns
fungal_cols <- unique(unlist(lapply(prefixes, function(pref) {
  grep(paste0("^", pref), colnames(expr), value = TRUE)
})))

if (length(fungal_cols) == 0) {
  stop("No fungal genes found with the provided prefixes. Check prefixes vs column names.")
}

fungi_expression_matrix <- expr %>% dplyr::select(X, all_of(fungal_cols))
write.table(fungi_expression_matrix, file = "fungi_expression_matrix.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

## (Optional) save non-fungal genes too (handy for host-only later)
non_fungal_cols <- setdiff(colnames(expr), c("X", fungal_cols))
non_prefix_expression_matrix <- expr %>% dplyr::select(X, all_of(non_fungal_cols))
write.table(non_prefix_expression_matrix, file = "non_prefix_expression_matrix.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Saved fungi_expression_matrix.tsv with ", length(fungal_cols), " fungal columns.")

## =========================
## 3) Fungal QC maps (UMI + detected genes)
## =========================
# UMI per spot = rowSum across fungal genes
fungi_counts <- rowSums(fungi_expression_matrix[, fungal_cols, drop = FALSE], na.rm = TRUE)
fungi_genes_detected <- rowSums(fungi_expression_matrix[, fungal_cols, drop = FALSE] != 0)

qc_df <- fungi_expression_matrix %>%
  tidyr::separate(X, into = c("A", "B"), sep = "x", convert = TRUE) %>%
  mutate(count = fungi_counts,
         gene_count = fungi_genes_detected)

# robust upper limit: if median=0 use max
upper_from_median <- function(x, mult = 8) {
  med <- median(x, na.rm = TRUE)
  up <- med * mult
  if (!is.finite(up) || up <= 0) up <- max(x, na.rm = TRUE)
  if (!is.finite(up) || up <= 0) up <- 1
  up
}

# ---- Fungi UMI heatmap ----
pdf("Fungi_UMI_heatmap.pdf", width = 8.6, height = 8.6)
p1 <- ggplot(qc_df, aes(x = as.numeric(A), y = as.numeric(B), color = count)) +
  geom_point(shape = 16, size = 4) +
  scale_color_viridis_c(option = "inferno",
                        limits = c(0, upper_from_median(qc_df$count, mult = 8)),
                        oob = scales::squish) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  ggtitle("Fungal UMI") +
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
print(add_bg(p1))
dev.off()

# ---- Fungi detected gene heatmap ----
pdf("Fungi_Gene_heatmap.pdf", width = 8.6, height = 8.6)
p2 <- ggplot(qc_df, aes(x = as.numeric(A), y = as.numeric(B), color = gene_count)) +
  geom_point(shape = 16, size = 4) +
  scale_color_viridis_c(option = "inferno",
                        limits = c(0, upper_from_median(qc_df$gene_count, mult = 5)),
                        oob = scales::squish) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  ggtitle("Fungal detected genes") +
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
print(add_bg(p2))
dev.off()

## =========================
## 4) Seurat clustering on fungal genes only
## =========================
# Make spots as cells: genes x spots
dat <- fungi_expression_matrix
rownames(dat) <- dat$X
dat$X <- NULL

mat <- Matrix(as.matrix(t(dat)), sparse = TRUE)  # features x cells
fungi <- CreateSeuratObject(mat, min.cells = 1, min.features = 1, project = "fungi_only")

# For fungal-only matrix, mitochondrial % usually not meaningful.
# Keep SCTransform without mt regression (safer).
fungi <- SCTransform(fungi, verbose = FALSE)
fungi <- RunPCA(fungi, verbose = FALSE)
fungi <- RunUMAP(fungi, dims = 1:10, verbose = FALSE)
fungi <- FindNeighbors(fungi, dims = 1:10, verbose = FALSE)
fungi <- FindClusters(fungi, resolution = 0.8, verbose = FALSE)

pdf("Dimplot_fungi.pdf", width = 8.6, height = 8.6)
print(DimPlot(fungi, label = TRUE) + NoLegend())
dev.off()

## markers + heatmap
fungi.markers <- FindAllMarkers(fungi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(fungi.markers, file = "fungi_markers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

top10 <- fungi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("Heatmap_fungi_top10.pdf", width = 10, height = 8)
print(DoHeatmap(fungi, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red")))
dev.off()

## =========================
## 5) Spatial cluster map (fungi clusters)
## =========================
ident <- as.character(Idents(fungi))
cl_df <- data.frame(
  X = names(ident),
  cluster = ident,
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(X, into = c("A", "B"), sep = "x", convert = TRUE)

# dynamic color mapping: only use colors for clusters present
present_clusters <- sort(unique(cl_df$cluster))
pal <- cluster_colors[present_clusters]
# if some clusters exceed provided colors, fallback to ggplot hue
if (any(is.na(pal))) {
  pal <- setNames(scales::hue_pal()(length(present_clusters)), present_clusters)
}

pdf("clustering_SCT_fungi_spatial.pdf", width = 8.6, height = 8.6)
pcl <- ggplot(cl_df, aes(x = as.numeric(A), y = as.numeric(B), color = cluster)) +
  geom_point(shape = 16, size = 4) +
  scale_color_manual(values = pal) +
  ggtitle("Fungal clusters") +
  coord_equal(xlim = xlim_grid, ylim = ylim_grid) +
  scale_y_reverse(expand = expansion(mult = c(-0.013, 0.008))) +
  scale_x_continuous(expand = expansion(mult = c(-0.013, -0.013))) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )
print(add_bg(pcl))
dev.off()

message("Done. Outputs: fungi_expression_matrix.tsv, Fungi_UMI_heatmap.pdf, Fungi_Gene_heatmap.pdf, Dimplot_fungi.pdf, Heatmap_fungi_top10.pdf, clustering_SCT_fungi_spatial.pdf, fungi_markers.tsv")
