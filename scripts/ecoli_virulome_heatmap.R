# E. coli iron-acquisition/virulence gene heatmap, annotated with ST,
# phylogroup, pathotype, O:H serotype, and K-locus metadata.
#
# NOTE: output filename changed from the original "Hongkong_iron_heatmap.pdf"
# to "Ghana_ecoli_iron_heatmap.pdf" to match the project's actual setting
# (Southern Ghana) -- see NOTES_FOR_REVIEW.md.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 1. Load ABRicate data
raw_data <- read_tsv("08_virulence/vfdb_ecoli_matrix.tsv")

cleaned_data <- raw_data %>%
  mutate(Isolate = gsub("^EC_", "", gsub("\\.fasta$", "", basename(`#FILE`)))) %>%
  select(-`#FILE`, -NUM_FOUND) %>%
  column_to_rownames("Isolate")

# 2. Filter iron genes
iron_regex_pattern <- "^(chu|ent|fep|fes|iuc|iut|irp|fyu|iro|sit|ybt|hma)"
iron_gene_columns <- grep(iron_regex_pattern, colnames(cleaned_data), ignore.case = TRUE, value = TRUE)
if ("clbA" %in% colnames(cleaned_data)) {
  iron_gene_columns <- c(iron_gene_columns, "clbA")
}

df_iron <- cleaned_data[, iron_gene_columns, drop = FALSE]
mat_gene <- as.matrix(df_iron)
mat_gene[mat_gene == "."] <- "0"
mat_gene[mat_gene != "0"] <- "1"
class(mat_gene) <- "numeric"

collapse_prefixes <- c("chu", "ent", "fep", "iro", "iuc", "sit", "ybt")
mat_collapsed <- as.data.frame(mat_gene)

for (locus in collapse_prefixes) {
  locus_cols <- grep(paste0("^", locus), colnames(mat_collapsed), value = TRUE)
  if (length(locus_cols) > 0) {
    mat_collapsed[[locus]] <- ifelse(rowSums(mat_collapsed[, locus_cols, drop = FALSE]) > 0, 1, 0)
    mat_collapsed <- mat_collapsed[, !(colnames(mat_collapsed) %in% locus_cols), drop = FALSE]
  }
}
mat_gene <- as.matrix(mat_collapsed)

# 3. Group genes by function
system_list <- list(
  "Haem uptake and utilisation" = c(grep("^chu", colnames(mat_gene), value=TRUE), grep("^hma", colnames(mat_gene), value=TRUE)),
  "Colibactin"                  = grep("^clb", colnames(mat_gene), value=TRUE),
  "Enterobactin system"         = grep("^(ent|fep|fes)", colnames(mat_gene), value=TRUE),
  "Salmochelin system"          = grep("^iro", colnames(mat_gene), value=TRUE),
  "Aerobactin system"           = grep("^(iuc|iut)", colnames(mat_gene), value=TRUE),
  "Sit system"                  = grep("^sit", colnames(mat_gene), value=TRUE),
  "Yersiniabactin system"       = c(grep("^ybt", colnames(mat_gene), value=TRUE), grep("^irp", colnames(mat_gene), value=TRUE), grep("^fyu", colnames(mat_gene), value=TRUE))
)

ordered_cols <- c()
col_groups <- c()
for (sys_name in names(system_list)) {
  genes_in_sys <- intersect(system_list[[sys_name]], colnames(mat_gene))
  genes_in_sys <- sort(genes_in_sys)
  if (length(genes_in_sys) > 0) {
    ordered_cols <- c(ordered_cols, genes_in_sys)
    col_groups <- c(col_groups, rep(sys_name, length(genes_in_sys)))
  }
}
mat_final <- mat_gene[, ordered_cols, drop = FALSE]
group_factor <- factor(col_groups, levels = names(system_list))

# 4. Define metadata & custom order
metadata_raw <- data.frame(
  Isolate  = c("C11", "C43", "C52", "C55", "C56", "E17", "E18", "E47", "E55", "E66", "K34", "K59", "R18", "R5", "R8"),
  ST       = c("ST410", "ST131", "ST131", "ST167", "ST127", "ST648", "ST648", "ST12",  "ST131", "ST131", "ST68", "ST569", "ST127", "ST1642", "ST127"),
  Phylogroup = c("C", "B2", "B2", "A", "B2", "F", "F", "B2", "B2", "B2", "D", "B2", "B2", "B1", "B2"),
  Pathotype = c("ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "ExPEC", "unknown*", "ExPEC"),
  O_type   = c("O8", "O25", "O153", "O89", "O6",  "O8",  "O8",  "O4",  "O25", "O25", "O51", "O46", "O6",  "O8/O71", "O6"),
  H_type   = c("H9", "H4", "H4",  "H10", "H31", "H4", "H4", "H1", "H4", "H4", "H6",  "H31", "H31", "H7", "H31"),
  K_locus  = c("KL302", "KL20", "KL2", "KL767", "KL2", "KL301", "KL301", "KL12", "KL114", "KL114", "KL123", "KL1", "KL2", "KL301", "KL14"),
  stringsAsFactors = FALSE
) %>% column_to_rownames("Isolate")

# Custom order (better organisation)
custom_order <- c("E47", "R8", "C56", "R18", "K34", "K59", "E17", "E18", "E66", "E55", "C43", "C52", "C11", "C55", "R5")

mat_final <- mat_final[custom_order, , drop = FALSE]
metadata_ordered <- metadata_raw[custom_order, ] %>%
  rownames_to_column("Isolate") %>%
  rename(`Sequence type` = ST, `O-type` = O_type, `H-type` = H_type, `K locus` = K_locus)

# 5. Generate heatmaps
heatmap_colors <- c("0" = "#F7F7F7", "1" = "#967BB6")

ht_metadata <- Heatmap(
  as.matrix(metadata_ordered),
  name = "metadata",
  col = c("dummy" = "white"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  show_column_names = FALSE,
  column_split = factor(colnames(metadata_ordered), levels = colnames(metadata_ordered)),
  column_title_gp = gpar(fontsize = 10, fontface = "plain"),
  column_title_rot = 90,
  column_gap = unit(5, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
    grid.text(as.matrix(metadata_ordered)[i, j], x, y, gp = gpar(fontsize = 9))
  },
  width = unit(12, "cm")
)

ht_genes <- Heatmap(
  mat_final,
  name = " ",
  col = heatmap_colors,
  cluster_rows = FALSE,
  show_row_dend = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_split = group_factor,
  column_title_rot = 90,
  column_gap = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 10, fontface = "plain"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  width = unit(10, "cm"),
  height = unit(10, "cm"),
  column_names_gp = gpar(fontsize = 9, fontface = "italic", rot = 90, hjust = 1),
  column_names_side = "top",
  heatmap_legend_param = list(labels = c("Gene absent", "Gene present"), at = c(0, 1))
)

ht_list <- ht_metadata + ht_genes

pdf("10_figures/Ghana_ecoli_iron_heatmap.pdf", width = 18, height = 8)
draw(ht_list)

for (i in 1:ncol(metadata_ordered)) {
  seekViewport(paste0("metadata_column_title_", i))
  grid.lines(x = c(0, 1), y = c(0, 0), gp = gpar(lwd = 1))
  grid.lines(x = c(0, 0), y = c(0, -0.025), gp = gpar(lwd = 1))
  grid.lines(x = c(1, 1), y = c(0, -0.025), gp = gpar(lwd = 1))
}
for (i in 1:length(levels(group_factor))) {
  seekViewport(paste0(" _column_title_", i))
  grid.lines(x = c(0, 1), y = c(0, 0), gp = gpar(lwd = 1))
  grid.lines(x = c(0, 0), y = c(0, -0.025), gp = gpar(lwd = 1))
  grid.lines(x = c(1, 1), y = c(0, -0.025), gp = gpar(lwd = 1))
}
dev.off()
