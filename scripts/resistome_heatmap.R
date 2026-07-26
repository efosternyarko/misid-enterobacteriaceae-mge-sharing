# Resistome heatmap: intrinsic vs chromosomal-mediated vs plasmid-mediated
# acquired AMR genes/mutations across the 23 misidentified isolates,
# annotated by species.
#
# NOTE: output filename changed from the original "Hongkong_combined_
# resistome.pdf" to "Ghana_combined_resistome.pdf" to match the project's
# actual setting (Southern Ghana) -- see NOTES_FOR_REVIEW.md.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 1. Define paths and output
dir_path <- "05_amr_filtered"
output_dir <- "10_figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
tsv_files <- list.files(path = dir_path, pattern = "\\.tsv$", full.names = TRUE)

# Extract all non-duplicate plasmid-borne ARGs
plasmid_amr_dir <- "07_mge_sharing/08_plasmid_amr"
plasmid_tsv_files <- list.files(path = plasmid_amr_dir, pattern = "_plasmid_amr\\.tsv$", full.names = TRUE)
plasmid_args <- map_dfr(plasmid_tsv_files, ~ {
  df <- read_tsv(.x, show_col_types = FALSE)
  if ("Element symbol" %in% colnames(df)) {
    df %>% select(`Element symbol`)
  } else {
    tibble(`Element symbol` = character())
  }
}) %>%
  drop_na(`Element symbol`) %>%
  pull(`Element symbol`) %>%
  unique()

custom_row_order <- c("E18", "E17", "C11", "E55", "E66", "C43", "C52", "K34", "K59",
                      "C55", "R5", "E47", "C56", "R18", "R8", "K28", "K31", "C22",
                      "C6", "A02", "K4", "E59", "E60")

custom_mge_col_order <- c("Aminoglycoside", "Beta-lactamase", "ESBL", "AmpC Beta-lactamase",
                          "Macrolide", "MLS", "Phenicol",
                          "Sulfonamide", "Tetracycline", "Trimethoprim",
                          "Quinolone", "Quaternary ammonium compound")

custom_chrom_col_order <- c("Multidrug efflux pump repressor", "Quinolone", "Nitrofuran", "Fosfomycin",
                            "Fosmidomycin", "Beta-lactamase", "Siderophore cephalosporin")

# 2. Broad group classification
classify_gene <- function(symbol) {
  if (is.na(symbol)) return("Unknown")
  if (str_detect(symbol, "^(?i)(blalen|blaokp|fosa|oqxa|oqxb)")) {
    return("Intrinsic")
  }
  if (symbol %in% plasmid_args) {
    return("Plasmid")
  }
  return("Chromosomal")
}

# 3. Isolate metadata
all_raw_names <- tools::file_path_sans_ext(basename(tsv_files))
meta_df <- tibble(raw_name = all_raw_names) %>%
  mutate(
    Isolate = map_chr(raw_name, ~ {
      parts <- unlist(strsplit(.x, "_"))
      if (length(parts) >= 2) parts[2] else .x
    }),
    Species = case_when(
      str_starts(raw_name, "EC") ~ "E. coli",
      str_starts(raw_name, "KV") ~ "K. variicola",
      str_starts(raw_name, "KQ") ~ "K. quasipneumoniae",
      TRUE ~ "Other"
    )
  ) %>%
  select(-raw_name) %>%
  distinct()

process_amr_file_all <- function(file_path) {
  raw_name <- tools::file_path_sans_ext(basename(file_path))
  parts <- unlist(strsplit(raw_name, "_"))
  iso_name <- if (length(parts) >= 2) parts[2] else raw_name

  df <- read_tsv(file_path, show_col_types = FALSE)
  if (all(c("Element symbol", "Class") %in% colnames(df))) {
    return(df %>%
             select(`Element symbol`, Class) %>%
             drop_na(Class, `Element symbol`) %>%
             mutate(Group = map_chr(`Element symbol`, classify_gene)) %>%
             separate_longer_delim(Class, delim = "/") %>%
             distinct(Group, Class) %>%
             mutate(Isolate = iso_name, Presence = 1))
  } else {
    return(tibble(Isolate = iso_name, Group = character(), Class = character(), Presence = numeric()))
  }
}

all_parsed_data <- map_dfr(tsv_files, process_amr_file_all)

# 5. Build group matrices
build_group_matrix <- function(data_df, target_group, col_order_plan = NULL) {
  group_data <- data_df %>%
    filter(Group == target_group, Isolate %in% custom_row_order)

  if (nrow(group_data) > 0) {
    matrix_df <- group_data %>%
      select(-Group) %>%
      pivot_wider(names_from = Class, values_from = Presence, values_fill = list(Presence = 0))
  } else {
    matrix_df <- tibble(Isolate = custom_row_order)
  }

  missing_iso <- setdiff(custom_row_order, matrix_df$Isolate)
  if (length(missing_iso) > 0) {
    matrix_df <- bind_rows(matrix_df, tibble(Isolate = missing_iso))
  }
  matrix_df <- matrix_df %>% mutate(across(-Isolate, ~replace_na(.x, 0)))

  matrix_df <- matrix_df %>%
    mutate(Isolate = factor(Isolate, levels = custom_row_order)) %>%
    arrange(Isolate) %>%
    mutate(Isolate = as.character(Isolate))

  mat <- as.matrix(matrix_df %>% select(-Isolate))
  rownames(mat) <- matrix_df$Isolate

  if (!is.null(col_order_plan) && ncol(mat) > 0) {
    matched <- intersect(col_order_plan, colnames(mat))
    unmatched <- setdiff(colnames(mat), col_order_plan)
    mat <- mat[, c(matched, unmatched), drop = FALSE]
  }
  return(mat)
}

mat_intrinsic   <- build_group_matrix(all_parsed_data, "Intrinsic")
mat_chromosomal <- build_group_matrix(all_parsed_data, "Chromosomal", custom_chrom_col_order)
mat_mge         <- build_group_matrix(all_parsed_data, "Plasmid", custom_mge_col_order)

# 6. Visual settings
colors <- colorRamp2(c(0, 1), c("white", "#F27A6D"))
species_colors <- c("E. coli" = "#6DC6F2", "K. variicola" = "#6D6DF2", "K. quasipneumoniae" = "#6D99F2", "Other" = "darkgray")

aligned_metadata <- tibble(Isolate = custom_row_order) %>%
  left_join(meta_df, by = "Isolate") %>%
  mutate(Species = replace_na(Species, "Other"))

row_ha <- rowAnnotation(
  Species = aligned_metadata$Species,
  col = list(Species = species_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Species = list(title = "Species", title_gp = gpar(fontface = "bold"), labels_gp = gpar(fontface = "italic")))
)

# 7. Define sub-heatmaps
title_headroom <- "\n\n\n\n\n\n"

ht_intrinsic <- Heatmap(
  mat_intrinsic, name = "Resistance", col = colors, rect_gp = gpar(col = "white", lwd = 0.5),
  width = ncol(mat_intrinsic) * unit(6, "mm"), height = nrow(mat_intrinsic) * unit(6, "mm"),
  column_title = title_headroom,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  left_annotation = row_ha, row_names_side = "left", show_row_names = TRUE,
  column_names_side = "top", column_names_rot = 90,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(at = c(0, 1), labels = c("Absent", "Present"), title = "Gene/mutation status", color_bar = "discrete")
)

ht_chromosomal <- Heatmap(
  mat_chromosomal, name = "Chromosomal_Heat", col = colors, rect_gp = gpar(col = "white", lwd = 0.5),
  width = ncol(mat_chromosomal) * unit(6, "mm"), height = nrow(mat_chromosomal) * unit(6, "mm"),
  column_title = title_headroom,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  show_row_names = FALSE, column_names_side = "top", column_names_rot = 90,
  show_heatmap_legend = FALSE
)

ht_mge <- Heatmap(
  mat_mge, name = "MGE_Heat", col = colors, rect_gp = gpar(col = "white", lwd = 0.5),
  width = ncol(mat_mge) * unit(6, "mm"), height = nrow(mat_mge) * unit(6, "mm"),
  column_title = title_headroom,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  show_row_names = FALSE, column_names_side = "top", column_names_rot = 90,
  show_heatmap_legend = FALSE
)

final_ht_list <- ht_intrinsic + ht_chromosomal + ht_mge

# 8. Output layout execution
combined_pdf_path <- file.path(output_dir, "Ghana_combined_resistome.pdf")
pdf(combined_pdf_path, width = 16, height = 12)

draw(final_ht_list, padding = unit(c(5, 5, 5, 5), "mm"), heatmap_legend_side = "right")

# --- Custom grid plot overlays ---
draw_inner_bracket <- function(heatmap_name, label_text) {
  decorate_column_title(heatmap_name, {
    grid.lines(x = c(0, 1), y = c(0.35, 0.35), gp = gpar(lwd = 1.0, col = "black"))
    grid.lines(x = c(0, 0), y = c(0.35, 0.10), gp = gpar(lwd = 1.0, col = "black"))
    grid.lines(x = c(1, 1), y = c(0.35, 0.10), gp = gpar(lwd = 1.0, col = "black"))
    grid.text(label = label_text, x = 0.5, y = 0.50, gp = gpar(fontface = "plain", fontsize = 10))
  })
}

draw_inner_bracket("Chromosomal_Heat", "Chromosomal-mediated")
draw_inner_bracket("MGE_Heat", "Plasmid-mediated")

draw_group_bracket <- function(heatmap_name, label_text) {
  decorate_column_title(heatmap_name, {
    grid.lines(x = c(0, 1), y = c(0.75, 0.75), gp = gpar(lwd = 1.5, col = "black"))
    grid.lines(x = c(0, 0), y = c(0.75, 0.65), gp = gpar(lwd = 1.5, col = "black"))
    grid.lines(x = c(1, 1), y = c(0.75, 0.65), gp = gpar(lwd = 1.5, col = "black"))
    grid.text(label = label_text, x = 0.5, y = 0.90, gp = gpar(fontface = "bold", fontsize = 10))
  })
}

draw_group_bracket("Resistance", "Intrinsic")

decorate_column_title("Chromosomal_Heat", {
  grid.lines(x = c(0, 1), y = c(0.75, 0.75), gp = gpar(lwd = 1.5, col = "#bd0026"))
  grid.lines(x = c(0, 0), y = c(0.75, 0.65), gp = gpar(lwd = 1.5, col = "#bd0026"))
})

decorate_column_title("MGE_Heat", {
  grid.lines(x = c(0, 1), y = c(0.75, 0.75), gp = gpar(lwd = 1.5, col = "#bd0026"))
  grid.lines(x = c(1, 1), y = c(0.75, 0.65), gp = gpar(lwd = 1.5, col = "#bd0026"))
})

decorate_column_title("Chromosomal_Heat", {
  grid.lines(x = c(1, 1.1), y = c(0.75, 0.75), gp = gpar(lwd = 1.5, col = "#bd0026"))
})

decorate_column_title("MGE_Heat", {
  grid.text(
    label = "                                Acquired resistance genes/mutations",
    x = -0.15,
    y = 0.90,
    gp = gpar(fontface = "bold", fontsize = 10, col = "#bd0026")
  )
})

dev.off()
