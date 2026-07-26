# Cross-reference confirmed plasmid-borne ARGs (from
# plasmid_borne_arg_identification.sh) against the inter-species plasmid
# clusters, to report which AMR genes are confirmed present specifically
# on the plasmids driving inter-species sharing.

library(tidyverse)

MOB_DIR       <- "07_mge_sharing/01_mob_recon"
PLASMID_AMR   <- "07_mge_sharing/08_plasmid_amr"
CLUSTER_TABLE <- "07_mge_sharing/06_network/inter_species_MGE_summary_table.csv"

# --- 1. Load mob_recon contig reports ---
# Maps each contig (sequence_id) to its primary_cluster_id and molecule_type
contig_reports <- list.files(
  MOB_DIR, pattern = "contig_report.txt",
  recursive = TRUE, full.names = TRUE
) |>
  map_dfr(\(f) {
    read_tsv(f, col_types = cols(.default = "c")) |>
      mutate(isolate = basename(dirname(f)))
  })

# --- 2. Load plasmid AMRFinderPlus results ---
plasmid_amr_raw <- list.files(
  PLASMID_AMR, pattern = "_plasmid_amr\\.tsv$",
  full.names = TRUE
) |>
  map_dfr(\(f) {
    read_tsv(f, col_types = cols(.default = "c")) |>
      mutate(isolate = str_remove(basename(f), "_plasmid_amr\\.tsv$"))
  })

# --- 3. Apply PARTIAL hit filters (consistent with full-genome AMRFinder filtering) ---
# PARTIALX             = frameshift detected -> gene likely non-functional -> exclude
# PARTIAL_CONTIG_ENDX  = contig truncation -> real gene, apply identity/coverage threshold
plasmid_amr_filtered <- plasmid_amr_raw |>
  filter(
    !str_detect(Method, "^PARTIALX"),
    !(str_detect(Method, "PARTIAL_CONTIG_END") &
        (as.numeric(`% Coverage of reference sequence`) < 70 |
         as.numeric(`% Identity to reference sequence`) < 90))
  )

# --- 4. Join with contig report to get primary_cluster_id ---
# AMRFinderPlus "Contig id" column matches mob_recon "sequence_id"
plasmid_amr_annotated <- plasmid_amr_filtered |>
  left_join(
    contig_reports |> select(isolate, sequence_id, primary_cluster_id, molecule_type),
    by = c("isolate", "Contig id" = "sequence_id")
  ) |>
  filter(molecule_type == "plasmid")

# --- 5. Filter to inter-species shared plasmids only ---
inter_clusters <- read_csv(CLUSTER_TABLE)

plasmid_amr_inter <- plasmid_amr_annotated |>
  filter(primary_cluster_id %in% inter_clusters$primary_cluster_id)

# --- 6. Summarise: confirmed plasmid ARGs per inter-species cluster ---
plasmid_arg_summary <- plasmid_amr_inter |>
  group_by(primary_cluster_id) |>
  summarise(
    n_isolates_with_plasmid_ARG = n_distinct(isolate),
    confirmed_plasmid_ARGs = paste(sort(unique(`Gene symbol`)), collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(inter_clusters, by = "primary_cluster_id") |>
  relocate(primary_cluster_id, species_list, n_isolate, confirmed_plasmid_ARGs)

write_tsv(
  plasmid_arg_summary,
  file.path(PLASMID_AMR, "inter_species_confirmed_plasmid_ARGs.tsv")
)
