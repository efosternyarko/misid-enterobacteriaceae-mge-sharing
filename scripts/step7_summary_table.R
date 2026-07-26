# Summary table of inter-species MGE-sharing clusters cross-referenced with
# their carried AMR genes.
#
# FIX APPLIED: the placeholder path "path/to/amrfinderplus_all_isolates.tsv"
# is replaced with the real combined AMR file produced by
# aggregate_amr_filtered.sh. Also corrected the AMRFinderPlus column names
# used in the join -- the real output columns are "Gene symbol" and
# "Sequence name" (with spaces, as used throughout Section 2), not the
# snake_case gene_symbol/sequence_name in the original placeholder code.

library(dplyr)
library(readr)

# `inter_species` comes from identify_inter_species_clusters.R (run earlier
# in the same session, or re-load its saved CSV output here if starting fresh):
# inter_species <- read_csv("~/shared-team/ronnie.dir/07_mge_sharing/02_mob_cluster/inter_species_clusters.csv")

amr <- read_tsv("07_mge_sharing/07_amr_summary/all_amr_filtered_combined.tsv")

summary_table <- inter_species %>%
    left_join(amr %>% select(isolate, `Gene symbol`, `Sequence name`),
              by = c("isolate_id" = "isolate")) %>%
    group_by(cluster_id, species_list) %>%
    summarise(
        n_isolates = n_distinct(isolate_id),
        n_species  = n_distinct(species),
        AMR_genes  = paste(sort(unique(`Gene symbol`)), collapse = "; "),
        .groups    = "drop"
    ) %>%
    arrange(desc(n_species), desc(n_isolates))

write_csv(summary_table, "03_MGE_sharing/06_network/inter_species_MGE_summary_table.csv")
