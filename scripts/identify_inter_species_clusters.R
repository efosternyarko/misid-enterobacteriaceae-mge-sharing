library(dplyr)
library(readr)

mobtyper <- read_tsv(
    "~/shared-team/ronnie.dir/07_mge_sharing/02_mob_cluster/all_mobtyper_results.tsv"
)

# Parse species from isolate ID prefix
clusters <- mobtyper %>%
    rename(isolate_id = isolate) %>%
    mutate(
        species = case_when(
            grepl("^EC_",  isolate_id) ~ "E. coli",
            grepl("^KQ_",  isolate_id) ~ "K. quasipneumoniae",
            grepl("^KV_",  isolate_id) ~ "K. variicola",
            grepl("Kpneumoniae", isolate_id) ~ "K. pneumoniae",
            TRUE ~ "Unknown"
        )
    ) %>%
    filter(!is.na(primary_cluster_id), primary_cluster_id != "-")

# Flag primary clusters spanning >1 species
inter_species <- clusters %>%
    group_by(primary_cluster_id) %>%
    summarise(
        n_plasmids   = n(),
        n_isolates   = n_distinct(isolate_id),
        n_species    = n_distinct(species),
        species_list = paste(sort(unique(species)), collapse = "; ")
    ) %>%
    filter(n_species > 1) %>%
    arrange(desc(n_isolates))

write_csv(inter_species,
    "~/shared-team/ronnie.dir/07_mge_sharing/02_mob_cluster/inter_species_clusters.csv")
print(inter_species)

# Check specifically for Mills et al. dominant K. pneumoniae clusters
mills_clusters <- c("AA274", "AA277", "AA553", "AA556")
clusters %>%
    filter(primary_cluster_id %in% mills_clusters) %>%
    group_by(primary_cluster_id, species) %>%
    summarise(n_isolates = n_distinct(isolate_id), .groups = "drop") %>%
    arrange(primary_cluster_id)
