# Build and plot the inter-species MGE sharing network (plasmid clusters +
# ICE BLAST hits) across all 126 isolates.
#
# FIXES APPLIED (see NOTES_FOR_REVIEW.md for detail):
#  1. Duplicate node bug: ICE-derived isolate labels retained suffix cruft
#     (e.g. "_chromosome_merged") that isolate labels from mob_recon/mobtyper
#     did not, creating phantom duplicate nodes for the same real isolate.
#     Fixed by normalising every isolate reference down to its base ID with
#     a shared helper that strips the species suffix and anything after it.
#  2. Legend label mismatch: `labels=` in scale_edge_colour_manual was an
#     unnamed positional vector, silently swapping which label went with
#     which colour. Fixed with a named vector keyed to the real edge_type
#     values ("plasmid", "ICE").

library(tidygraph)
library(ggraph)
library(igraph)
library(dplyr)
library(readr)
library(stringr)

# Normalise any isolate reference (from mobtyper OR from ICE FASTA headers)
# down to its base isolate ID, regardless of what suffix cruft a given
# pipeline stage attached.
normalise_isolate <- function(x) {
  str_remove(x, "_(Ecoli|Kquasipneumoniae|Kvariicola|Kpneumoniae).*")
}

# --- 1. Build edge list from MOB-suite clusters ---
mobtyper <- read_tsv(
    "~/shared-team/ronnie.dir/07_mge_sharing/02_mob_cluster/all_mobtyper_results.tsv"
)

clusters <- mobtyper %>%
    rename(isolate_id = isolate, cluster_id = primary_cluster_id) %>%
    filter(!is.na(cluster_id), cluster_id != "-") %>%
    mutate(
        species = case_when(
            grepl("^EC_",        isolate_id) ~ "E. coli",
            grepl("^KQ_",        isolate_id) ~ "K. quasipneumoniae",
            grepl("^KV_",        isolate_id) ~ "K. variicola",
            grepl("Kpneumoniae", isolate_id) ~ "K. pneumoniae",
            TRUE ~ "Unknown"
        ),
        isolate_id = normalise_isolate(isolate_id)
    )

# For each cluster, create all pairwise isolate edges
plasmid_edges <- clusters %>%
    group_by(cluster_id) %>%
    filter(n() > 1) %>%
    do({
        isolates <- .$isolate_id
        expand.grid(from = isolates, to = isolates, stringsAsFactors = FALSE) %>%
            filter(from < to) %>%
            mutate(cluster_id = .$cluster_id[1], edge_type = "plasmid")
    }) %>%
    ungroup()

# --- 2. Add ICE edges ---
ice_blast <- read_tsv("03_MGE_sharing/04_ICE_blast/ICE_blast_no_selfhits.tsv",
                      col_names = c("qseqid","sseqid","pident","qcovs","length",
                                    "mismatch","gapopen","qstart","qend",
                                    "sstart","send","evalue","bitscore")) %>%
    mutate(
        from = normalise_isolate(sub("__ICE.*", "", qseqid)),
        to   = normalise_isolate(sub("__ICE.*", "", sseqid))
    ) %>%
    filter(from < to) %>%
    mutate(edge_type = "ICE", cluster_id = paste0("ICE_", qseqid)) %>%
    select(from, to, cluster_id, edge_type, pident, qcovs)

# --- 3. Build node table ---
node_metadata <- clusters %>%
    select(isolate_id, species) %>%
    distinct(isolate_id, .keep_all = TRUE) %>%
    rename(name = isolate_id) %>%
    mutate(
        dataset  = if_else(grepl("Kpneumoniae", name), "Mills et al. 2024", "Misidentified"),
        hospital = sub(".*_H([0-9]+)_.*", "H\\1", name)   # adjust regex to your naming
    )

# --- 4. Combine edges ---
all_edges <- bind_rows(
    plasmid_edges %>% select(from, to, edge_type, cluster_id),
    ice_blast     %>% select(from, to, edge_type, cluster_id)
)

# --- 5. Build tidygraph object ---
g <- tbl_graph(nodes = node_metadata, edges = all_edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(degree = centrality_degree())

# --- 6. Plot ---
species_colours <- c(
    "E. coli"            = "#E41A1C",
    "K. pneumoniae"      = "#377EB8",
    "K. quasipneumoniae" = "#4DAF4A",
    "K. variicola"       = "#984EA3",
    "K. aerogenes"       = "#FF7F00"
)

edge_colours <- c("plasmid" = "#2166AC", "ICE" = "#D6604D")

set.seed(42)
p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(colour = edge_type), alpha = 0.5, width = 0.8) +
    geom_node_point(aes(colour = species, shape = dataset, size = degree)) +
    geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
    scale_colour_manual(values = species_colours) +
    scale_edge_colour_manual(values = edge_colours,
                             name = "MGE type",
                             labels = c("plasmid" = "Plasmid cluster",
                                        "ICE" = "ICE (≥95% identity)")) +
    scale_shape_manual(values = c("Mills et al. 2024" = 16, "Misidentified" = 17),
                       name = "Dataset") +
    scale_size_continuous(range = c(3, 10), name = "Shared MGE degree") +
    labs(
        title    = "Inter-species mobile genetic element sharing network",
        subtitle = "Southern Ghana tertiary hospitals | Misidentified isolates + Mills et al. 2024 K. pneumoniae",
        colour   = "Species",
        caption  = "Edges: shared plasmid cluster (MOB-suite) or shared ICE (BLASTn ≥95% identity, ≥80% coverage)"
    ) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "right")

ggsave("03_MGE_sharing/06_network/inter_species_MGE_network.pdf",
       p, width = 14, height = 10)
ggsave("03_MGE_sharing/06_network/inter_species_MGE_network.png",
       p, width = 14, height = 10, dpi = 300)
