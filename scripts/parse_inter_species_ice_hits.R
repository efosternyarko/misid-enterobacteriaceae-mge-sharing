library(dplyr)
library(readr)

blast <- read_tsv("03_MGE_sharing/04_ICE_blast/ICE_blast_no_selfhits.tsv",
                  col_names = c("qseqid","sseqid","pident","qcovs","length",
                                "mismatch","gapopen","qstart","qend",
                                "sstart","send","evalue","bitscore"))

blast <- blast %>%
    mutate(
        query_isolate   = sub("__ICE.*", "", qseqid),
        subject_isolate = sub("__ICE.*", "", sseqid),
        query_species   = case_when(
            grepl("Ecoli",            query_isolate) ~ "E. coli",
            grepl("Kpneumoniae",      query_isolate) ~ "K. pneumoniae",
            grepl("Kquasipneumoniae", query_isolate) ~ "K. quasipneumoniae",
            grepl("Kvariicola",       query_isolate) ~ "K. variicola",
            grepl("aerogenes",        query_isolate) ~ "K. aerogenes",
            TRUE ~ "Unknown"
        ),
        subject_species = case_when(
            grepl("Ecoli",            subject_isolate) ~ "E. coli",
            grepl("Kpneumoniae",      subject_isolate) ~ "K. pneumoniae",
            grepl("Kquasipneumoniae", subject_isolate) ~ "K. quasipneumoniae",
            grepl("Kvariicola",       subject_isolate) ~ "K. variicola",
            grepl("aerogenes",        subject_isolate) ~ "K. aerogenes",
            TRUE ~ "Unknown"
        )
    ) %>%
    filter(query_isolate != subject_isolate)

inter_species_ICE <- blast %>%
    filter(query_species != subject_species) %>%
    arrange(desc(pident))

write_csv(inter_species_ICE, "03_MGE_sharing/04_ICE_blast/inter_species_ICE_hits.csv")
