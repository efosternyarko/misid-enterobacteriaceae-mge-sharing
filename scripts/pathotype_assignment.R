# E. coli pathotype assignment from ABRicate VFDB screening results:
# intestinal-marker detection plus Johnson et al. (2018) and Spurbeck et
# al. (2012) ExPEC scoring criteria.

library(tidyverse)

# 1. Load data
vf <- read_tsv("08_virulence/vfdb_ecoli.tsv")
all_files <- list.files("00_assemblies/05_ecoli", pattern = "EC_.*\\.fasta", full.names = TRUE)

# 2. Define pattern regex
# Word boundaries prevent eae matching eaeH, estA matching estABC, etc.
inpec_markers <- "\\beae\\b|\\bbfpA\\b|\\bstx|\\baggR\\b|\\bipa|\\bestA\\b|\\bestB\\b|\\beltA\\b|\\beltB\\b"

# 3. Build master isolate list
isolates_master <- tibble(file_path = all_files) %>%
  mutate(sample_id = str_extract(basename(file_path), "EC_[^.]+")) %>%
  select(sample_id)

# 4. Profile isolates
isolate_profiles <- vf %>%
  mutate(sample_id = str_extract(basename(`#FILE`), "EC_[^.]+")) %>%
  group_by(sample_id) %>%
  summarise(
    Has_Intestinal_Markers = any(str_detect(GENE, inpec_markers)),

    # Criteria A: Johnson et al. (2018) ExPEC Score (threshold: >= 2)
    Johnson_Score = sum(c(
      any(str_detect(GENE, "\\b(papA|papC)\\b")),
      any(str_detect(GENE, "\\b(sfaD|sfaE|sfaS|focC|focI)\\b")),
      any(str_detect(GENE, "\\bafa")),          # prefix match: catches afaC, afaD, afaE etc.
      any(str_detect(GENE, "\\b(iutA|iroN)\\b")),
      any(str_detect(GENE, "\\b(kpsM|kpsT|kpsD)\\b")),
      any(str_detect(GENE, "\\b(sat|hlyA|cnf1)\\b"))
    ), na.rm = TRUE),

    # Criteria B: Spurbeck et al. (2012) UPEC Score (threshold: >= 2)
    Spurbeck_Score = sum(c(
      any(str_detect(GENE, "\\bfim")),          # prefix match: catches fimA, fimH, fimC etc.
      any(str_detect(GENE, "\\b(chuA)\\b")),
      any(str_detect(GENE, "\\b(fyuA)\\b")),
      any(str_detect(GENE, "\\b(vat)\\b")),
      any(str_detect(GENE, "\\b(yfcV)\\b"))
    ), na.rm = TRUE)
  ) %>%
  mutate(
    Is_ExPEC = (Johnson_Score >= 2) | (Spurbeck_Score >= 2)
  )

# 5. Report generation
pathotype_report <- isolates_master %>%
  left_join(isolate_profiles, by = "sample_id") %>%
  mutate(
    Has_Intestinal_Markers = replace_na(Has_Intestinal_Markers, FALSE),
    Johnson_Score          = replace_na(Johnson_Score, 0L),
    Spurbeck_Score         = replace_na(Spurbeck_Score, 0L),
    Is_ExPEC               = replace_na(Is_ExPEC, FALSE),

    Pathotype = case_when(
      Has_Intestinal_Markers & Is_ExPEC ~ "ExPEC with intestinal markers",
      Has_Intestinal_Markers            ~ "Intestinal",
      Is_ExPEC                          ~ "ExPEC",
      TRUE                              ~ "Commensal / Unknown"
    )
  )

write_tsv(pathotype_report, "08_virulence/pathotype_report.tsv")
print(pathotype_report, width = Inf)
