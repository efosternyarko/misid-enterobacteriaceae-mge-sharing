# Genomic Characterisation of Misidentified Clinical Enterobacteriaceae: Inter-Species Mobile Genetic Element Sharing

Analysis protocol and scripts for an MSc dissertation project characterising clinical Enterobacteriaceae isolates from tertiary hospitals in Southern Ghana that were misidentified at initial species assignment. After species correction (Kraken2 contig-level re-identification), the confirmed dataset comprises 15 *Escherichia coli*, 6 *Klebsiella quasipneumoniae*, and 2 *Klebsiella variicola* isolates, analysed alongside the 103 *K. pneumoniae* assemblies from Mills et al. 2024 (126 isolates total).

The central analysis determines whether plasmids and integrative conjugative elements (ICEs) are shared across species boundaries within the same hospital environment.

**Student:** Cheuk Kiu Liu (Ronnie)
**Supervisor:** Dr E. Foster-Nyarko

## Contents

- [`protocol.qmd`](protocol.qmd) — the full analysis protocol (Quarto document), covering assembly QC and species re-identification, typing and annotation, AMR/plasmid detection, the inter-species MGE-sharing network analysis, plasmid-borne ARG confirmation, and species-specific virulence characterisation. Render with `quarto render protocol.qmd`.
- [`scripts/`](scripts/) — standalone bash/R/Python scripts referenced by the protocol, runnable independently.
- [`NOTES_FOR_REVIEW.md`](NOTES_FOR_REVIEW.md) — a record of fixes applied and points flagged for verification before this protocol is treated as final.

## Citation

If citing this protocol in a dissertation Methods section, reference it as:

> Liu, C.K. Genomic characterisation of misidentified clinical Enterobacteriaceae isolates: inter-species mobile genetic element sharing. Analysis protocol. Available at: https://github.com/efosternyarko/misid-enterobacteriaceae-mge-sharing

(Update the URL/access date once this repository is published.)
