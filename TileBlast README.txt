# TILEBLAST

A lightweight R tool for discovering nucleotide sequence similarity using BLASTn and visualizing matched regions with checker-style tile plots.

---

## Overview

**TILEBLAST** automates remote BLASTn searches against the NCBI nucleotide database (`nt`), filters and retrieves the resulting sequences, and visualizes each match as a tile plot, with the alignment region highlighted.

It is ideal for:

* Exploring nucleotide-level similarity (e.g., potential homologues or conserved regions)
* Visually comparing matched sequences
* Generating layout-aware, publication-ready plots

---

## Features

* Accepts raw DNA or NCBI accession as input
* Uses remote BLASTn with customizable filters (identity, E-value, coverage, etc.)
* Retrieves full sequences for matched hits
* Produces one tile plot per match with alignment region highlighted
* Optional debug plots for reverse and complement comparison
* Saves FASTA files, metadata, and plot images
* Fully layout-sensitive (supports variable grid dimensions)

---

## Installation

This script requires R and a few core packages:

```r
install.packages(c("rentrez", "ggplot2", "gridExtra", "Biostrings"))
```

If hosted via GitHub in the future:

```r
devtools::install_github("yourusername/tileblast")
```

---

## Usage

```r
result <- run_blast_and_retrieve(
  input_seq = "NM_205518.2",
  max_hits = 500,
  exclude_predicted = TRUE,
  filter_refseq = TRUE,
  user_defined_columns = 100
)
```

### Output

* `*.png` tile plots
* Debug plots (reverse/complement)
* `.fasta` files of retrieved hits
* `.meta` metadata per plot
* Summary log of filtering and hits

---

## Parameters (Selected)

| Parameter              | Description                                | Default |
| ---------------------- | ------------------------------------------ | ------- |
| `input_seq`            | Raw DNA string or NCBI accession           | —       |
| `max_hits`             | Number of hits to retain *after filtering* | `100`   |
| `exclude_predicted`    | Remove hypothetical/predicted sequences    | `FALSE` |
| `filter_refseq`        | Limit to RefSeq accessions only            | `FALSE` |
| `coverage_threshold`   | Minimum query coverage required            | `0.1`   |
| `user_defined_columns` | Number of columns in tile plot             | `100`   |

---

## Output Structure

```
BLASTn_Plots/
  └── <accession>_<query_name>/
        ├── <hit>.png
        ├── <hit>_RC.png       # (reverse debug)
        ├── <hit>_COMP.png     # (complement debug)
        ├── <hit>.fasta
        ├── <hit>.meta
        └── filter_log.txt
```

---

## Citation & Acknowledgments

This tool uses:

* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/)
  *Altschul et al., 1990, J. Mol. Biol.*
  [https://doi.org/10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836%2805%2980360-2)

* [rentrez](https://cran.r-project.org/package=rentrez) — D. Winter (2017)
  [https://journal.r-project.org/archive/2017/RJ-2017-058/](https://journal.r-project.org/archive/2017/RJ-2017-058/)

* [ggplot2](https://ggplot2.tidyverse.org) — H. Wickham (2016)

---

## License

Released under the [MIT License](LICENSE).

---

## Author & Project

TILEBLAST is part of the broader [TileSuit project](https://zenodo.org/doi/...).
Developed by Jason Fajardo, 2024–2025.
