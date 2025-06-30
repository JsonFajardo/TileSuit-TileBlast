# ================================================================
# TILEBLAST: Remote BLASTn Search + Tile Plot Visualization Tool
# ================================================================
#
# Purpose:
# This script automates the process of detecting similar nucleotide sequences
# using BLASTn against the NCBI database, applying quality filters, retrieving
# the full matched sequences, and generating visual tile plots that highlight
# the aligned region.
#
# Key Features:
# - Accepts either an NCBI accession number or raw DNA sequence as input.
# - Performs a remote BLASTn search and optionally enforces a strict quota
#   of high-quality hits (based on identity, coverage, E-value, etc.).
# - Downloads and stores full sequences for matched hits.
# - Generates tile plots for each hit with optional debug views for reverse
#   complement comparison.
#
# Goals:
# - Discover nucleotide sequences that are similar to a given input.
# - Enable visual analysis of sequence structure and alignment context.
# - Support exploratory and comparative sequence analysis via structured,
#   layout-sensitive plots.
#
# Usage Notes:
# - The number of returned hits depends on both NCBI's internal scoring and
#   user-defined filters. To retrieve more curated results:
#     • Set `filter_refseq = TRUE` to limit to RefSeq accessions.
#     • Set `exclude_predicted = TRUE` to remove hypothetical/predicted entries.
#     • Set `enforce_hit_quota = TRUE` to ensure the final number of hits matches
#       `max_hits`, even if many must be scanned.
# - The default BLAST database is `"nt"` but can be changed internally in
#   `run_remote_blast()`.
#
# Output:
# - PNG tile plots per hit, saved in a structured output directory.
# - Metadata summary file and optional FASTA files for each hit.
# - Optional debug plots (_RC and _COMP) for visual strand comparison.
#
# Recommended Use Case:
# Use this tool to search for nucleotide-level sequence similarity (e.g., potential
# homologues, conserved regions, or structural patterns) and to visualize the results
# in a spatially interpretable, layout-aware format.

# --------------------------------------------------
# References and Tools Cited
# --------------------------------------------------
# NCBI BLAST+ — Altschul et al. (1990)
# https://doi.org/10.1016/S0022-2836(05)80360-2
#
# rentrez — Winter (2017), The R Journal
# https://journal.r-project.org/archive/2017/RJ-2017-058/
#
# ggplot2 — Wickham (2016), Springer
# https://ggplot2.tidyverse.org


library(httr)
library(XML)
library(rentrez)
library(ggplot2)
library(reshape2)

# --- Detect sequence type ---

# --- detect_seq_type(input) ---
# Purpose:
# Determines whether the input is a valid nucleotide accession or a raw DNA sequence,
# and extracts associated metadata if possible.
#
# Behavior:
# - Accepts either:
#   - An NCBI accession string (e.g. "NM_000518.5"), or
#   - A raw nucleotide sequence (e.g. "ATGCGTA...")
# - Returns a list with:
#   - type: "nucleotide" or "unknown"
#   - seq: the cleaned sequence
#   - accession: detected accession or "RawSeq"
#   - name: gene name (if present in header)
#   - header: full FASTA header (if available)
#
# Logic:
# - Fetches FASTA if input looks like an accession.
# - Validates sequence (A/T/G/C/U/N only).
# - Rejects invalid or ambiguous input.
#
# Used in:
# run_blast_and_retrieve() (initial input processing)

detect_seq_type <- function(input) {
  clean <- gsub("\\s", "", input)
  if (grepl("^[A-Z]{1,2}_?\\d+", input)) {
    res <- tryCatch(GET(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", input, "&rettype=fasta")), error=function(e) NULL)
    txt <- rawToChar(res$content)
    if (grepl("^>", txt)) {
      lines <- strsplit(txt, "\\n")[[1]]
      header_line <- sub("^>", "", lines[1])
      accession <- sub("^(\\S+).*", "\\1", header_line)
      name <- sub("^\\S+\\s+", "", header_line)
      seq <- paste(tail(lines, -1), collapse = "")
      if (!grepl("[^ATGCUN]", seq, ignore.case=TRUE)) {
        return(list(type = "nucleotide", seq = seq, accession = accession, name = name, header = header_line))
      }
    }
    return(list(type="unknown", seq=NA, accession=NA, name=NA, header=NA))
  }
  if (nchar(clean) < 10 || grepl("[^A-Za-z]", clean)) return(list(type="unknown", seq=NA, accession=NA, name=NA, header=NA))
  if (!grepl("[^ATGCUN]", clean, ignore.case=TRUE)) return(list(type="nucleotide", seq=clean, accession="RawSeq", name="", header="Raw nucleotide"))
  return(list(type="unknown", seq=NA, accession=NA, name=NA, header=NA))
}

# --- Complement only ---

# --- complement_only(seq) ---
# Purpose:
# Returns the DNA complement of a given nucleotide sequence (A <-> T, G <-> C).
#
# Behavior:
# - Accepts a DNA sequence as a string.
# - Converts it to uppercase and replaces each nucleotide with its complement:
#   - A → T
#   - T → A
#   - G → C
#   - C → G
#
# Logic:
# - Uses chartr() for fast character-level substitution.
# - Does not reverse the sequence.
#
# Used in:
# - Debugging visualizations (e.g., _COMP plots) to show the complementary region in the forward strand.

complement_only <- function(seq) {
  chartr("ATGC", "TACG", toupper(seq))
}

# --- Run BLAST remotely and get XML ---

# --- NOTE on BLAST database selection ---
# The BLAST search uses the "nt" (nucleotide collection) database by default.
# If you want to limit results to manually curated sequences, consider changing the
# database in `run_remote_blast()` from "nt" to one of the following:
# 
#   - "refseq_rna"     : curated transcript records (e.g., NM_ or NR_)
#   - "refseq_genomic" : curated genomic records
#   - "est"            : expressed sequence tags
#   - "gss"            : genome survey sequences
#
# To do so, edit the `POST(...)` call inside the `run_remote_blast()` function.

# --- run_remote_blast(seq_info, max_attempts = 30, retry_delay = 3) ---
# Purpose:
# Submits a BLASTn query to NCBI and retrieves the XML result when ready.
#
# Behavior:
# - Posts the sequence to NCBI BLAST and retrieves the Request ID (RID).
# - Polls the server periodically until the result is ready or max_attempts is reached.
# - Parses and returns the XML content upon success.
#
# Parameters:
# - seq_info: A list from detect_seq_type() containing the query sequence.
# - max_attempts: Max number of polling attempts (default = 30).
# - retry_delay: Seconds to wait between attempts or before initial poll.
#
# Logic:
# - Sends POST request with query sequence to initiate BLAST.
# - Waits for `RTOE` seconds or `retry_delay`, whichever is valid.
# - Polls using GET with RID until XML is available or timeout.
# - Handles temporary HTML errors with retry logic.
#
# Used in:
# - run_blast_and_retrieve() to initiate the remote BLAST step.


run_remote_blast <- function(seq_info, max_attempts = 30, retry_delay = 3, hitlist_size = 1000) {
  blast_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
  
  post <- POST(blast_url, body = list(
    CMD = "Put",
    PROGRAM = "blastn",
    DATABASE = "nt",
    HITLIST_SIZE = hitlist_size,
    QUERY = seq_info$seq
  ), encode = "form")
  
  post_text <- rawToChar(post$content)
  rid_match <- regmatches(post_text, regexpr("RID = \\S+", post_text))
  rtoe_match <- regmatches(post_text, regexpr("RTOE = \\d+", post_text))
  
  rid <- sub("RID = ", "", rid_match)
  rtoe <- if (length(rtoe_match) > 0) as.numeric(sub("RTOE = ", "", rtoe_match)) else retry_delay
  if (is.na(rtoe) || rtoe < 0) rtoe <- retry_delay
  Sys.sleep(rtoe)
  
  attempts <- 0
  repeat {
    res <- GET(blast_url, query = list(
      CMD = "Get",
      FORMAT_TYPE = "XML",
      RID = rid,
      HITLIST_SIZE = hitlist_size
    ))
    
    content <- rawToChar(res$content)
    if (startsWith(content, "<?xml")) {
      xml <- tryCatch(xmlParse(content), error = function(e) NULL)
      if (!is.null(xml)) break
    } else {
      warning("Non-XML response from BLAST (e.g. HTML error). Retrying...")
    }
    
    attempts <- attempts + 1
    if (attempts > max_attempts) stop(paste0("BLAST polling failed after ", max_attempts, " attempts."))
    Sys.sleep(3)
  }
  
  xml
}


# --- Filter results ---

# --- filter_hits(xml, seq_info, min_pid, max_evalue, min_length, coverage_threshold, max_hits, exclude_predicted, filter_refseq) ---
# Purpose:
# Filters raw BLAST hits based on user-defined criteria and returns a curated list of matches.
#
# Behavior:
# - Parses BLAST XML results.
# - Applies identity, E-value, alignment length, and coverage filters.
# - Optionally excludes predicted entries or non-RefSeq accessions.
# - Logs a detailed breakdown of inclusion/exclusion statistics.
# - Ensures the final number of hits equals max_hits if enough pass the filters.
#
# Parameters:
# - xml: Parsed BLAST result XML.
# - seq_info: Output of detect_seq_type(), includes query sequence and metadata.
# - min_pid: Minimum percent identity (e.g., 60).
# - max_evalue: Maximum acceptable E-value (e.g., 1e-2).
# - min_length: Minimum alignment length.
# - coverage_threshold: Minimum fraction of query length covered.
# - max_hits: Desired number of final hits after filtering.
# - exclude_predicted: Exclude hits with "PREDICTED", "hypothetical", etc.
# - filter_refseq: If TRUE, only includes hits from RefSeq (accessions with underscore).
#
# Used in:
# - run_blast_and_retrieve() to curate results before downloading sequences.

filter_hits <- function(xml, seq_info, min_pid, max_evalue, min_length, coverage_threshold, max_hits, exclude_predicted, filter_refseq) {
  n_total <- 0
  n_excluded_predicted <- 0
  n_excluded_pid <- 0
  n_excluded_eval <- 0
  n_excluded_length <- 0
  n_excluded_coverage <- 0
  n_passed <- 0
  
  hits <- list()
  all_hits <- getNodeSet(xml, "//Hit")
  
  for (hit in all_hits) {
    n_total <- n_total + 1
    hsp <- hit[["Hit_hsps"]][["Hsp"]]
    pid <- as.numeric(xmlValue(hsp[["Hsp_identity"]])) * 100 / as.numeric(xmlValue(hsp[["Hsp_align-len"]]))
    eval <- as.numeric(xmlValue(hsp[["Hsp_evalue"]]))
    len <- as.numeric(xmlValue(hsp[["Hsp_align-len"]]))
    hstart <- as.numeric(xmlValue(hsp[["Hsp_hit-from"]]))
    hend <- as.numeric(xmlValue(hsp[["Hsp_hit-to"]]))
    cov <- len / nchar(seq_info$seq)
    def <- xmlValue(hit[["Hit_def"]])
    acc <- xmlValue(hit[["Hit_accession"]])
    
    reason <- ""
    
    if (filter_refseq && !grepl("^(NM_|NR_|NG_)", acc)) {
      next
    }
    if (exclude_predicted && grepl("(PREDICTED|unnamed|hypothetical|LOW QUALITY)", def, ignore.case = TRUE)) {
      n_excluded_predicted <- n_excluded_predicted + 1
      reason <- "predicted"
    } else if (pid < min_pid) {
      n_excluded_pid <- n_excluded_pid + 1
      reason <- "pid"
    } else if (eval > max_evalue) {
      n_excluded_eval <- n_excluded_eval + 1
      reason <- "evalue"
    } else if (len < min_length) {
      n_excluded_length <- n_excluded_length + 1
      reason <- "length"
    } else if (cov < coverage_threshold) {
      n_excluded_coverage <- n_excluded_coverage + 1
      reason <- "coverage"
    } else {
      n_passed <- n_passed + 1
      hits[[length(hits) + 1]] <- list(
        title = def,
        accession = acc,
        start = hstart,
        end = hend,
        pid = pid,
        eval = eval,
        len = len
      )
      if (n_passed >= max_hits) break
      next
    }
    
    cat(sprintf("Excluded hit: %s | PID=%.1f | Eval=%.2g | Len=%d | Cov=%.2f | Reason: %s\n",
                acc, pid, eval, len, cov, reason))
  }
  
  summary <- list(
    total = n_total,
    excluded_predicted = n_excluded_predicted,
    excluded_pid = n_excluded_pid,
    excluded_eval = n_excluded_eval,
    excluded_length = n_excluded_length,
    excluded_coverage = n_excluded_coverage,
    passed = n_passed
  )
  
  # Save summary log to file
  name_clean <- gsub("[^A-Za-z0-9._-]", "_", seq_info$name)
  log_dir <- file.path("BLASTn_Plots", paste0(seq_info$accession, "_", name_clean))
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  log_file <- file.path(log_dir, "filter_log.txt")
  writeLines(c(
    paste0("Total hits processed: ", n_total),
    paste0("Excluded (predicted): ", n_excluded_predicted),
    paste0("Excluded (low identity): ", n_excluded_pid),
    paste0("Excluded (high E-value): ", n_excluded_eval),
    paste0("Excluded (short alignment): ", n_excluded_length),
    paste0("Excluded (low coverage): ", n_excluded_coverage),
    paste0("Final hits retained: ", n_passed)
    ), log_file)
  
  return(list(hits = hits, summary = summary))
}


# --- Retrieve full sequences ---

# --- retrieve_full_sequences(hits, save_fasta = TRUE, output_dir = "BLASTn_Plots") ---
# Purpose:
# Downloads the full nucleotide sequences for a list of filtered BLAST hits from NCBI.
#
# Behavior:
# - Uses each hit’s accession to fetch its full FASTA record via the Entrez API.
# - Parses and stores each sequence along with its header metadata.
# - Optionally saves each FASTA to a file in the specified output directory.
#
# Parameters:
# - hits: List of filtered hit metadata (from filter_hits()).
# - save_fasta: Whether to save the retrieved sequences to disk (default = TRUE).
# - output_dir: Directory where FASTA files are stored.
#
# Logic:
# - Skips hits that fail to retrieve or parse.
# - Creates output directory if it doesn't exist.
# - Saves each FASTA with safe filenames derived from accession and name.
#
# Used in:
# - run_blast_and_retrieve() to obtain sequences for tile plotting.

retrieve_full_sequences <- function(hits, save_fasta = TRUE, output_dir = "BLASTn_Plots") {
  seq_list <- list()
  for (hit in hits) {
    fasta <- tryCatch(entrez_fetch(db="nucleotide", id=hit$accession, rettype="fasta", retmode="text"), error=function(e) NA)
    if (!is.na(fasta)) {
      lines <- strsplit(fasta, "\n")[[1]]
      header <- gsub("^>", "", lines[1])
      accession <- sub("^(\\S+).*", "\\1", header)
      name <- sub("^\\S+\\s+", "", header)
      hit$accession <- accession
      hit$name <- name
      hit$header <- header
      hit$sequence <- paste(lines[-1], collapse = "")
      
      if (save_fasta) {
        fasta_lines <- c(
          paste0(">", header),
          gsub("(.{1,70})", "\\1\n", hit$sequence, perl = TRUE)
        )
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        safe_accession <- gsub("[[:space:]/\\\\:*?\"<>|]", "_", hit$accession)
        safe_name <- gsub("[[:space:]/\\\\:*?\"<>|]", "_", hit$name)
        fasta_file <- file.path(output_dir, paste0(safe_accession, "_", safe_name, ".fasta"))
        writeLines(fasta_lines, fasta_file)
      }
      seq_list[[length(seq_list) + 1]] <- hit
    }
  }
  seq_list
}

# --- Helper to reverse complement a DNA sequence ---

# --- reverse_complement(seq) ---
# Purpose:
# Computes the reverse complement of a nucleotide sequence.
#
# Behavior:
# - Translates each base to its complement using standard DNA pairing:
#   A <-> T, G <-> C
# - Reverses the resulting complemented sequence.
#
# Parameters:
# - seq: A character string representing a DNA sequence.
#
# Returns:
# - A new character string that is the reverse complement of the input.
#
# Used in:
# - run_blast_and_retrieve() for generating reverse strand comparison plots.

reverse_complement <- function(seq) {
  comp <- chartr("ATGC", "TACG", toupper(seq))
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

# --- Compute Highlight Thickness ---

# --- compute_highlight_thickness(num_rows, num_cols) ---
# Purpose:
# Dynamically adjusts highlight border thickness based on plot size and shape.
#
# Behavior:
# - Computes a scaling factor based on the total grid size and shape ratio.
# - Penalizes highly unbalanced or very large plots with thinner lines.
#
# Parameters:
# - num_rows: Number of rows in the tile plot grid.
# - num_cols: Number of columns in the tile plot grid.
#
# Returns:
# - A numeric value representing line thickness for `geom_tile()` borders.
#   (Clipped between 0.003 and 0.4)
#
# Used in:
# - create_dna_plot() to determine highlight border width for alignment region.

compute_highlight_thickness <- function(num_rows, num_cols) {
  total_dim <- num_rows + num_cols
  ratio <- max(num_rows / num_cols, num_cols / num_rows)
  falloff <- exp(-0.001 * total_dim) * (1 / (1 + (ratio - 1)^1.5))
  thickness <- 0.9 * falloff
  return(min(max(thickness, 0.05), 0.9))
}


# --- Plotting Function ---

### create_dna_plot(sequence, header, name, accession, num_columns, metadata, output_dir, highlight_range)
# Purpose:
# Creates a colored tile plot of a DNA sequence, highlighting a specific alignment region.
#
# Behavior:
# - Converts a DNA sequence into a 2D grid of tiles (A, T, G, C, N) with color mapping.
# - Highlights a specific range in the sequence (e.g. alignment region).
# - Saves the plot as a PNG image in a structured output directory.
#
# Parameters:
# - sequence: Nucleotide sequence string (A/T/G/C/N).
# - header: FASTA header string (used for metadata and file naming).
# - name: Gene or sequence name (used in file naming).
# - accession: Accession number (used in file naming).
# - num_columns: Number of columns in the tile grid.
# - metadata: Named list of metadata (e.g., Query, Strand, PID, etc.) for title display.
# - output_dir: Directory where the plot will be saved.
# - highlight_range: Numeric vector of length 2 (start, end) specifying the alignment to highlight.
#
# Returns:
# - Saves a tile plot to the output directory. Does not return a value.
#
# Used in:
# - `run_blast_and_retrieve()` to visualize each matched sequence.

create_dna_plot <- function(sequence, header, name, accession, num_columns = NULL, metadata = list(), output_dir = "BLASTn_Plots", highlight_range = NULL) {
  dna_sequence <- gsub("\\s", "", sequence)
  dna_sequence <- toupper(dna_sequence)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  seq_length <- nchar(dna_sequence)
  grid_width <- if (is.null(num_columns)) ceiling(sqrt(seq_length)) else num_columns
  grid_height <- ceiling(seq_length / grid_width)
  
  padding <- grid_width * grid_height - seq_length
  nucleotides <- unlist(strsplit(dna_sequence, ""))
  nucleotides_padded <- c(nucleotides, rep(NA, padding))
  
  safe_accession <- gsub("[[:space:]/\\\\:*?\"<>|]", "_", accession)
  safe_name <- gsub("[[:space:]/\\\\:*?\"<>|]", "_", name)
  
  filename <- file.path(output_dir, paste0(safe_accession, "_", safe_name, "_plot.png"))
  nucleotide_matrix <- matrix(nucleotides_padded, nrow = grid_height, ncol = grid_width, byrow = TRUE)
  grid_df <- melt(nucleotide_matrix, varnames = c("row", "col"), value.name = "nucleotide")
  grid_df <- grid_df[!is.na(grid_df$nucleotide), ]
  nucleotide_colors <- c("A" = "red", "T" = "blue", "C" = "blue3", "G" = "red3", "N" = "grey")
  grid_df$nucleotide <- factor(grid_df$nucleotide, levels = c("A", "T", "C", "G", "N"))
  
  grid_df$row_offset <- grid_df$row
  grid_df$col_offset <- grid_df$col
  
  title_lines <- c(
    paste0("Query Accession: ", metadata[["Query"]]),
    paste0("Hit Accession: ", accession),
    paste0("Columns: ", grid_width, " | Rows: ", grid_height)
  )
  
  if (!is.null(metadata[["Strand"]])) {
    title_lines <- c(title_lines, paste0("Strand: ", metadata[["Strand"]]))
  }
  if (length(metadata) > 0) {
    meta_lines <- unlist(lapply(setdiff(names(metadata), c("Query", "Strand")), function(k) paste0(k, ": ", metadata[[k]])))
    title_lines <- c(title_lines, meta_lines)
  }
  
  highlight_thickness <- compute_highlight_thickness(grid_height, grid_width)
  highlight_layer <- {
    coords <- sort(highlight_range)
    coords <- coords[!is.na(coords)]
    if (length(coords) == 2 && all(coords > 0)) {
      highlight_idx <- coords[1]:coords[2]
      h_rows <- ((highlight_idx - 1) %/% grid_width) + 1
      h_cols <- ((highlight_idx - 1) %% grid_width) + 1
      h_df <- data.frame(row_offset = h_rows, col_offset = h_cols)
      geom_tile(data = h_df, aes(x = col_offset, y = row_offset), fill = NA, color = "black", linewidth = highlight_thickness)
    } else {
      NULL
    }
  }
  
  plot <- ggplot(grid_df, aes(x = col_offset, y = row_offset, fill = nucleotide)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(
      values = nucleotide_colors,
      na.value = "black",
      name = "Nucleotide") +
    theme(legend.position = "right") +
    highlight_layer +
    scale_y_reverse() +
    theme_minimal() +
    labs(title = paste(title_lines, collapse = "\n"), x = "Column No.", y = "Row No.") +
    coord_fixed()
  
  ggsave(filename, plot, width = 10, height = 10, bg = "white")
  cat("Saved plot:", filename, "\n")
}

# --- Main function ---

### run_blast_and_retrieve(input_seq,
#                          min_pid = 60,
#                          max_evalue = 1e-2,
#                          min_length = 20,
#                          coverage_threshold = 0.1,
#                          max_hits = 100,
#                          exclude_predicted = FALSE,
#                          user_defined_columns = 100,
#                          output_dir = NULL,
#                          save_fasta = TRUE,
#                          retry_delay = 3,
#                          max_attempts = 600,
#                          debug_reverse = TRUE,
#                          filter_refseq = FALSE,
#                          hitlist_size = 3000,
#                          enforce_hit_quota = TRUE)
#
# Purpose:
# Full wrapper function to run BLASTn against NCBI remotely, filter and retrieve hits,
# download matching sequences, and generate tile plots with alignment highlights.
#
# Behavior:
# - Determines if the input is an accession or raw sequence.
# - Runs a BLAST search remotely and polls until results are available.
# - Filters raw hits using identity, E-value, length, coverage, prediction status, and RefSeq filters.
# - Downloads matching full sequences.
# - Generates tile plots with alignment region highlighted, including reverse/complement debugging if enabled.
#
# Parameters:
# - input_seq: Accession number or raw DNA string.
# - min_pid: Minimum % identity to accept (default = 60).
# - max_evalue: Max E-value threshold (default = 1e-2).
# - min_length: Minimum alignment length (default = 20).
# - coverage_threshold: Min fraction of query covered (default = 0.1).
# - max_hits: Number of final hits to retain (after filtering).
# - exclude_predicted: If TRUE, removes predicted/hypothetical sequences.
# - user_defined_columns: Number of columns in the tile grid.
# - output_dir: Folder to store output files (auto-named if NULL).
# - save_fasta: Whether to save retrieved sequences as FASTA files.
# - retry_delay: Seconds between retry attempts for polling.
# - max_attempts: Max number of poll attempts for BLAST completion.
# - debug_reverse: If TRUE, generate RC and COMP debug plots for minus strand hits.
# - filter_refseq: If TRUE, only allow RefSeq accessions (e.g., containing underscore).
# - hitlist_size: Number of hits requested from NCBI BLAST (default = 3000).
# - enforce_hit_quota: If TRUE, keeps filtering until `max_hits` non-excluded hits are retained.
#
# Returns:
# - A list with metadata about the processed run and results.
#
# Output:
# - PNG tile plots with highlights.
# - FASTA files (optional).
# - Metadata files (log + .meta files per plot).
# - Optional debug plots (_RC and _COMP).

run_blast_and_retrieve <- function(input_seq,
                                   min_pid = 50,
                                   max_evalue = 1e-2,
                                   min_length = 20,
                                   coverage_threshold = 0.1,
                                   max_hits = 1000,
                                   exclude_predicted = FALSE,
                                   user_defined_columns = 100,
                                   output_dir = NULL,
                                   save_fasta = TRUE,
                                   retry_delay = 3,
                                   max_attempts = 600,
                                   debug_reverse = TRUE,
                                   filter_refseq = FALSE,
                                   hitlist_size = 3000) {
  
  seq_info <- detect_seq_type(input_seq)
  if (is.null(output_dir)) {
    name_clean <- gsub("[^A-Za-z0-9._-]", "_", seq_info$name)
    output_dir <- file.path("BLASTn_Plots", paste0(seq_info$accession, "_", name_clean))
  }
  if (seq_info$type != "nucleotide") return(cat("Not a valid nucleotide sequence or accession.\n"))
  
  xml <- run_remote_blast(seq_info, max_attempts, retry_delay, hitlist_size = hitlist_size)
  if (is.null(xml)) return(cat("BLAST failed or timed out.\n"))
  
  filter_result <- filter_hits(xml, seq_info, min_pid, max_evalue, min_length, coverage_threshold, max_hits, exclude_predicted, filter_refseq)
  hits <- filter_result$hits
  summary <- filter_result$summary
  
  cat("\nFiltering summary:\n")
  print(summary)
  
  # Write detailed metadata for retained hits
  log_file <- file.path(output_dir, "filter_log.txt")
  if (length(hits) > 0) {
    meta_lines <- c("\nMetadata of retained hits:")
    for (i in seq_along(hits)) {
      h <- hits[[i]]
      strand <- if (h$start > h$end) "minus" else "plus"
      meta_lines <- c(meta_lines, sprintf(
        "%d. Accession: %s | PID: %.1f%% | E-value: %g | Align len: %d | Range: %d-%d | Strand: %s",
        i, h$accession, h$pid, h$eval, h$len, h$start, h$end, strand
      ))
    }
    write(meta_lines, file = log_file, append = TRUE)
  }
  
  if (length(hits) == 0) return(cat("No hits passed filters.\n"))
  full_sequences <- retrieve_full_sequences(hits, save_fasta, output_dir)
  cat("Retrieved", length(full_sequences), "full sequences.\n")
  
  for (seq_entry in full_sequences) {
    strand <- if (seq_entry$start > seq_entry$end) "minus" else "plus"
    metadata <- list(
      "Query" = seq_info$accession,
      "PID" = sprintf("%.1f%%", seq_entry$pid),
      "E-value" = seq_entry$eval,
      "Align length" = seq_entry$len,
      "Alignment range" = paste0(seq_entry$start, "-", seq_entry$end),
      "Strand" = strand
    )
    
    create_dna_plot(
      sequence = seq_entry$sequence,
      header = seq_entry$header,
      name = seq_entry$name,
      accession = seq_entry$accession,
      num_columns = user_defined_columns,
      metadata = metadata,
      output_dir = output_dir,
      highlight_range = c(seq_entry$start, seq_entry$end)
    )
    
    if (debug_reverse && strand == "minus") {
      debug_reverse_plots(seq_entry, metadata, user_defined_columns, output_dir)
    }
  }
  
  return(full_sequences)
}

# --- Reverse Complement Utility ---
reverse_complement <- function(seq) {
  comp <- chartr("ATGC", "TACG", toupper(seq))
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

# --- Debug Reverse Block ---
debug_reverse_plots <- function(seq_entry, metadata, user_defined_columns, output_dir) {
  metadata$alignment_range <- paste0(seq_entry$start, "-", seq_entry$end, " of minus strand")
  aln_start <- min(seq_entry$start, seq_entry$end)
  aln_end <- max(seq_entry$start, seq_entry$end)
  aln_len <- abs(seq_entry$end - seq_entry$start) + 1
  rc_highlight <- c(nchar(seq_entry$sequence) - aln_end + 1, nchar(seq_entry$sequence) - aln_start + 1)
  
  create_dna_plot(
    sequence = reverse_complement(seq_entry$sequence),
    header = seq_entry$header,
    name = paste0(seq_entry$name, "_RC"),
    accession = seq_entry$accession,
    num_columns = user_defined_columns,
    metadata = metadata,
    output_dir = file.path(output_dir, "debug_reverse"),
    highlight_range = rc_highlight
  )
  
  create_dna_plot(
    sequence = chartr("ATGC", "TACG", toupper(substr(seq_entry$sequence, aln_start, aln_end))),
    header = seq_entry$header,
    name = paste0(seq_entry$name, "_COMP"),
    accession = seq_entry$accession,
    num_columns = user_defined_columns,
    metadata = metadata,
    output_dir = file.path(output_dir, "debug_complement"),
    highlight_range = c(1, aln_len)
  )
}


# Example: Run the full process with plotting
result <- run_blast_and_retrieve("NM_205518.2")

result <- run_blast_and_retrieve("NM_205518.2", max_hits = 100, filter_refseq = FALSE, exclude_predicted = TRUE, user_defined_columns = 150)

run_blast_and_retrieve("NM_000518.5")  # human beta-globin
