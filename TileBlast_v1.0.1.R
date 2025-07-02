library(httr)
library(XML)
library(rentrez)
library(ggplot2)
library(reshape2)

# --- Detect sequence type ---

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

complement_only <- function(seq) {
  chartr("ATGC", "TACG", toupper(seq))
}


# --- Run BLAST remotely and get XML ---

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
      if (!is.null(xml)) {
        raw_hits <- getNodeSet(xml, "//Hit")
        cat("Total raw hits returned by BLAST:", length(raw_hits), "
")
        break
      }
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
    hit_len <- as.numeric(xmlValue(hit[["Hit_len"]]))
    hstart <- as.numeric(xmlValue(hsp[["Hsp_hit-from"]]))
    hend <- as.numeric(xmlValue(hsp[["Hsp_hit-to"]]))
    cov <- len / nchar(seq_info$seq)
    def <- xmlValue(hit[["Hit_def"]])
    acc <- xmlValue(hit[["Hit_accession"]])
    
    reason <- ""
    
    if (filter_refseq && !grepl("^(NM_|NR_|NG_)", acc)) {
      next
    }
    if (exclude_predicted && (
      grepl("(PREDICTED|unnamed|hypothetical|LOW QUALITY)", def, ignore.case = TRUE) ||
      grepl("^(XM_|XR_|XP_|YP_|ZP_)", acc)
    )) {
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
        len = len,
        hit_length = hit_len
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


reverse_complement <- function(seq) {
  comp <- chartr("ATGC", "TACG", toupper(seq))
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

# --- Compute Highlight Thickness ---

compute_highlight_thickness <- function(num_rows, num_cols) {
  total_dim <- num_rows + num_cols
  ratio <- max(num_rows / num_cols, num_cols / num_rows)
  falloff <- exp(-0.001 * total_dim) * (1 / (1 + (ratio - 1)^1.5))
  thickness <- 0.9 * falloff
  return(min(max(thickness, 0.05), 0.9))
}


# --- Helper to generate y-axis labels dynamically ---

generate_y_axis_labels <- function(grid_height, desired_tick_count = 10, label_x = 0, label_size = 3) {
  tick_step <- ceiling(grid_height / desired_tick_count)
  y_vals <- which(seq_len(grid_height) %% tick_step == 0)
  geom_text(
    data = data.frame(y = y_vals, x = label_x),
    aes(x = x, y = grid_height - y + 1, label = y),
    inherit.aes = FALSE,
    hjust = 1,
    size = label_size
  )
}

# --- Modified create_dna_plot with label-only inversion ---

create_dna_plot <- function(sequence, header, name, accession, num_columns = NULL, metadata = list(), output_dir = "BLASTn_Plots", highlight_range = NULL, desired_tick_count = 20) {
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
  grid_df <- reshape2::melt(nucleotide_matrix, varnames = c("row", "col"), value.name = "nucleotide")
  grid_df <- grid_df[!is.na(grid_df$nucleotide), ]
  grid_df$nucleotide <- factor(grid_df$nucleotide, levels = c("A", "T", "C", "G", "N"))
  
  grid_df$row_offset <- grid_height - grid_df$row + 1
  grid_df$col_offset <- grid_df$col
  
  title_lines <- c(
    paste0("Query Accession: ", metadata[["Query"]], " (", nchar(sequence), " bp)"),
    paste0("Hit Accession: ", accession, " (", if (!is.null(metadata[["Hit length"]])) metadata[["Hit length"]] else "?", " bp)"),
    paste0("Columns: ", grid_width, " | Rows: ", grid_height)
  )
  if (!is.null(metadata[["Strand"]])) {
    title_lines <- c(title_lines, paste0("Strand: ", metadata[["Strand"]]))
  }
  if (length(metadata) > 0) {
    meta_lines <- unlist(lapply(setdiff(names(metadata), c("Query", "Strand", "Hit length")), function(k) paste0(k, ": ", metadata[[k]])))
    title_lines <- c(title_lines, meta_lines)
  }
  
  highlight_thickness <- compute_highlight_thickness(grid_height, grid_width)
  highlight_layer <- {
    coords <- sort(highlight_range)
    coords <- coords[!is.na(coords)]
    if (length(coords) == 2 && all(coords > 0)) {
      highlight_idx <- coords[1]:coords[2]
      highlight_idx <- highlight_idx[highlight_idx <= seq_length]  # avoid padding
      h_rows <- ((highlight_idx - 1) %/% grid_width) + 1
      h_cols <- ((highlight_idx - 1) %% grid_width) + 1
      h_row_offset <- grid_height - h_rows + 1  # correct for 5' to 3' layout
      h_df <- data.frame(row_offset = h_row_offset, col_offset = h_cols)
      geom_tile(data = h_df, aes(x = col_offset, y = row_offset), fill = NA, color = "black", linewidth = highlight_thickness)
    } else {
      NULL
    }
  }
  
  nucleotide_colors <- c("A" = "red", "T" = "blue", "C" = "blue3", "G" = "red3", "N" = "grey")
  
  plot <- ggplot(grid_df, aes(x = col_offset, y = row_offset, fill = nucleotide)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(values = nucleotide_colors, na.value = "black", name = "Nucleotide") +
    NULL +
    highlight_layer +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(title = paste(title_lines, collapse = "\n"), x = "Column No.", y = "Row No.") +
    generate_y_axis_labels(grid_height, desired_tick_count) +
    coord_fixed()
  
  ggsave(filename, plot, width = 10, height = 10, bg = "white")
  cat("Saved plot:", filename, "\n")
}

# --- Main function ---

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
                                   hitlist_size = 3000,
                                   desired_tick_count = 20, slice_start = NULL, slice_end = NULL) {
  
  seq_info <- detect_seq_type(input_seq)
  if (!is.null(slice_start) || !is.null(slice_end)) {
    if (!is.null(slice_start) && (!is.numeric(slice_start) || slice_start < 1)) stop("Invalid slice_start")
    if (!is.null(slice_end) && (!is.numeric(slice_end) || slice_end < 1)) stop("Invalid slice_end")
    slice_start <- ifelse(is.null(slice_start), 1, slice_start)
    slice_end <- ifelse(is.null(slice_end), nchar(seq_info$seq), slice_end)
    if (slice_end > nchar(seq_info$seq)) warning("slice_end exceeds sequence length; will be capped")
    seq_info$seq <- substr(seq_info$seq, slice_start, min(slice_end, nchar(seq_info$seq)))
    seq_info$name <- paste0(seq_info$name, "_slice", slice_start, "to", slice_end)
  }
  if (is.null(output_dir)) {
    name_clean <- gsub("[^A-Za-z0-9._-]", "_", seq_info$name)
    output_dir <- file.path("BLASTn_Plots", paste0(seq_info$accession, "_", name_clean))
  }
  if (seq_info$type != "nucleotide") return(cat("Not a valid nucleotide sequence or accession.
"))
  
  xml <- run_remote_blast(seq_info, max_attempts, retry_delay, hitlist_size = hitlist_size)
  if (is.null(xml)) return(cat("BLAST failed or timed out.
"))
  
  filter_result <- filter_hits(xml, seq_info, min_pid, max_evalue, min_length, coverage_threshold, max_hits, exclude_predicted, filter_refseq)
  hits <- filter_result$hits
  summary <- filter_result$summary
  
  cat("
Filtering summary:
")
  print(summary)
  
  if (length(hits) == 0) return(cat("No hits passed filters.
"))
  
  full_sequences <- retrieve_full_sequences(hits, save_fasta, output_dir)
  cat("Retrieved", length(full_sequences), "full sequences.
")
  
  # Write metadata log with updated accessions
  log_file <- file.path(output_dir, "filter_log.txt")
  if (length(full_sequences) > 0) {
    meta_lines <- c("
Metadata of retained hits:")
    for (i in seq_along(full_sequences)) {
      h <- full_sequences[[i]]
      strand <- if (h$start > h$end) "minus" else "plus"
      meta_lines <- c(meta_lines, sprintf(
        "%d. Accession: %s | PID: %.1f%% | E-value: %g | Align len: %d | Range: %d-%d | Strand: %s",
        i, h$accession, h$pid, h$eval, h$len, h$start, h$end, strand
      ))
    }
    write(meta_lines, file = log_file, append = TRUE)
  }
  
  for (seq_entry in full_sequences) {
    strand <- if (seq_entry$start > seq_entry$end) "minus" else "plus"
    metadata <- list(
      "Query" = if (!is.null(slice_start) || !is.null(slice_end)) paste0(seq_info$accession, " (", slice_start, ":", slice_end, ")") else seq_info$accession,
      "PID" = sprintf("%.1f%%", seq_entry$pid),
      "E-value" = seq_entry$eval,
      "Align length" = seq_entry$len,
      "Alignment range" = paste0(seq_entry$start, "-", seq_entry$end),
      "Strand" = strand,
      "Hit length" = seq_entry$hit_length
    )
    
    create_dna_plot(
      sequence = seq_entry$sequence,
      header = seq_entry$header,
      name = seq_entry$name,
      accession = seq_entry$accession,
      num_columns = user_defined_columns,
      metadata = metadata,
      output_dir = output_dir,
      highlight_range = c(seq_entry$start, seq_entry$end),
      desired_tick_count = desired_tick_count
    )
    
    if (debug_reverse && strand == "minus") {
      debug_reverse_plots(seq_entry, metadata, user_defined_columns, output_dir, desired_tick_count)
    }
  }
  
  return(full_sequences)
}

# --- Reverse Complement Utility ---
reverse_complement <- function(seq) {
  comp <- chartr("ATGC", "TACG", toupper(seq))
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

debug_reverse_plots <- function(seq_entry, metadata, user_defined_columns, output_dir, desired_tick_count = 20) {
  metadata$alignment_range <- paste0(seq_entry$start, "-", seq_entry$end, " of minus strand")
  
  # Extract the aligned region (always low to high)
  aln_start <- min(seq_entry$start, seq_entry$end)
  aln_end <- max(seq_entry$start, seq_entry$end)
  aln_len <- aln_end - aln_start + 1
  
  # Reverse complement the aligned region so it matches the query visually
  aligned_substring <- substr(seq_entry$sequence, aln_start, aln_end)
  rc_seq <- reverse_complement(aligned_substring)
  
  rc_metadata <- metadata
  rc_metadata[["Strand"]] <- "reverse complement of aligned region (5′→3′ visual match to query)"
  
  create_dna_plot(
    sequence = rc_seq,
    header = seq_entry$header,
    name = paste0(seq_entry$name, "_RC"),
    accession = seq_entry$accession,
    num_columns = user_defined_columns,
    metadata = rc_metadata,
    output_dir = file.path(output_dir, "debug_reverse_complement"),
    highlight_range = c(1, aln_len),
    desired_tick_count = desired_tick_count
  )
}


# Example: Run the full process with plotting
#result <- run_blast_and_retrieve("NM_205518.2")

#result <- run_blast_and_retrieve("NM_001257893.1", max_hits = 10, filter_refseq = FALSE, exclude_predicted = TRUE, user_defined_columns = 100)

#run_blast_and_retrieve("NM_000518.5")  # human beta-globin


result <- run_blast_and_retrieve("NM_205518.2",
                                 slice_start = NULL,
                                 slice_end = NULL,
                                 max_hits = 100, 
                                 filter_refseq = FALSE, 
                                 exclude_predicted = TRUE, 
                                 user_defined_columns = 70, 
                                 desired_tick_count = 10,  
                                 hitlist_size = 1000)



