############################################################
# Project: Comparative Genomic Analysis of
#          Mycobacterium tuberculosis Strains
#
# Script: 01_download_data.R
# Source: NCBI Nucleotide Database
# Method: NCBI E-utilities API
# Requirements: Base R only
# Run from the repository root: Rscript code/01_download_data.R
############################################################

# 1. PROJECT SETTINGS
############################################################
options(timeout = 600)                       # 10 minute timeout
if (!dir.exists("data")) dir.create("data")  # create data folder

base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
number_of_genomes <- 5

fasta_file  <- "data/mtb_genomes.fasta"
ids_file    <- "data/mtb_ncbi_ids.csv"
failed_file <- "data/failed_downloads.txt"

# 2. SEARCH NCBI FOR GENOME RECORDS
############################################################
cat("MYCOBACTERIUM TUBERCULOSIS DATA DOWNLOAD\n")
cat("Step 1: Searching NCBI...\n")

search_term <- paste(
  'txid1773[Organism:exp]',
  'AND "complete genome"[Title]',
  'AND 3500000:5500000[SLEN]'
)

search_url <- paste0(
  base_url, "esearch.fcgi?db=nuccore&term=",
  URLencode(search_term, reserved = TRUE),
  "&retmax=", number_of_genomes, "&retmode=xml"
)

search_file <- tempfile(fileext = ".xml")
download.file(search_url, destfile = search_file, mode = "wb", method = "libcurl")
xml_text <- paste(readLines(search_file, warn = FALSE), collapse = "")
unlink(search_file)

if (grepl("<ERROR>", xml_text, fixed = TRUE)) {
  stop("NCBI returned an error during the search.")
}

# Extract NCBI record IDs
ids <- regmatches(xml_text, gregexpr("<Id>[0-9]+</Id>", xml_text))[[1]]
ids <- unique(gsub("</?Id>", "", ids))
if (length(ids) == 0) stop("No genome records found. Check the NCBI search.")
cat("Records found:", length(ids), "\n\n")

# 3. SAVE NCBI RECORD IDS
############################################################
write.csv(data.frame(nuccore_uid = ids), ids_file, row.names = FALSE)
cat("NCBI record IDs saved to:", ids_file, "\n\n")

# 4. PREPARE FASTA OUTPUT
############################################################
if (file.exists(fasta_file))  file.remove(fasta_file)
if (file.exists(failed_file)) file.remove(failed_file)

# 5. DOWNLOAD GENOMES ONE BY ONE
############################################################
cat("Step 2: Downloading genome sequences...\n\n")

successful_downloads <- 0
failed_ids <- character()

for (i in seq_along(ids)) {
  current_id <- ids[i]
  cat("Downloading genome", i, "of", length(ids), "| NCBI ID:", current_id, "\n")

  fasta_url <- paste0(
    base_url, "efetch.fcgi?db=nuccore&id=", current_id,
    "&rettype=fasta&retmode=text"
  )

  success <- FALSE

  # Try downloading up to 3 times
  for (attempt in 1:3) {
    temp_file <- tempfile(fileext = ".fasta")
    result <- try(
      download.file(fasta_url, destfile = temp_file, mode = "wb",
                    method = "libcurl", quiet = TRUE),
      silent = TRUE
    )

    if (!inherits(result, "try-error") &&
        file.exists(temp_file) &&
        file.info(temp_file)$size > 0) {

      fasta_lines <- try(readLines(temp_file, warn = FALSE), silent = TRUE)

      if (!inherits(fasta_lines, "try-error") &&
          length(fasta_lines) > 1 &&
          startsWith(fasta_lines[1], ">")) {

        sequence <- paste(fasta_lines[!startsWith(fasta_lines, ">")], collapse = "")

        # Check expected whole-genome length
        if (nchar(sequence) >= 3500000 && nchar(sequence) <= 5500000) {
          success <- TRUE
          break
        }
      }
    }

    unlink(temp_file)
    cat("  Attempt", attempt, "failed.\n")
    if (attempt < 3) Sys.sleep(5 * attempt)
  }

  # Save successfully downloaded FASTA
  if (success) {
    if (!file.exists(fasta_file)) {
      file.copy(temp_file, fasta_file)
    } else {
      file.append(fasta_file, temp_file)
    }
    successful_downloads <- successful_downloads + 1
    cat("  Download successful!\n")
  } else {
    failed_ids <- c(failed_ids, current_id)
    cat("  Download failed for ID:", current_id, "\n")
  }

  if (exists("temp_file") && file.exists(temp_file)) unlink(temp_file)
  Sys.sleep(0.5)   # be polite to NCBI servers
}

# 6. SAVE FAILED RECORD IDS
############################################################
if (length(failed_ids) > 0) {
  writeLines(failed_ids, failed_file)
  cat("\nSome records could not be downloaded.\n")
  cat("Failed IDs saved to:", failed_file, "\n")
}

# 7. VERIFY DOWNLOADED DATASET
############################################################
cat("\nDOWNLOAD SUMMARY\n")
if (file.exists(fasta_file)) {
  final_lines <- readLines(fasta_file, warn = FALSE)
  cat("Successful downloads:", successful_downloads, "\n")
  cat("FASTA records:", sum(startsWith(final_lines, ">")), "\n")
  cat("File size (MB):", round(file.info(fasta_file)$size / (1024^2), 2), "\n")
} else {
  cat("No FASTA file was created.\n")
}
cat("Failed downloads:", length(failed_ids), "\n")
cat("Download process finished.\n")
