############################################################
# Script: 02_import_qc_eda.R
# Project: Comparative Genomic Analysis of M. tuberculosis Strains
# Purpose: (1) import data, (2) quality control and cleaning,
#          (3) exploratory data analysis (Figures 1-3)
# Run from the repository root: Rscript code/02_import_qc_eda.R
############################################################

# 1. PACKAGES
############################################################
# Install any missing packages (Biostrings comes from Bioconductor)
for (p in c("ggplot2", "pheatmap")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(Biostrings)  # DNA sequence handling
  library(ggplot2)     # plots
  library(pheatmap)    # heatmaps
})

# Output folders
for (d in c("figures", "results")) if (!dir.exists(d)) dir.create(d)

# 2. DATA IMPORT
############################################################
fasta_file <- "data/mtb_genomes.fasta"
if (!file.exists(fasta_file)) {
  stop("FASTA file not found. Run code/01_download_data.R first.")
}

genomes <- readDNAStringSet(fasta_file)

# Extract accession and strain name from the FASTA header, e.g.
# ">NZ_CP154696.1 Mycobacterium tuberculosis strain CT2018-00088 chromosome, complete genome"
headers   <- names(genomes)
accession <- sub(" .*", "", headers)
strain    <- sub(".* strain (.*) chromosome.*", "\\1", headers)
names(genomes) <- strain

cat("Genomes imported:", length(genomes), "\n")

# 3. QUALITY CONTROL AND CLEANING
############################################################
# Count A, C, G, T and "other" (ambiguous) bases per genome
base_counts <- alphabetFrequency(genomes, baseOnly = TRUE)

qc <- data.frame(
  accession = accession,
  strain    = strain,
  length_bp = width(genomes),
  A = base_counts[, "A"],
  C = base_counts[, "C"],
  G = base_counts[, "G"],
  T = base_counts[, "T"],
  other = base_counts[, "other"],
  row.names = NULL
)
qc$gc_percent    <- (qc$G + qc$C) / qc$length_bp * 100
qc$other_percent <- qc$other / qc$length_bp * 100

# QC checks
cat("Duplicated accessions:", sum(duplicated(accession)), "\n")
cat("Duplicated sequences :", sum(duplicated(as.character(genomes))), "\n")
cat("Ambiguous bases (N etc.):", sum(qc$other), "\n")

# Cleaning: keep complete-genome-sized, unambiguous, non-duplicate records
keep <- qc$length_bp >= 3.5e6 & qc$length_bp <= 5.5e6 &
        qc$other_percent < 1 & !duplicated(accession)
cat("Genomes removed by QC:", sum(!keep), "\n")
genomes <- genomes[keep]
qc      <- qc[keep, ]

# Save QC table and print summary statistics
write.csv(qc[, c("accession", "strain", "length_bp", "gc_percent", "other")],
          "results/qc_summary.csv", row.names = FALSE)
print(qc[, c("accession", "strain", "length_bp", "gc_percent", "other")])
cat("\nGenome length summary (bp):\n"); print(summary(qc$length_bp))
cat("\nGC content summary (%):\n");     print(summary(qc$gc_percent))

# 4. EXPLORATORY DATA ANALYSIS
############################################################
# Figure 1: genome length per strain
p1 <- ggplot(qc, aes(x = reorder(strain, length_bp), y = length_bp / 1e6)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = format(length_bp, big.mark = ",")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Genome length of M. tuberculosis strains",
       x = "Strain", y = "Genome length (Mb)") +
  theme_minimal()
ggsave("figures/fig1_genome_length.png", p1, width = 7, height = 4, dpi = 300)

# Figure 2: GC content per strain (y-axis zoomed because differences are tiny)
p2 <- ggplot(qc, aes(x = strain, y = gc_percent)) +
  geom_point(size = 4, colour = "firebrick") +
  scale_y_continuous(limits = c(65.5, 65.75)) +
  labs(title = "GC content of M. tuberculosis strains",
       x = "Strain", y = "GC content (%)") +
  theme_minimal()
ggsave("figures/fig2_gc_content.png", p2, width = 7, height = 4, dpi = 300)

# Figure 3: dinucleotide frequencies (rows scaled so small differences show)
dinuc <- dinucleotideFrequency(genomes, as.prob = TRUE) * 100
rownames(dinuc) <- names(genomes)
write.csv(round(dinuc, 4), "results/dinucleotide_frequencies.csv")
pheatmap(t(dinuc), scale = "row",
         main = "Dinucleotide frequencies (row-scaled)",
         filename = "figures/fig3_dinucleotide_heatmap.png",
         width = 7, height = 6)

cat("\nEDA finished. Figures saved in figures/\n")
