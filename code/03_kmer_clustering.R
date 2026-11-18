############################################################
# Script: 03_kmer_clustering.R
# Project: Comparative Genomic Analysis of M. tuberculosis Strains
# Purpose: Main analysis - alignment-free comparison of strains using
#          4-mer frequency profiles, hierarchical clustering (Figures 4-5)
# Run from the repository root: Rscript code/03_kmer_clustering.R
############################################################

# 1. PACKAGES
############################################################
suppressPackageStartupMessages({
  library(Biostrings)
  library(pheatmap)
})
for (d in c("figures", "results")) if (!dir.exists(d)) dir.create(d)

# 2. IMPORT AND CLEAN (same rules as 02_import_qc_eda.R)
############################################################
genomes <- readDNAStringSet("data/mtb_genomes.fasta")
accession <- sub(" .*", "", names(genomes))
names(genomes) <- sub(".* strain (.*) chromosome.*", "\\1", names(genomes))

other_pct <- alphabetFrequency(genomes, baseOnly = TRUE)[, "other"] / width(genomes) * 100
keep <- width(genomes) >= 3.5e6 & width(genomes) <= 5.5e6 &
        other_pct < 1 & !duplicated(accession)
genomes <- genomes[keep]
cat("Genomes used:", length(genomes), "\n")

# 3. K-MER FREQUENCY PROFILES
############################################################
# Relative frequency of each of the 256 possible 4-mers per genome
kmer_freq <- oligonucleotideFrequency(genomes, width = 4, as.prob = TRUE)
rownames(kmer_freq) <- names(genomes)
cat("k-mer matrix dimensions:", dim(kmer_freq), "\n")

# 4. DISTANCE MATRIX AND HIERARCHICAL CLUSTERING
############################################################
dist_mat <- dist(kmer_freq, method = "euclidean")
write.csv(as.matrix(dist_mat), "results/kmer4_distance_matrix.csv")
print(signif(as.matrix(dist_mat), 3))

hc <- hclust(dist_mat, method = "average")   # UPGMA
cat("\nMerge heights:\n"); print(signif(hc$height, 3))

# Figure 4: dendrogram
png("figures/fig4_kmer_dendrogram.png", width = 2000, height = 1400, res = 300)
plot(hc, hang = -1, xlab = "", sub = "",
     ylab = "Euclidean distance (4-mer frequencies)",
     main = "Clustering of M. tuberculosis strains")
dev.off()

# Figure 5: distance heatmap
pheatmap(as.matrix(dist_mat),
         clustering_distance_rows = dist_mat,
         clustering_distance_cols = dist_mat,
         clustering_method = "average",
         display_numbers = TRUE, number_format = "%.2e",
         main = "Pairwise 4-mer distances",
         filename = "figures/fig5_distance_heatmap.png",
         width = 7, height = 6)

# 5. SESSION INFORMATION (for reproducibility)
############################################################
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("\nAnalysis finished. Figures saved in figures/\n")
