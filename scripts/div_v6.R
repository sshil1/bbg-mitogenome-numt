#!/usr/bin/env Rscript
# D-loop diversity/neutrality analysis for the 22 BBG animals.
# Input: /home/gmbluser/mt_remap/trees_v6/dloop_v2_aligned.fa (73 taxa: 22 BBG + refs)
# Output: bbg_diversity_summary.csv, bbg_mismatch_distribution.png

library(ape)
library(pegas)

aln_path <- "/home/gmbluser/mt_remap/trees_v6/dloop_v2_aligned.fa"
aln_full <- read.dna(aln_path, format = "fasta")

bbg_names <- c(
  "116f", "189f", "803m", "1buck_m_", "2183f", "246f",
  "2869f", "3237m", "369m", "4_Female", "4lb__f_", "584f",
  "7ps__f_", "814m", "836m", "849m", "855m", "9ks__f_",
  "Bk_6__m_", "Hir_13_m_", "Hir_5f", "Kh_5__m_"
)

missing <- setdiff(bbg_names, rownames(aln_full))
if (length(missing) > 0) {
  stop("These BBG names were not found in the alignment: ",
       paste(missing, collapse = ", "))
}

bbg <- aln_full[bbg_names, ]
cat("Subset dimensions (animals x sites):", nrow(bbg), "x", ncol(bbg), "\n")

pi <- nuc.div(bbg, pairwise.deletion = TRUE)
hd <- hap.div(bbg)

haps <- haplotype(bbg)
n_haplotypes <- nrow(haps)

seg <- seg.sites(bbg)
n_seg_sites <- length(seg)

taj <- tajima.test(bbg)

results <- data.frame(
  Statistic = c(
    "N_animals",
    "Alignment_length_bp",
    "Nucleotide_diversity_pi",
    "Haplotype_diversity_Hd",
    "N_haplotypes",
    "N_polymorphic_sites",
    "Tajimas_D",
    "Tajimas_D_Pvalue_beta",
    "Tajimas_D_Pvalue_normal"
  ),
  Value = c(
    nrow(bbg),
    ncol(bbg),
    pi,
    hd,
    n_haplotypes,
    n_seg_sites,
    taj$D,
    taj$Pval.beta,
    taj$Pval.normal
  )
)

write.csv(results, "bbg_diversity_summary.csv", row.names = FALSE)
cat("\nWrote bbg_diversity_summary.csv\n")
print(results)

# Mismatch distribution: pairwise raw distances (proportion of differing
# sites) converted to absolute nucleotide-difference counts.
d <- dist.dna(bbg, model = "raw", pairwise.deletion = TRUE)
mismatches <- as.vector(d) * ncol(bbg)

png("bbg_mismatch_distribution.png", width = 900, height = 700, res = 120)
hist(mismatches,
     breaks = seq(0, max(mismatches) + 1, by = 1),
     main = "Mismatch distribution - 22 BBG D-loop sequences",
     xlab = "Pairwise nucleotide differences",
     ylab = "Frequency",
     col = "steelblue", border = "white")
dev.off()
cat("Wrote bbg_mismatch_distribution.png\n")
