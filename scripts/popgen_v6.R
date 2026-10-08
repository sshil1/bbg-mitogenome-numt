suppressMessages({library(ape);library(pegas)})
st <- function(m,label){
  x2 <- as.DNAbin(m); keep <- apply(m,2,function(c) all(c!="N")); mc <- as.DNAbin(m[,keep,drop=FALSE])
  td <- tajima.test(mc)
  cat(sprintf("%-10s sites=%5d complete=%5d S=%3d haplotypes=%2d Hd=%.4f pi=%.5f TajimaD=%.3f (p=%.3f)\n",
    label, ncol(m), sum(keep), length(seg.sites(mc)), nrow(haplotype(mc)), hap.div(mc),
    nuc.div(x2, pairwise.deletion=TRUE), td$D, td$Pval.normal))
}
for (f in c("/data/mt/bbg_22S_mt_consensus.fasta","/home/gmbluser/mt_remap/bbg_22S_mt_consensus_v6.fasta")) {
  cat("\n##", basename(f), "\n")
  m <- toupper(as.character(read.dna(f,"fasta"))); m[!(m %in% c("A","C","G","T"))] <- "N"
  st(m,"full mt"); st(m[,15434:16643],"D-loop"); st(m[,1:15433],"non-D-loop")
}
