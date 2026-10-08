suppressMessages(library(ape))
bbg <- c("116f","189f","1buck_m_","2183f","246f","2869f","3237m","369m","4_Female","4lb__f_","584f","7ps__f_","803m","814m","836m","849m","855m","9ks__f_","Bk_6__m_","Hir_13_m_","Hir_5f","Kh_5__m_")
for (f in c("/data/mt/trees/bbg_world_ML.treefile","/home/gmbluser/mt_remap/trees_v6/world_nomoth.treefile")) {
  t <- read.tree(f); cat("\n##",f,"| tips:",Ntip(t),"\n")
  b <- intersect(bbg,t$tip.label); cat("BBG tips found:",length(b),"\n")
  cat("BBG monophyletic:", is.monophyletic(t,b),"\n")
  m <- getMRCA(t,b); inside <- setdiff(extract.clade(t,m)$tip.label,bbg)
  cat("Non-BBG tips inside the BBG MRCA clade:",length(inside),"\n"); print(inside)
  d <- cophenetic(t)[b, setdiff(t$tip.label,bbg), drop=FALSE]
  cat("Closest references (min patristic distance):\n"); print(round(head(sort(apply(d,2,min)),6),5))
}
