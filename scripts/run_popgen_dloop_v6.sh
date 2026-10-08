set -u
RS=${RSCRIPT:?set RSCRIPT}
R=$HOME/mt_remap/results_v6; mkdir -p $R/old $R/v2
NEW=$HOME/mt_remap/trees_v6/dloop_v2_aligned.fa
OLD=/data/mt/alignments/bbg_dloop_aligned.fasta
for tag in old v2; do
  if [ $tag = old ]; then A=$OLD; else A=$NEW; fi
  cd $R/$tag
  sed -e "s#/data/mt/alignments/bbg_dloop_aligned.fasta#$A#" \
      -e 's/pi <- nuc.div(bbg)$/pi <- nuc.div(bbg, pairwise.deletion = TRUE)/' \
      $HOME/diversity_neutrality.R | sed 's/\r$//' > div.R
  echo "##### $tag: nuc.div lines"; grep -n "nuc.div" div.R
  echo "##### $tag: diversity_neutrality"; $RS div.R 2>&1 | tail -40
  echo "##### $tag: fus_fs"; $RS $HOME/fus_fs_FINAL.R $A 2>&1 | tail -40
done
echo ALL_DONE
