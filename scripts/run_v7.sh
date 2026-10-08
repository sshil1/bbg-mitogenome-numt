#!/usr/bin/env bash
set -uo pipefail
export MAFFT=/opt/conda/bin/mafft RSCRIPT=/opt/conda/bin/Rscript
IQ=/opt/conda/envs/variants/bin/iqtree3; T=8; cd ~
GT=$HOME/mt_remap/mitos_v2/803m/803m_gene_table_v2_new.csv
echo "### 1 derive v7 scripts"
for f in rebuild_v6.sh run_popgen_dloop_v6.sh popgen_v6.R compare_world_v6.R nearest_refs_v6.R; do sed 's/v6/v7/g' ~/$f > ~/${f/v6/v7}; done
awk '1; /^PY$/{exit}' ~/rebuild_v7.sh > ~/rebuild_v7a.sh
echo "### 2 raw inputs"; bash ~/rebuild_v7a.sh || { echo FAIL_RAW; exit 1; }
cd ~/mt_remap/trees_v7
awk '/^>/{skip=($1==">KP677508.1")} !skip' world_v2_raw.fa > world_nomoth_raw.fa
echo "world_nomoth taxa: $(grep -c '>' world_nomoth_raw.fa) (expect 101)"
for s in complete dloop; do $MAFFT --auto --thread $T ${s}_v2_raw.fa > ${s}_v2_aligned.fa 2> ${s}_mafft.log; done
$MAFFT --auto --thread $T world_nomoth_raw.fa > world_nomoth_aligned.fa 2> world_mafft.log
echo "### 3 trees"
$IQ -s complete_v2_aligned.fa -m MFP -B 1000 --alrt 1000 -T $T --prefix complete_v2 > /dev/null
$IQ -s dloop_v2_aligned.fa -m MFP -B 1000 --alrt 1000 -T $T --prefix dloop_v2 > /dev/null
$IQ -s world_nomoth_aligned.fa -m MFP -B 1000 --alrt 1000 -T $T --prefix world_nomoth > /dev/null
for p in complete_v2 dloop_v2 world_nomoth; do echo "== $p"; grep -E "Best-fit model|Total tree length" $p.iqtree; done
echo "### 4 popgen / world / nearest refs"
RSCRIPT=$RSCRIPT bash ~/run_popgen_dloop_v7.sh > ~/mt_remap/results_v7_dloop.log 2>&1
$RSCRIPT ~/popgen_v7.R > ~/mt_remap/popgen_v7.log 2>&1
$RSCRIPT ~/compare_world_v7.R > ~/mt_remap/compare_world_v7.log 2>&1
$RSCRIPT ~/nearest_refs_v7.R > ~/mt_remap/nearest_refs_v7.log 2>&1
echo "### 5 dN/dS and codons"
D=~/mt_remap/results_v7; mkdir -p $D/dnds
python3 ~/dnds_analysis.py $GT ~/mt_remap/fasta_v7/ $D/dnds/dnds_v7.csv > $D/dnds_v7.log 2>&1
python3 ~/dnds_analysis.py $GT ~/mt_remap/fasta_v7core/ $D/dnds/dnds_v7core.csv > $D/dnds_v7core.log 2>&1
python3 ~/variant_codons.py ~/mt_remap/bbg_22S_mt_consensus_v7.fasta $GT 2 > $D/codons_v7.txt 2>&1
python3 ~/variant_codons.py ~/mt_remap/bbg_22S_mt_consensus_v7core.fasta $GT 2 > $D/codons_v7core.txt 2>&1
echo "### DONE"
