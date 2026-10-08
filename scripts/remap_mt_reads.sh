#!/usr/bin/env bash
# remap_mt_reads.sh SAMPLE [SAMPLE ...]
# NUMT-aware mitochondrial read recovery. Read-only on /data/bam and the reference; writes only under ~/mt_remap/.
#  0) (once) find NUMT loci by aligning NC_005044.2 to the nuclear genome with minimap2 -> numt_regions.bed
#  1) pull every read overlapping chrM or a NUMT locus from /data/bam/sorted/SAMPLE.sorted.bam
#  2) re-map those reads to the mitochondrial reference ALONE (so mtDNA and NUMT reads compete only with each other;
#     mtDNA is far more abundant, so the majority allele is the true mt allele)
#  3) output: ~/mt_remap/SAMPLE/SAMPLE_mt.bam  (name chosen so majority_consensus.py labels it SAMPLE)
set -euo pipefail
REF=${REF:-/data/db/reference/goat_ARS1.fna}
BAMDIR=${BAMDIR:-/data/bam/sorted}
W=$HOME/mt_remap; THREADS=${THREADS:-8}
mkdir -p "$W"
MT=$W/ref_mt.fa
[ -s "$MT" ] || { samtools faidx "$REF" NC_005044.2 > "$MT"; samtools faidx "$MT"; }
BED=$W/numt_regions.bed
if [ ! -s "$BED" ]; then
  echo "[0] locating NUMT loci (one-off; minimap2 index of the nuclear genome takes a few minutes)"
  minimap2 -x asm20 -t "$THREADS" -N 200 -p 0.1 -c "$REF" "$MT" > "$W/mt_vs_genome.paf" 2> "$W/minimap2_numt.log"
  awk -F'\t' '$6!="NC_005044.2" && $11>=100 && $10/$11>=0.80 { s=$8-300; if(s<0)s=0; print $6"\t"s"\t"($9+300) }' "$W/mt_vs_genome.paf" \
   | sort -k1,1 -k2,2n \
   | awk 'BEGIN{OFS="\t"} NR==1{c=$1;s=$2;e=$3;next} $1==c && $2<=e {if($3>e)e=$3; next} {print c,s,e; c=$1;s=$2;e=$3} END{if(NR)print c,s,e}' > "$W/numt_nuclear.bed"
  { printf "NC_005044.2\t0\t16643\n"; cat "$W/numt_nuclear.bed"; } > "$BED"
  echo "    NUMT loci: $(wc -l < "$W/numt_nuclear.bed")  (regions file: $BED)"
fi
for S in "$@"; do
  BAM=$BAMDIR/$S.sorted.bam; D=$W/$S; mkdir -p "$D"
  echo "[1] $S: extracting reads from chrM + NUMT loci"
  samtools view -@ "$THREADS" -b -M -L "$BED" "$BAM" > "$D/$S.region.bam"
  samtools fastq -@ "$THREADS" -F 0x900 -n "$D/$S.region.bam" > "$D/$S.reads.fq" 2> "$D/fastq.log"
  echo "    reads recovered: $(( $(wc -l < "$D/$S.reads.fq") / 4 ))"
  echo "[2] $S: re-mapping to the mitochondrial reference only"
  minimap2 -ax sr -t "$THREADS" "$MT" "$D/$S.reads.fq" 2> "$D/minimap2_remap.log" \
   | samtools sort -@ 4 -o "$D/${S}_mt.bam" -
  samtools index "$D/${S}_mt.bam"
  echo "    done: $D/${S}_mt.bam   (mean depth: $(samtools depth -a "$D/${S}_mt.bam" | awk '{s+=$3} END{printf "%.0f", s/NR}'))"
  rm -f "$D/$S.region.bam" "$D/$S.reads.fq"
done
