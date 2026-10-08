#!/usr/bin/env bash
# remap_rot.sh SAMPLE [SAMPLE ...]
# Origin-edge-safe version of remap_mt_reads.sh: identical read recovery (chrM + NUMT loci from the original WGS BAM),
# but reads are re-mapped to a ROTATED copy of NC_005044.2 so that the linear-reference ends (positions 1-300 and the
# end of the control region) lie in the middle of the contig and keep full mitochondrial coverage.
# Requires remap_mt_reads.sh to have been run once (needs ~/mt_remap/ref_mt.fa and ~/mt_remap/numt_regions.bed).
set -euo pipefail
SHIFT=${SHIFT:-8000}; THREADS=${THREADS:-8}; BAMDIR=${BAMDIR:-/data/bam/sorted}
W=$HOME/mt_remap; R=$HOME/mt_remap_rot; mkdir -p "$R"
BED=$W/numt_regions.bed; MT0=$W/ref_mt.fa; ROT=$R/ref_mt_rot.fa
[ -s "$BED" ] && [ -s "$MT0" ] || { echo "need $BED and $MT0 (run remap_mt_reads.sh first)"; exit 1; }
if [ ! -s "$ROT" ]; then
  python3 - "$MT0" "$ROT" "$SHIFT" <<'PY'
import sys
src,dst,sh=sys.argv[1],sys.argv[2],int(sys.argv[3])
s=''.join(l.strip() for l in open(src) if not l.startswith('>')).upper()
r=s[sh:]+s[:sh]
open(dst,'w').write('>NC_005044.2\n'+'\n'.join(r[i:i+80] for i in range(0,len(r),80))+'\n')
print('rotated reference written:',len(r),'bp, shift',sh)
PY
  samtools faidx "$ROT"
fi
for S in "$@"; do
  BAM=$BAMDIR/$S.sorted.bam; D=$R/$S; mkdir -p "$D"
  echo "[1] $S: extracting reads from chrM + NUMT loci"
  samtools view -@ "$THREADS" -b -M -L "$BED" "$BAM" > "$D/$S.region.bam"
  samtools fastq -@ "$THREADS" -F 0x900 -n "$D/$S.region.bam" > "$D/$S.reads.fq" 2> "$D/fastq.log"
  echo "    reads recovered: $(( $(wc -l < "$D/$S.reads.fq") / 4 ))"
  echo "[2] $S: re-mapping to the ROTATED mitochondrial reference"
  minimap2 -ax sr -t "$THREADS" "$ROT" "$D/$S.reads.fq" 2> "$D/minimap2_remap.log" | samtools sort -@ 4 -o "$D/${S}_mt.bam" -
  samtools index "$D/${S}_mt.bam"
  echo "    done: $D/${S}_mt.bam (mean depth: $(samtools depth -a "$D/${S}_mt.bam" | awk '{s+=$3} END{printf "%.0f", s/NR}'))"
  rm -f "$D/$S.region.bam" "$D/$S.reads.fq"
done
