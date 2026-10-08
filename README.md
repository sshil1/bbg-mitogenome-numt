# Black Bengal goat mitogenomes: NUMT-aware reconstruction (v7)

Scripts, final consensus sequences and results behind the manuscript
"NUMT-aware reconstruction of 22 Black Bengal goat (*Capra hircus*) mitogenomes places the breed
within haplogroup A" (target: *Mitochondrial DNA Part A*).

Archive: GitHub `sshil1/bbg-mitogenome-numt`; Zenodo DOI: 10.5281/zenodo.23236016. GenBank accessions for the 22 consensus mitogenomes will be added on acceptance.

## What is here

| Folder | Contents |
|---|---|
| `data/consensus_v7/` | **Final** consensus: 22 sequences, 16,643 bp, 19 IUPAC codes, 205 N (`*v7core.fasta`: positions 4, 14526 and 16495 set to N) |
| `data/per_animal_v7/` | The same sequences, one FASTA per animal (names follow the sampling IDs) |
| `data/per_animal_v7core/` | As above with positions 4, 14526 and 16495 set to N (used for variant and dN/dS analyses) |
| `data/provenance/` | Earlier consensus versions (v2: first NUMT-aware remap; v5: after majority resolution; v6: 11 curated positions) |
| `data/hifi/` | The two PacBio HiFi mitogenomes: 803M (GenBank PZ809991; used as the allele reference) and D863F (GenBank PZ809992). GenBank records are released on publication; see `HIFI_PROVENANCE.txt` |
| `results/trees/` | IQ-TREE 3.1.3 outputs for complete (56 taxa), D-loop (73) and world (101) sets |
| `results/alignments/` | MAFFT alignments after rotating all references to the NC_005044.2 origin |
| `results/diversity/` | pegas summary for v7 and for the conventional consensus; `popgen_v7.log`; Fu's Fs log |
| `results/selection/` | dN/dS tables (`*_v7core.csv` is the reported one), codon counts and logs |
| `results/network/` | PopART-ready inputs (`bbg_dloop_popart_DATA_block.nex`, a CRLF variant, FASTA, PHYLIP, `name_map.tsv`) and the original NEXUS inputs for PopART (strict D-loop columns) and haplotype membership; file names keep `v6` because the D-loop columns are identical in v6 and v7 (positions 203 and 211 lie outside the D-loop) |
| `results/curation/` | `curated_positions_v7.tsv` (Supplementary Table S1: 13 curated positions and 3 masked positions, calls at each stage) and per-position review tables |
| `results/annotation/` | MITOS2 gene table of the 803M consensus (Table 2) and RSCU of the 13 protein-coding genes (`rscu_803m_v7.tsv`) |
| `results/validation/` | HiFi comparison of the v7 consensus (803M and D863F) |
| `results/world_tree_checks/` | BBG monophyly test and nearest-reference tables |
| `figures/` | Manuscript figures: Fig. 1A/B (ML trees), Fig. 2 (mismatch distribution), Fig. 3 (PopART median-joining network, SVG and PNG), Fig. 4 (read-haplotype phasing) |
| `scripts/` | All analysis scripts |

## Method in one paragraph

Reads overlapping chrM or NUMT loci were extracted (`remap_mt_reads.sh`), remapped to the mitochondrial
reference alone, and a majority-rule consensus was called with all reads counted equally
(`majority_consensus.py`). Mixed positions that recur across animals were shown to be a paralog signal by
read phasing (`phase_check.py`) and set to the allele of the PacBio HiFi mitogenome (`make_v4.py`, `make_v5.py`,
`make_v6.py`; `make_v7.py` adds positions 203 and 211). Trees: `rebuild_v6.sh` / `run_v7.sh` (full v7 pipeline). Diversity: `run_popgen_dloop_v6.sh`, `popgen_v6.R`.
Selection: `variant_codons.py`, `dnds_analysis.py`. Sanity checks: `verify_release.py`, `tree_checks.py`.

## Reproducing

`run_v7.sh` regenerates the v7 scripts from the v6 ones by replacing `v6` with `v7` (so both sets must be present). Scripts were written for the authors' server and contain absolute paths under `$HOME/mt_remap`.
Edit the path block at the top of each script before running. Software: minimap2, samtools, MAFFT, IQ-TREE 3.1.3,
MITOS2, R (ape, pegas), Python 3 (dendropy for `tree_checks.py`).
Raw reads are not included (see the manuscript's Data availability statement).

## Known limitations

- Positions 1-300 contain a nuclear paralog; the mitochondrial allele there is set from the HiFi mitogenome of 803M
  at 13 near-universal positions (11 from the first curation; 203 and 211 added after comparison with the HiFi, because phasing placed C on the HiFi-matching read haplotype). Positions that vary among animals are not independently confirmed.
- A 6-24 bp N-run in the control-region repeat (about 15,958-15,990) remains in 18 of 22 sequences.
- A poly-A run at position 11,559 (tRNA cluster between ND4 and ND5) has 9 A in the Illumina consensus and 8 in the 803M HiFi; the correct length is unresolved.
- Fu's Fs is close to its threshold and depends on how ambiguous positions are treated.
- Some script and file names keep an internal `v2` label (for example `dloop_v2_aligned.fa` inside the scripts);
  the files under `results/` are renamed. `results/curation/review_all22_*` is the v2 remap review.

## Provenance

HiFi mitogenomes: GenBank PZ809991 (803M) and PZ809992 (D863F). After curation, the v7 consensus of 803M matches its HiFi at all 16,633 comparable sites (`results/validation/`). WGS reads: BioProject PRJNA1502507 (SRA status: to be confirmed).
