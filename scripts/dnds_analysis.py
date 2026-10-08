#!/usr/bin/env python3
"""
dN/dS (Nei-Gojobori 1986 method) per protein-coding gene, across all 22 BBG
individuals -- tests whether each gene shows the purifying selection (dN/dS < 1)
expected for mitochondrial PCGs, or an unusual signal (dN/dS >= 1, suggesting
relaxed constraint or positive selection).

Usage:
    python3 dnds_analysis.py gene_table.csv fasta_dir/ out.csv
"""
import glob
import math
import os
import sys
import itertools
from collections import defaultdict

bases = ["T", "C", "A", "G"]
_aa_std = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON_TABLE[b1 + b2 + b3] = _aa_std[i]
            i += 1
CODON_TABLE["AGA"] = "*"
CODON_TABLE["AGG"] = "*"
CODON_TABLE["ATA"] = "M"
CODON_TABLE["TGA"] = "W"


def translate(codon):
    return CODON_TABLE.get(codon)


def synonymous_sites(codon):
    """Fraction of synonymous sites in this codon (Nei-Gojobori style):
    for each position, of the 3 possible single-base changes, what fraction
    leave the amino acid unchanged."""
    aa0 = translate(codon)
    if aa0 is None or aa0 == "*":
        return None
    s = 0.0
    for pos in range(3):
        syn_count = 0
        for b in bases:
            if b == codon[pos]:
                continue
            new_codon = codon[:pos] + b + codon[pos+1:]
            aa_new = translate(new_codon)
            if aa_new == aa0:
                syn_count += 1
        s += syn_count / 3.0
    return s  # out of 3 possible "site-equivalents" for this codon


def codon_diff_counts(codon1, codon2):
    """Nei-Gojobori pairwise comparison between two codons that differ.
    Enumerates all mutational paths between them (all orderings of the
    differing positions), classifies each single-step substitution as
    synonymous or nonsynonymous, and averages Sd/Nd across all equally
    parsimonious paths."""
    diffs = [i for i in range(3) if codon1[i] != codon2[i]]
    if not diffs:
        return 0.0, 0.0
    if len(diffs) == 1:
        aa1, aa2 = translate(codon1), translate(codon2)
        if aa1 is None or aa2 is None:
            return None, None
        return (0.0, 1.0) if aa1 != aa2 else (1.0, 0.0)

    total_sd, total_nd, n_paths = 0.0, 0.0, 0
    for order in itertools.permutations(diffs):
        path_codon = codon1
        sd, nd = 0.0, 0.0
        valid = True
        for pos in order:
            next_codon = path_codon[:pos] + codon2[pos] + path_codon[pos+1:]
            aa_before = translate(path_codon)
            aa_after = translate(next_codon)
            if aa_before is None or aa_after is None:
                valid = False
                break
            if aa_after == "*" and next_codon != codon2:
                # a premature stop mid-path -- still count it, NG86 doesn't
                # discard these paths by default in the simple implementation
                pass
            if aa_after == aa_before:
                sd += 1
            else:
                nd += 1
            path_codon = next_codon
        if valid:
            total_sd += sd
            total_nd += nd
            n_paths += 1
    if n_paths == 0:
        return None, None
    return total_sd / n_paths, total_nd / n_paths


def jukes_cantor_correct(p):
    if p is None or p >= 0.75:
        return None
    try:
        return -0.75 * math.log(1 - (4.0/3.0) * p)
    except ValueError:
        return None


def read_single_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    gene_table_path, fasta_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    genes = []
    with open(gene_table_path) as fh:
        import csv
        for row in csv.DictReader(fh):
            if row["Gene"] in ("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
                                "COX1","COX2","COX3","ATP6","ATP8","CYTB"):
                genes.append({"name": row["Gene"], "start": int(row["Start"]),
                              "end": int(row["Stop"]), "strand": row["Strand"]})

    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")) +
                          glob.glob(os.path.join(fasta_dir, "*.fa")))
    if not fasta_files:
        print(f"No FASTA files in {fasta_dir}")
        sys.exit(1)

    animals = {}
    for f in fasta_files:
        name = os.path.splitext(os.path.basename(f))[0]
        animals[name] = read_single_fasta(f)

    lens = set(len(s) for s in animals.values())
    if len(lens) > 1:
        print(f"WARNING: sequence lengths differ: {lens} -- results below may be unreliable "
              f"if this reflects real indels rather than trailing whitespace/truncation.")

    results = []
    for g in genes:
        gene_seqs = {}
        for name, s in animals.items():
            if len(s) < g["end"]:
                continue
            sub = s[g["start"]-1:g["end"]]
            if g["strand"] == "L":
                sub = revcomp(sub)
            gene_seqs[name] = sub

        if len(gene_seqs) < 2:
            continue

        gene_len = len(next(iter(gene_seqs.values())))
        n_codons = gene_len // 3

        # average synonymous/nonsynonymous SITE counts across all animals
        S_total, N_total, n_valid_seqs = 0.0, 0.0, 0
        for name, seq in gene_seqs.items():
            s_sum = 0.0
            valid = True
            for c in range(n_codons):
                codon = seq[c*3:c*3+3]
                if any(ch not in "ACGT" for ch in codon):
                    continue
                s = synonymous_sites(codon)
                if s is None:
                    continue
                s_sum += s
            S_total += s_sum
            N_total += (n_codons * 3) - s_sum  # remainder are nonsyn sites (rough per-seq total, ignoring skipped codons' small effect)
            n_valid_seqs += 1
        if n_valid_seqs == 0:
            continue
        S_avg = S_total / n_valid_seqs
        N_avg = N_total / n_valid_seqs

        # average pairwise Sd, Nd across all animal pairs
        names = list(gene_seqs.keys())
        Sd_sum, Nd_sum, n_pairs = 0.0, 0.0, 0
        for a, b in itertools.combinations(names, 2):
            seq_a, seq_b = gene_seqs[a], gene_seqs[b]
            pair_sd, pair_nd = 0.0, 0.0
            skip_pair = False
            for c in range(n_codons):
                ca = seq_a[c*3:c*3+3]
                cb = seq_b[c*3:c*3+3]
                if any(ch not in "ACGT" for ch in ca) or any(ch not in "ACGT" for ch in cb):
                    continue
                if ca == cb:
                    continue
                sd, nd = codon_diff_counts(ca, cb)
                if sd is None:
                    continue
                pair_sd += sd
                pair_nd += nd
            Sd_sum += pair_sd
            Nd_sum += pair_nd
            n_pairs += 1

        if n_pairs == 0:
            continue
        Sd_avg = Sd_sum / n_pairs
        Nd_avg = Nd_sum / n_pairs

        pS = Sd_avg / S_avg if S_avg > 0 else None
        pN = Nd_avg / N_avg if N_avg > 0 else None
        dS = jukes_cantor_correct(pS) if pS is not None else None
        dN = jukes_cantor_correct(pN) if pN is not None else None

        ratio = None
        note = ""
        if dS is not None and dN is not None:
            if dS == 0:
                note = "dS=0 -- ratio undefined (no synonymous divergence observed)"
            else:
                ratio = dN / dS

        results.append({
            "Gene": g["name"], "N_animals": n_valid_seqs, "N_codons": n_codons,
            "S_sites": round(S_avg, 2), "N_sites": round(N_avg, 2),
            "pS": round(pS, 5) if pS is not None else "NA",
            "pN": round(pN, 5) if pN is not None else "NA",
            "dS": round(dS, 5) if dS is not None else "NA",
            "dN": round(dN, 5) if dN is not None else "NA",
            "dN_dS_ratio": round(ratio, 4) if ratio is not None else "NA",
            "Note": note,
        })

    with open(out_path, "w") as out:
        import csv
        writer = csv.DictWriter(out, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {out_path}\n")
    print(f"{'Gene':<8}{'dN':<10}{'dS':<10}{'dN/dS':<10}{'Note'}")
    for r in results:
        print(f"{r['Gene']:<8}{str(r['dN']):<10}{str(r['dS']):<10}{str(r['dN_dS_ratio']):<10}{r['Note']}")
    print("\ndN/dS << 1 = purifying (stabilizing) selection, the norm for mtDNA PCGs.")
    print("dN/dS close to or above 1 would suggest relaxed constraint or positive selection.")


if __name__ == "__main__":
    main()
