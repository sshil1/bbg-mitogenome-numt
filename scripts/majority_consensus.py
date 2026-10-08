#!/usr/bin/env python3
"""majority_consensus.py -- haploid majority-rule consensus (IUPAC code where minor allele >= --mixed-frac) from mitochondrial BAMs, with NO mapping-quality weighting.
Rationale: reads from the true mt allele that are identical to a NUMT get MAPQ~0 and are discounted by bcftools/vcf2fq,
so the NUMT allele can win. Counting every read equally lets the (much more abundant) mtDNA allele win.
Writes: FASTA (uppercase; N where depth < --min-depth) and a TSV of sites needing review (depth low or major allele < --review-frac).
Usage: python3 majority_consensus.py --ref REF.fna --out cons.fasta --report review.tsv [--contig NC_005044.2 --length 16643] BAM [BAM ...]
Needs samtools; read-only on BAMs. Sample name = BAM basename minus '_mt.bam'.
"""
import argparse, subprocess, re, os
from collections import Counter
IUPAC = {frozenset('AC'):'M', frozenset('AG'):'R', frozenset('AT'):'W', frozenset('CG'):'S', frozenset('CT'):'Y', frozenset('GT'):'K'}

def parse_bases(bases, refbase):
    out, i = [], 0
    while i < len(bases):
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            m = re.match(r'[+-](\d+)', bases[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        if c in '.,': out.append(refbase.upper())
        elif c in '*#': out.append('*')
        elif c in '<>': pass
        else: out.append(c.upper())
        i += 1
    return out

def consensus(bam, ref, contig, length, min_depth, review_frac, minq, mixed_frac):
    cmd = ['samtools', 'mpileup', '-f', ref, '-Q', str(minq), '-q', '0', '-B', '-d', '1000000', '-r', contig, bam]
    seq = ['N'] * length; review = []
    for line in subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines():
        f = line.split('\t')
        if len(f) < 5: continue
        pos = int(f[1]); cnt = Counter(b for b in parse_bases(f[4], f[2]) if b in 'ACGT')
        depth = sum(cnt.values())
        if depth == 0: continue
        base, n = cnt.most_common(1)[0]; frac = n / depth
        if depth >= min_depth:
            seq[pos-1] = base
            if len(cnt) > 1:
                (b2, n2) = cnt.most_common(2)[1]
                if n2 / depth >= mixed_frac: seq[pos-1] = IUPAC[frozenset((base, b2))]   # genuinely mixed site -> ambiguity code
        if depth < min_depth or frac < review_frac:
            review.append((pos, depth, cnt['A'], cnt['C'], cnt['G'], cnt['T'], seq[pos-1] if depth >= min_depth else 'N', round(frac, 2)))
    return ''.join(seq), review

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True); ap.add_argument('--out', required=True); ap.add_argument('--report', required=True)
    ap.add_argument('--contig', default='NC_005044.2'); ap.add_argument('--length', type=int, default=16643)
    ap.add_argument('--min-depth', type=int, default=5); ap.add_argument('--review-frac', type=float, default=0.8)
    ap.add_argument('--minq', type=int, default=20, help='minimum BASE quality (not mapping quality)')
    ap.add_argument('--mixed-frac', type=float, default=0.25, help='minor-allele fraction at/above which an IUPAC code is written instead of the majority base')
    ap.add_argument('bams', nargs='+'); a = ap.parse_args()
    with open(a.out, 'w') as fo, open(a.report, 'w') as fr:
        fr.write('sample\tpos\tdepth\tA\tC\tG\tT\tcalled\tmajor_frac\n')
        for bam in a.bams:
            s = os.path.basename(bam).replace('_mt.bam', '')
            seq, rev = consensus(bam, a.ref, a.contig, a.length, a.min_depth, a.review_frac, a.minq, a.mixed_frac)
            fo.write(f'>{s}\n{seq}\n')
            for r in rev: fr.write(s + '\t' + '\t'.join(map(str, r)) + '\n')
            print(f'{s}: N={seq.count("N")}  sites to review={len(rev)}')
