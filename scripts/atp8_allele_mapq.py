#!/usr/bin/env python3
"""atp8_allele_mapq.py -- allele counts and mean MAPQ per allele at one mitochondrial position.
Read-only. Needs only samtools (>=1.9). No reference FASTA needed: '.'/',' = matches the reference base.
Usage: python3 atp8_allele_mapq.py [--contig NC_005044.2] [--pos 7968] [--ref A] BAM [BAM ...]
"""
import argparse, subprocess, re, os, sys
from collections import defaultdict

def parse_pileup(bases, mapqs):
    """Return list of (allele, mapq) for each read in one mpileup line (samtools mpileup -s)."""
    out, i, r = [], 0, 0
    while i < len(bases):
        c = bases[i]
        if c == '^': i += 2; continue            # read start + mapq char (not a new read)
        if c == '$': i += 1; continue            # read end
        if c in '+-':                            # indel attached to previous read
            m = re.match(r'[+-](\d+)', bases[i:]); n = int(m.group(1)); i += len(m.group(0)) + n; continue
        a = 'ref' if c in '.,' else c.upper()
        out.append((a, ord(mapqs[r]) - 33)); r += 1; i += 1
    return out

def run(bam, contig, pos):
    cmd = ['samtools', 'mpileup', '-Q', '0', '-q', '0', '-B', '-d', '1000000', '-s',
           '-r', f'{contig}:{pos}-{pos}', bam]
    p = subprocess.run(cmd, capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split('\t')
        if len(f) >= 7 and int(f[1]) == pos:
            return parse_pileup(f[4], f[6])
    return []

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--contig', default='NC_005044.2'); ap.add_argument('--pos', type=int, default=7968)
    ap.add_argument('--ref', default='A', help='reference base at the position (used only if samtools prints ./,)')
    ap.add_argument('bams', nargs='+'); a = ap.parse_args()
    print('\t'.join(['sample', 'depth', 'A(n/meanMQ)', 'C(n/meanMQ)', 'G(n/meanMQ)', 'T(n/meanMQ)',
                     'G_frac', 'reads_MAPQ0_frac', 'reads_MAPQ>=20_frac']))
    for bam in a.bams:
        reads = run(bam, a.contig, a.pos); s = os.path.basename(bam).replace('_mt.bam', '')
        by = defaultdict(list)
        for al, q in reads: by[a.ref.upper() if al == 'ref' else al].append(q)
        cell = lambda b: f"{len(by[b])}/{sum(by[b])/len(by[b]):.1f}" if by[b] else '0/NA'
        n = len(reads); allq = [q for _, q in reads]
        f = lambda x: f'{x:.2f}' if n else 'NA'
        print('\t'.join(map(str, [s, n, cell('A'), cell('C'), cell('G'), cell('T'),
              f(len(by['G'])/n if n else 0), f(sum(q == 0 for q in allq)/n if n else 0), f(sum(q >= 20 for q in allq)/n if n else 0)])))
