#!/usr/bin/env python3
"""compare_pacbio_to_consensus.py -- (1) which of the 22 Illumina consensus mitogenomes is each PacBio mitogenome closest to?
(2) what base does each sequence carry at reference position --pos (default NC_005044.2:7968, the ATP8 stop codon)?
Read-only. Needs mafft on PATH; no Biopython needed. Accepts FASTA or GenBank (.gb) for the PacBio files.
Usage: python3 compare_pacbio_to_consensus.py CONSENSUS.fasta LABEL=PATH [LABEL=PATH ...]
  e.g. python3 compare_pacbio_to_consensus.py /data/mt/bbg_22S_mt_consensus.fasta \
         C01_dir=/data/gmbluser/mitohifi_C01/final_mitogenome.fasta PZ809991=/data/gmbluser/mitohifi_C01/FINAL/PZ809991.gb
"""
import sys, subprocess, tempfile, os, re

def read_fasta(path):
    d, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith('>'): name = line[1:].split()[0]; d[name] = []
        elif name: d[name].append(line.strip())
    return {k: ''.join(v) for k, v in d.items()}


COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
def revcomp(x): return x.translate(COMP)[::-1]

def orient_to_ref(seq, ref, k=25):
    """Rotate (and reverse-complement if needed) a circular mitogenome so it starts like `ref`. Returns (seq, note)."""
    ref = ref.upper(); L = len(seq)
    for strand, s in (('+', seq), ('-', revcomp(seq))):
        dbl = s.upper() * 2
        for p in range(100, len(ref) - k, 1500):
            kmer = ref[p:p+k]
            if set(kmer) <= set('ACGT'):
                q = dbl.find(kmer)
                if q != -1:
                    start = (q - p) % L
                    return s[start:] + s[:start], f'strand {strand}, rotated by {start}'
    return seq, 'WARNING: no anchor found; left as is'

def read_gb(path):
    seq, on = [], False
    for line in open(path):
        if line.startswith('ORIGIN'): on = True; continue
        if line.startswith('//'): break
        if on: seq.append(re.sub(r'[^A-Za-z]', '', line))
    return ''.join(seq)

def load_pacbio(label, path):
    if path.lower().endswith(('.gb', '.gbk', '.genbank')): return {label: read_gb(path)}
    f = read_fasta(path)
    return {label + ('' if len(f) == 1 else ':' + k): v for k, v in f.items()}

def main():
    pos = 7968
    args = sys.argv[1:]
    if '--pos' in args: i = args.index('--pos'); pos = int(args[i+1]); del args[i:i+2]
    cons = read_fasta(args[0]); pb = {}
    for a in args[1:]:
        label, path = a.split('=', 1); pb.update(load_pacbio(label, path))
    refseq = next(v for v in cons.values() if set(v[100:2000].upper()) <= set('ACGT')) if any(set(v[100:2000].upper()) <= set('ACGT') for v in cons.values()) else next(iter(cons.values()))
    for k in list(pb):
        pb[k], note = orient_to_ref(pb[k], refseq); print(f'{k}: {note}')
    with tempfile.TemporaryDirectory() as t:
        inp = os.path.join(t, 'in.fa')
        with open(inp, 'w') as f:
            for k, v in {**cons, **pb}.items(): f.write(f'>{k}\n{v}\n')
        aln = subprocess.run(['mafft', '--auto', '--quiet', inp], capture_output=True, text=True, check=True).stdout
    tmp = os.path.join(tempfile.gettempdir(), '_aln.fa'); open(tmp, 'w').write(aln); al = read_fasta(tmp)
    ref = next(iter(cons)); col, n = None, 0               # alignment column of reference position `pos` (consensus is in reference coordinates)
    for c, ch in enumerate(al[ref]):
        if ch != '-': n += 1
        if n == pos: col = c; break
    ok = set('ACGT')
    def diffs(a, b):
        m = t = 0
        for x, y in zip(a.upper(), b.upper()):
            if x in ok and y in ok: t += 1; m += (x != y)
        return m, t
    print(f'alignment column for position {pos}: {col+1}')
    print('\nBase at position', pos, '(context +/-5 in alignment):')
    for k in list(cons) + list(pb): print(f'  {k:28s} {al[k][col-5:col]}[{al[k][col]}]{al[k][col+1:col+6]}')
    for p in pb:
        res = sorted((diffs(al[p], al[c]) + (c,) for c in cons), key=lambda r: r[0])
        print(f'\nClosest Illumina consensus animals to {p} (mismatches over comparable ACGT sites):')
        for m, t, c in res[:5]: print(f'  {c:12s} {m:4d} mismatches / {t} sites')
main()
