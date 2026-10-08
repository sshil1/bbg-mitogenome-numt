#!/usr/bin/env python3
"""mismatch_sites.py -- where does a HiFi mitogenome differ from the Illumina consensus of the same-named animal, and why?
For every mismatch it reports: reference coordinate, Illumina base (case kept: lowercase = low-confidence in vcf2fq),
HiFi base, and the reference (NC_005044.2) base, then summarises how many Illumina calls equal the reference / are lowercase.
Usage: python3 mismatch_sites.py CONSENSUS.fasta SAMPLE HIFI.fasta REF_MT.fasta   (needs mafft; read-only)
"""
import sys, subprocess, tempfile, os
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
def revcomp(x): return x.translate(COMP)[::-1]
def read_fasta(path):
    d, n = {}, None
    for l in open(path):
        l = l.rstrip()
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
def orient_to_ref(seq, ref, k=25):
    ref = ref.upper(); L = len(seq)
    for strand, s in (('+', seq), ('-', revcomp(seq))):
        dbl = s.upper() * 2
        for p in range(100, len(ref) - k, 1500):
            kmer = ref[p:p+k]
            if set(kmer) <= set('ACGT'):
                q = dbl.find(kmer)
                if q != -1: st = (q - p) % L; return s[st:] + s[:st]
    sys.exit('could not anchor HiFi sequence to the reference')
def analyze(ill, hifi, ref):
    """aligned, equal-length strings; returns (rows, comparable_sites, lowercase_total)"""
    rows, pos, comp, lower_all = [], 0, 0, 0
    for a, h, r in zip(ill, hifi, ref):
        if a != '-': pos += 1
        if a == '-' or h == '-': continue
        if a.upper() in 'ACGT' and h.upper() in 'ACGT':
            comp += 1; lower_all += a.islower()
            if a.upper() != h.upper(): rows.append((pos, a, h.upper(), r.upper() if r != '-' else '-'))
    return rows, comp, lower_all
def main():
    cons_f, sample, hifi_f, ref_f = sys.argv[1:5]
    cons = read_fasta(cons_f)[sample]; hifi = next(iter(read_fasta(hifi_f).values())); ref = next(iter(read_fasta(ref_f).values()))
    hifi = orient_to_ref(hifi, cons if set(cons[100:2000].upper()) <= set('ACGT') else ref)
    with tempfile.TemporaryDirectory() as t:
        p = os.path.join(t, 'in.fa'); open(p, 'w').write(f'>ill\n{cons}\n>hifi\n{hifi}\n>ref\n{ref}\n')
        out = subprocess.run(['mafft', '--auto', '--quiet', '--preservecase', p], capture_output=True, text=True, check=True).stdout
    q = os.path.join(tempfile.gettempdir(), '_ms.fa'); open(q, 'w').write(out); al = read_fasta(q)
    rows, comp, low_all = analyze(al['ill'], al['hifi'], al['ref'])
    n = len(rows); eq_ref = sum(a.upper() == r for _, a, h, r in rows); low = sum(a.islower() for _, a, h, r in rows)
    hifi_ref = sum(h == r for _, a, h, r in rows); cr = sum(pos > 15430 for pos, *_ in rows)
    print(f'{sample}: {n} mismatches over {comp} comparable sites')
    print(f'  Illumina base == reference base : {eq_ref}/{n}')
    print(f'  HiFi base == reference base     : {hifi_ref}/{n}')
    print(f'  Illumina base lowercase         : {low}/{n}   (background: {low_all}/{comp} = {100*low_all/comp:.0f}% of all comparable sites)')
    print(f'  in control region (>15,430)     : {cr}/{n}')
    print('pos\tIllumina\tHiFi\tref')
    for r in rows: print('\t'.join(map(str, r)))
if __name__ == '__main__': main()
