# Run from the package root: python3 -I scripts/verify_release.py
import os,glob,collections
def rd(p):
    k=None;d={}
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {a:''.join(b).upper() for a,b in d.items()}
IU=set('RYKMSWBDHV')
v2=rd('data/provenance/bbg_22S_mt_consensus_v2.fasta');v5=rd('data/provenance/bbg_22S_mt_consensus_v5.fasta')
v6=rd('data/provenance/bbg_22S_mt_consensus_v6.fasta');v7=rd('data/consensus_v7/bbg_22S_mt_consensus_v7.fasta');c7=rd('data/consensus_v7/bbg_22S_mt_consensus_v7core.fasta')
print('names equal v2/v5/v6/v7/v7core:',list(v2)==list(v5)==list(v6)==list(v7)==list(c7),len(v7))
print('lengths',{len(s) for s in v7.values()},' charset',sorted(set(''.join(v7.values()))))
print('IUPAC v2/v5/v6/v7:',[sum(c in IU for s in d.values() for c in s) for d in (v2,v5,v6,v7)],' N:',[sum(s.count('N') for s in d.values()) for d in (v2,v5,v6,v7)])
chg=collections.Counter((i+1,a,b) for k in v7 for i,(a,b) in enumerate(zip(v6[k],v7[k])) if a!=b)
print('v6->v7 changes:',sorted(chg.items()))
print('ATP8 TAA',sum(s[7965:7968]=='TAA' for s in v7.values()),'COX2 TAA',sum(s[7696:7699]=='TAA' for s in v7.values()))
print('positions 203/211:',collections.Counter(s[202] for s in v7.values()),collections.Counter(s[210] for s in v7.values()))
bad=0
for k in v7:
    r=rd('data/per_animal_v7/%s_mt.fasta'%k)
    if r[k]!=v7[k]: bad+=1
    r=rd('data/per_animal_v7core/%s_mt.fasta'%k)
    d=[i+1 for i,(a,b) in enumerate(zip(v7[k],r[k])) if a!=b]
    if set(d)-{4,14526,16495}: bad+=1
print('per-animal fasta mismatches:',bad,len(glob.glob('data/per_animal_v7/*')),len(glob.glob('data/per_animal_v7core/*')))
for f in sorted(glob.glob('results/alignments/*.fa')):
    r=rd(f); L={len(s) for s in r.values()}
    print(f'{os.path.basename(f):28s} taxa={len(r):3d} len={L} moth={"KP677508.1" in r} KR059186={"KR059186.1" in r} KR059195={"KR059195.1" in r}')
