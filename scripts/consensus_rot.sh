#!/usr/bin/env bash
# consensus_rot.sh : majority consensus on the rotated reference, rotated back to NC_005044.2 coordinates, compared with v2
set -euo pipefail
SHIFT=${SHIFT:-8000}; R=$HOME/mt_remap_rot
python3 $HOME/majority_consensus.py --ref "$R/ref_mt_rot.fa" --out "$R/consensus_rot.fasta" --report "$R/review_rot.tsv" "$R"/*/*_mt.bam
python3 - "$R" "$SHIFT" "$HOME/mt_remap/bbg_22S_mt_consensus_v2.fasta" <<'PY'
import sys,csv,re,collections
R,SH,V2=sys.argv[1],int(sys.argv[2]),sys.argv[3]
def rd(p):
    d={};k=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=[]
        elif k is not None: d[k].append(l)
    return {a:''.join(b).upper() for a,b in d.items()}
rot=rd(R+'/consensus_rot.fasta'); v2=rd(V2); L=16643
v3={n:s[L-SH:]+s[:L-SH] for n,s in rot.items()}
assert all(len(s)==L for s in v3.values()), "unexpected length"
with open(R+'/bbg_22S_mt_consensus_v3.fasta','w') as f:
    for n in v2:
        s=v3[n]; f.write('>'+n+'\n'+'\n'.join(s[i:i+80] for i in range(0,L,80))+'\n')
rows=list(csv.reader(open(R+'/review_rot.tsv'),delimiter='\t')); h=rows[0]; pi=h.index('pos')
out=[h]
for r in rows[1:]:
    r=list(r); r[pi]=str((int(r[pi])-1+SH)%L+1); out.append(r)
out=[out[0]]+sorted(out[1:],key=lambda r:(r[0],int(r[pi])))
csv.writer(open(R+'/review_v3.tsv','w'),delimiter='\t',lineterminator='\n').writerows(out)
A=set('ACGT')
print("animal   diff_vs_v2(unambig)  v3_IUPAC(1-300/other)  v3_N   ATP8  v2_IUPAC(1-300)")
tot=collections.Counter()
for n in v2:
    a,b=v2[n],v3[n]
    d=sum(1 for x,y in zip(a,b) if x in A and y in A and x!=y)
    iu=[i for i,c in enumerate(b) if c not in A and c!='N']; nn=b.count('N')
    iu2=[i for i,c in enumerate(a) if c not in A and c!='N']
    print(f"{n:9s} {d:6d} {sum(1 for i in iu if i<300):12d}/{sum(1 for i in iu if i>=300):<4d} {nn:10d}   {b[7965:7968]}  {sum(1 for i in iu2 if i<300):8d}")
    tot['v3_iupac_start']+=sum(1 for i in iu if i<300); tot['v3_iupac_other']+=sum(1 for i in iu if i>=300)
    tot['v2_iupac_start']+=sum(1 for i in iu2 if i<300); tot['v2_iupac_other']+=sum(1 for i in iu2 if i>=300)
print(dict(tot)); print("wrote",R+"/bbg_22S_mt_consensus_v3.fasta")
PY
