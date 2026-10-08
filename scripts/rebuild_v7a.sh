#!/usr/bin/env bash
set -euo pipefail
MAFFT=${MAFFT:?set MAFFT}; IQ=/opt/conda/envs/variants/bin/iqtree3; T=8
OUT=$HOME/mt_remap/trees_v7; mkdir -p $OUT; cd $OUT
python3 - <<'PY'
import os,sys
H=os.path.expanduser('~'); AL='/data/mt/alignments/'
def rd(p):
    d={};k=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): k=l[1:]; d[k]=[]
        elif k is not None: d[k].append(l)
    return {a:''.join(b) for a,b in d.items()}
def first(d): return {k.split()[0]:s.upper() for k,s in d.items()}
def wr(p,items):
    with open(p,'w') as f:
        for k,s in items: f.write('>'+k+'\n'+'\n'.join(s[i:i+80] for i in range(0,len(s),80))+'\n')
COMP=str.maketrans('ACGT','TGCA')
def rc(s): return s.translate(COMP)[::-1]
v2=first(rd(H+'/mt_remap/bbg_22S_mt_consensus_v7.fasta')); old=first(rd('/data/mt/bbg_22S_mt_consensus.fasta'))
names=list(v2); assert len(names)==22
REF=v2['803m']   # NC_005044.2 coordinates
def orient(s):
    for strand,t in (('+',s),('-',rc(s))):
        for o in range(0,len(REF)-20,500):
            a=REF[o:o+20]; p=t.find(a)
            if p>=0 and t.find(a,p+1)<0:
                st=(p-o)%len(t); return strand,st,t[st:]+t[:st]
    return None,0,s
def others(d,tag):
    out=[]
    for k,s in d.items():
        if k.split()[0] in v2: continue
        s=s.replace('-','').upper(); strand,st,s2=orient(s)
        if strand is None: print(f'[{tag}] NO ANCHOR (left as is): {k.split()[0]} len={len(s)}')
        elif strand=='-' or (st>20 and len(s)-st>20): print(f'[{tag}] rotated/flipped: {k.split()[0]} strand={strand} shift={st}')
        else: s2=s if strand=='+' else s2
        out.append((k,s2))
    return out
def check(d,tag):
    got={k.split()[0] for k in d}; miss=[n for n in names if n not in got]
    if miss: sys.exit(f'{tag}: BBG names missing {miss}')
c=rd(AL+'bbg_complete_raw.fasta'); check(c,'complete'); wr('complete_v2_raw.fa',[(n,v2[n]) for n in names]+others(c,'complete'))
d=rd(AL+'bbg_dloop_raw.fasta'); check(d,'dloop'); dd=first(d); items=[]
for n in names:
    seg=dd[n].replace('-',''); o=old[n]; a=o.find(seg[:30]); b=o.find(seg[-30:])
    if a<0 or b<0: sys.exit(f'dloop anchor not found for {n}')
    items.append((n,v2[n][a:b+30]))
wr('dloop_v2_raw.fa',items+[(k,s.replace('-','').upper()) for k,s in d.items() if k.split()[0] not in v2])
w=rd(AL+'bbg_world_aligned.fasta'); check(w,'world'); wr('world_v2_raw.fa',[(n,v2[n]) for n in names]+others(w,'world'))
for f in ('complete','dloop','world'): print(f, sum(1 for l in open(f+'_v2_raw.fa') if l[0]=='>'),'sequences')
PY
