import os, collections
def rd(p):
    k=None;d={};h={}
    for l in open(os.path.expanduser(p)):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; h[k]=l[1:]; d[k]=[]
        elif k: d[k].append(l)
    return {a:''.join(b).upper() for a,b in d.items()},h
v4,hd=rd('~/mt_remap/bbg_22S_mt_consensus_v4.fasta')
v2,_=rd('~/mt_remap/bbg_22S_mt_consensus_v2.fasta')
IU=set('RYKMSWBDHV')
rows=[l.rstrip('\n').split('\t') for l in open(os.path.expanduser('~/mt_remap/review_all22.tsv'))][1:]
rev={(r[0],int(r[1])):r for r in rows}
out={};ch=0;keep=collections.Counter();miss=[]
for k,s in v4.items():
    s=list(s)
    for i,c in enumerate(s):
        if c in IU:
            r=rev.get((k,i+1))
            if not r: miss.append((k,i+1)); continue
            cnt=dict(zip('ACGT',map(int,r[3:7]))); mf=float(r[8]); mb=max(cnt,key=cnt.get)
            if mf>0.60: s[i]=mb; ch+=1
            else: keep[('1-300' if i<300 else 'Dloop' if i>=15400 else 'other')]+=1
    out[k]=''.join(s)
with open(os.path.expanduser('~/mt_remap/bbg_22S_mt_consensus_v5.fasta'),'w') as w:
    for k,s in out.items(): w.write('>'+hd[k]+'\n'+s+'\n')
print('IUPAC->majority:',ch,' kept IUPAC:',dict(keep),' not found in review:',miss[:5],len(miss))
# validation
L={len(s) for s in out.values()}; print('lengths:',L)
amb=lambda a,b:a!=b and a in 'ACGT' and b in 'ACGT'
d=sum(amb(a,b) for k in out for a,b in zip(v2[k],out[k]))
print('unambiguous diffs v2->v5:',d)
print('ATP8 TAA:',sum(s[7965:7968]=='TAA' for s in out.values()),'/22   COX2 stop TAA:',sum(s[7696:7699]=='TAA' for s in out.values()),'/22')
print('N total:',sum(s.count('N') for s in out.values()),' IUPAC total:',sum(c in IU for s in out.values() for c in s))
