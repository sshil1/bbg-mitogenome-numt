import subprocess, tempfile, os, itertools, collections
def rd(p):
    k=None;d={};h={}
    for l in open(os.path.expanduser(p)):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; h[k]=l[1:]; d[k]=[]
        elif k: d[k].append(l)
    return {a:''.join(b).upper() for a,b in d.items()},h
rc=lambda s:s[::-1].translate(str.maketrans('ACGTRYKMSWN','TGCAYRMKSWN'))
ref=rd('~/mt_remap/ref_mt.fa')[0]['NC_005044.2']
v5,hd=rd('~/mt_remap/bbg_22S_mt_consensus_v5.fasta')
hf=list(rd('/data/gmbluser/mitohifi_C01/final_mitogenome.fasta')[0].values())[0]
for s in (hf,rc(hf)):
    i=(s+s).find(ref[:20])
    if 0<=i<len(s): hf=(s+s)[i:i+len(s)]; break
t=tempfile.mkdtemp(); f=t+'/a.fa'; open(f,'w').write(f'>ref\n{ref}\n>hifi\n{hf}\n')
a=subprocess.run(['/opt/conda/bin/mafft','--auto','--quiet',f],capture_output=True,text=True).stdout
q={};k=None
for l in a.splitlines():
    if l.startswith('>'): k=l[1:]; q[k]=[]
    else: q[k].append(l.upper())
q={x:''.join(y) for x,y in q.items()}
n=0;m={}
for i,c in enumerate(q['ref']):
    if c!='-': n+=1; m[n]=i
SITES=[51,56,100,105,152,161,193,201,227,283,289]
chg=[];out={}
for k,s in v5.items():
    s=list(s)
    for p in SITES:
        h=q['hifi'][m[p]]
        if s[p-1]!=h: chg.append((k,p,s[p-1],h)); s[p-1]=h
    out[k]=''.join(s)
with open(os.path.expanduser('~/mt_remap/bbg_22S_mt_consensus_v6.fasta'),'w') as w:
    for k,s in out.items(): w.write('>'+hd[k]+'\n'+s+'\n')
print('forced to HiFi allele (animal,pos,old,new):',chg)
def stats(d,cols):
    S=list(d.values()); seg=0; tot=0; np_=0
    for i in cols:
        if len({s[i] for s in S if s[i] in 'ACGT'})>1: seg+=1
    for a,b in itertools.combinations(S,2):
        v=[(a[i],b[i]) for i in cols if a[i] in 'ACGT' and b[i] in 'ACGT']
        tot+=sum(x!=y for x,y in v)/len(v); np_+=1
    cc=[i for i in cols if all(s[i] in 'ACGT' for s in S)]
    hap=len({''.join(s[i] for i in cc) for s in S})
    return seg, tot/np_, hap, len(cols)
allc=range(16643)
print('v6 all columns        S/pi/hap/len:',stats(out,allc))
print('v6 excl 1-300         S/pi/hap/len:',stats(out,range(300,16643)))
print('v6 excl 1-300,>=15400 S/pi/hap/len:',stats(out,range(300,15400)))
rem=collections.Counter()
for s in out.values():
    for i in range(300):
        pass
vs=[i+1 for i in range(300) if len({s[i] for s in out.values() if s[i] in 'ACGT'})>1]
print('variable sites left in 1-300:',vs)
print('ATP8 TAA:',sum(s[7965:7968]=='TAA' for s in out.values()),'/22')
