import subprocess, tempfile, os, collections
def rd(p):
    k=None;d={}
    for l in open(os.path.expanduser(p)):
        l=l.strip()
        if l.startswith('>'): k=l[1:]; d[k]=[]
        elif k: d[k].append(l)
    return {a:''.join(b).upper() for a,b in d.items()}
rc=lambda s:s[::-1].translate(str.maketrans('ACGTRYKMSWN','TGCAYRMKSWN'))
ref=rd('~/mt_remap/ref_mt.fa')['NC_005044.2']
v2=rd('~/mt_remap/bbg_22S_mt_consensus_v2.fasta')
hf=list(rd('/data/gmbluser/mitohifi_C01/final_mitogenome.fasta').values())[0]
anc=ref[:20]
for s in (hf,rc(hf)):
    i=(s+s).find(anc)
    if 0<=i<len(s): hf=(s+s)[i:i+len(s)]; break
t=tempfile.mkdtemp(); f=t+'/a.fa'
open(f,'w').write(f'>ref\n{ref}\n>hifi\n{hf}\n')
a=subprocess.run(['/opt/conda/bin/mafft','--auto','--quiet',f],capture_output=True,text=True).stdout
q={};k=None
for l in a.splitlines():
    if l.startswith('>'): k=l[1:]; q[k]=[]
    else: q[k].append(l.upper())
q={x:''.join(y) for x,y in q.items()}
n=0;m={}
for i,c in enumerate(q['ref']):
    if c!='-': n+=1; m[n]=i
IU={'R':'AG','Y':'CT','K':'GT','M':'AC','S':'CG','W':'AT'}
SITES=[51,56,100,105,152,161,193,201,227,283,289]
out={};bad=[];fixed=0
for name,s in v2.items():
    s=list(s)
    for p in SITES:
        c=s[p-1]; h=q['hifi'][m[p]]
        if c in IU:
            if h in IU[c]: s[p-1]=h; fixed+=1
            else: bad.append((name,p,c,h))
    out[name]=''.join(s)
w=open(os.path.expanduser('~/mt_remap/bbg_22S_mt_consensus_v4.fasta'),'w')
for name,s in out.items(): w.write('>'+name+'\n'+s+'\n')
w.close()
print('resolved calls:',fixed,' incompatible with HiFi:',bad)
rem=collections.Counter(); oth=0
for s in out.values():
    for i,c in enumerate(s):
        if c in IU:
            if i<300: rem[i+1]+=1
            else: oth+=1
print('remaining IUPAC in 1-300 by site:',dict(sorted(rem.items())))
print('remaining IUPAC elsewhere (total):',oth)
