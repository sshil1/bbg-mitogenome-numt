import csv,collections,itertools
exec(open('gene_table_codon_checks.py').read().split("# codon check per animal")[0])  # genes, rd, v7, rc, PCG
# vertebrate mitochondrial code (table 2)
b='TCAG'; aa='FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'
code={a+b2+c:aa[i] for i,(a,b2,c) in enumerate(itertools.product(b,b,b))}
code['AGA']='*';code['AGG']='*'; code['ATA']='M'; code['TGA']='W'
fam=collections.defaultdict(list)
for c,a in code.items():
    if a!='*': fam[a].append(c)
def cds(s,g):
    seq=s[int(g['Start'])-1:int(g['Stop'])]
    if g['Strand']=='L': seq=rc(seq)
    if g['Stop_codon']=='T-': seq=seq+'AA'   # polyadenylation completes TAA
    return seq
def usage(s):
    cnt=collections.Counter(); third=collections.Counter()
    for g in genes:
        if g['Gene'] not in PCG: continue
        q=cds(s,g)
        for i in range(3,len(q)-3,3) if False else range(0,len(q)-3,3):
            c=q[i:i+3]
            if i==0: continue  # skip start codon
            if c in code and code[c]!='*': cnt[c]+=1; third[c[2]]+=1
    return cnt,third
cnt,third=usage(v7['803m'])
tot=sum(cnt.values()); print('codons counted (803m, excluding start/stop):',tot)
rows=[]
for a,cs in sorted(fam.items()):
    n=sum(cnt[c] for c in cs)
    for c in sorted(cs):
        rscu=(cnt[c]*len(cs)/n) if n else float('nan')
        rows.append((a,c,cnt[c],round(rscu,2)))
with open('rscu_803m_v7.tsv','w') as f:
    f.write('amino_acid\tcodon\tcount\tRSCU\n'); [f.write('\t'.join(map(str,r))+'\n') for r in rows]
t=sum(third.values()); print('third position %:',{k:round(v/t*100,1) for k,v in sorted(third.items())})
print('A+T at 3rd position %:',round((third['A']+third['T'])/t*100,1))
top=sorted([r for r in rows if len(fam[r[0]])>1],key=lambda r:-r[3])[:8]; print('highest RSCU:',top)
low=sorted([r for r in rows if len(fam[r[0]])>1 and r[2]>=0],key=lambda r:r[3])[:8]; print('lowest RSCU:',low)
print('unused codons (count 0):',[r[1] for r in rows if r[2]==0])
# same for all 22: total codon counts range and identity
tots=[sum(usage(s)[0].values()) for s in v7.values()]; print('codons counted per animal min/max',min(tots),max(tots))
# per-animal genome composition range (ACGT only)
rg=[]
for a,s in v7.items():
    c=collections.Counter(x for x in s if x in 'ACGT'); n=sum(c.values()); rg.append(((c['A']-c['T'])/(c['A']+c['T']),(c['G']-c['C'])/(c['G']+c['C']),(c['A']+c['T'])/n*100))
print('AT skew range',round(min(r[0] for r in rg),3),round(max(r[0] for r in rg),3),'GC skew',round(min(r[1] for r in rg),3),round(max(r[1] for r in rg),3),'AT%',round(min(r[2] for r in rg),2),round(max(r[2] for r in rg),2))
