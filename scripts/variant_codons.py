import csv,sys,collections
FA,GT,MAXMISS=sys.argv[1],sys.argv[2],int(sys.argv[3])
seqs={};k=None
for l in open(FA):
    l=l.strip()
    if l.startswith('>'): k=l[1:].split()[0]; seqs[k]=[]
    elif k: seqs[k].append(l.upper())
seqs={a:''.join(b) for a,b in seqs.items()}; names=list(seqs); ACGT=set('ACGT')
order='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
T={a+b+c:aa[16*i+4*j+k] for i,a in enumerate(order) for j,b in enumerate(order) for k,c in enumerate(order)}
T['AGA']='*';T['AGG']='*';T['ATA']='M';T['TGA']='W'
comp=str.maketrans('ACGT','TGCA')
PCG="ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB".split()
res={g:[0,0,0] for g in PCG}   # syn, nonsyn, excluded
for r in csv.DictReader(open(GT)):
    g=r['Gene']
    if g not in res: continue
    s,e,st=int(r['Start']),int(r['Stop']),r['Strand']
    idx=list(range(s-1,e)) if st=='H' else list(range(e-1,s-2,-1))
    for c in range(0,len(idx)-2,3):
        pos=idx[c:c+3]; cods=[]
        for n in names:
            x=''.join(seqs[n][p] for p in pos)
            cods.append(x if st=='H' else x.translate(comp))
        valid=[x for x in cods if all(y in ACGT for y in x)]
        if len(cods)-len(valid)>MAXMISS: res[g][2]+=1; continue
        if len(set(valid))>1:
            res[g][0 if len({T[x] for x in valid})==1 else 1]+=1
print(f'max missing animals per codon = {MAXMISS}')
print('gene    syn nonsyn total  codons_excluded')
S=N=0
for g in PCG:
    a,b,x=res[g]; S+=a; N+=b; print(f'{g:6s} {a:4d} {b:5d} {a+b:5d}  {x}')
print(f'TOTAL  {S:4d} {N:5d} {S+N:5d}')
