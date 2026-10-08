import os,glob
H=os.path.expanduser('~/mt_remap/')
def rd(p):
    n=None;s=[]
    for l in open(p):
        if l.startswith('>'):
            if n: yield n,''.join(s)
            n=l[1:].strip();s=[]
        else: s.append(l.strip())
    yield n,''.join(s)
def patch(s,tag):
    s=list(s)
    for p in (203,211):
        if s[p-1].upper() not in 'TYC': print('UNEXPECTED',tag,p,s[p-1])
        s[p-1]='C'
    return ''.join(s)
for suf in ('v6','v6core'):
    new=suf.replace('v6','v7')
    recs=list(rd(H+'bbg_22S_mt_consensus_%s.fasta'%suf))
    os.makedirs(H+'fasta_%s'%new,exist_ok=True)
    with open(H+'bbg_22S_mt_consensus_%s.fasta'%new,'w') as o:
        for n,s in recs:
            t=patch(s,n)
            o.write('>%s\n%s\n'%(n,t))
            animal=n.split()[0]
            open(H+'fasta_%s/%s_mt.fasta'%(new,animal),'w').write('>%s\n%s\n'%(n,t))
    print(suf,'->',new,len(recs),'records')
