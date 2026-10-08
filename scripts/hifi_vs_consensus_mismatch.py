import os,subprocess
def rd(p):
    n=None;s=[]
    for l in open(os.path.expanduser(p)):
        if l.startswith('>'):
            if n: yield n,''.join(s)
            n=l[1:].split()[0];s=[]
        else: s.append(l.strip())
    yield n,''.join(s)
cons=dict(rd('~/mt_remap/bbg_22S_mt_consensus_v6.fasta'))
h=next(rd('~/mt_remap/803M_HiFi.fasta'))[1].upper()
def rc(s): return s[::-1].translate(str.maketrans('ACGT','TGCA'))
L=len(cons['803m'])
for animal in ['803m','814m']:
    c=cons[animal].upper()
    anchor=None
    for st in range(300,2000,7):
        if set(c[st:st+30])<=set('ACGT'): anchor=st;break
    a=c[anchor:anchor+30]
    hh=h;i=hh.find(a)
    if i<0: hh=rc(h);i=hh.find(a)
    if i<0: print(animal,'anchor not found');continue
    cr=c[anchor:]+c[:anchor]; hr=hh[i:]+hh[:i]
    fa=os.path.expanduser('~/mt_remap/_hm_%s.fa'%animal)
    open(fa,'w').write('>c\n%s\n>h\n%s\n'%(cr,hr))
    ali=subprocess.run(['/opt/conda/bin/mafft','--quiet','--auto',fa],capture_output=True,text=True).stdout
    seqs={};n=None
    for l in ali.splitlines():
        if l.startswith('>'): n=l[1];seqs[n]=[]
        else: seqs[n].append(l.strip())
    ca=''.join(seqs['c']).upper();ha=''.join(seqs['h']).upper()
    pos=0;print('==',animal,'anchor',anchor+1,'HiFi len',len(h),'cons len',L)
    for x,y in zip(ca,ha):
        if x!='-': pos+=1
        orig=(anchor+pos-1)%L+1
        if x==y: continue
        kind='MISMATCH' if (x in 'ACGT' and y in 'ACGT') else ('indel' if '-' in (x,y) else 'ambig')
        print(kind,'v6pos',orig,'cons',x,'HiFi',y)
