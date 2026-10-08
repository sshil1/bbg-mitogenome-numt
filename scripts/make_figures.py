import dendropy, itertools, collections, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
BL='#0072B2'; OR='#D55E00'; GR='#8a8a8a'; SK='#56B4E9'; YL='#E69F00'
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
V='v6x/mt_remap/'
bbg=[l[1:].split()[0] for l in open(V+'bbg_22S_mt_consensus_v6.fasta') if l.startswith('>')]
B=set(bbg)

def draw_tree(path,out,title,figh):
    t=dendropy.Tree.get(path=path,schema='newick',preserve_underscores=True)
    t.reroot_at_midpoint(update_bipartitions=True)
    # layout
    y={};x={}
    tips=[l for l in t.leaf_node_iter()]
    for i,l in enumerate(tips): y[l]=i
    for n in t.postorder_node_iter():
        if not n.is_leaf(): y[n]=sum(y[c] for c in n.child_node_iter())/len(n.child_nodes())
    h={}
    for n in t.postorder_node_iter():
        h[n]=0 if n.is_leaf() else 1+max(h[c] for c in n.child_node_iter())
    for n in t.preorder_node_iter(): x[n]=-h[n]
    x={k:v+max(h.values()) for k,v in x.items()}
    fig,ax=plt.subplots(figsize=(8,figh))
    for n in t.preorder_node_iter():
        if n.parent_node is not None:
            p=n.parent_node
            ax.plot([x[p],x[n]],[y[n],y[n]],color='#444',lw=0.8,solid_capstyle='butt')
        if not n.is_leaf():
            ys=[y[c] for c in n.child_node_iter()]
            ax.plot([x[n],x[n]],[min(ys),max(ys)],color='#444',lw=0.8)
            lab=n.label
            if lab:
                try:
                    v=float(lab.split('/')[-1])
                    if v>=70 and n is not t.seed_node: ax.text(x[n],y[n]-0.15,f'{v:.0f}',fontsize=4.8,ha='right',va='bottom',color='#555')
                except: pass
    xm=max(x.values())
    for l in tips:
        nm=l.taxon.label; isb=nm in B
        ax.text(x[l]+xm*0.012,y[l],nm,fontsize=6 if isb else 5,va='center',color=OR if isb else '#555',fontweight='bold' if isb else 'normal')
        if isb: ax.plot(x[l],y[l],'o',color=OR,ms=3.2)
    ax.set_xlim(-xm*0.02,xm*1.25); ax.set_ylim(len(tips),-1)
    ax.set_yticks([]); ax.set_xticks([]); ax.spines['bottom'].set_visible(False); ax.set_xlabel('Cladogram: branch lengths not to scale'); ax.set_title(title,fontsize=10,loc='left')
    ax.plot([],[],'o',color=OR,label='Black Bengal goat (n=22)'); ax.plot([],[],color=GR,label='Published goat mitogenomes')
    ax.legend(frameon=False,fontsize=7,loc='lower left',bbox_to_anchor=(0.0,0.02))
    ax.spines['left'].set_visible(False)
    plt.tight_layout(); plt.savefig(out,dpi=220); plt.close()
draw_tree(V+'trees_v6/world_nomoth.contree','figs/fig1a_world_tree.png','A. World set (101 taxa, 16,916 columns), ML tree (IQ-TREE 3.1.3), midpoint-rooted; node labels = bootstrap ≥70',15)
draw_tree(V+'trees_v6/complete_v2.contree','figs/fig1b_complete_tree.png','B. Complete mitogenomes (56 taxa, 16,919 columns), ML tree, midpoint-rooted',9)

# Fig 2 strict mismatch
txt=open('/root/.claude/uploads/2a9e653b-162e-5cc3-a048-85753dd8b966/890e776c-1791293587587_bbg_dloop_only22_v6.nex').read()
seq={l.split()[0]:l.split()[1] for l in txt.split('MATRIX')[1].split(';')[0].strip().splitlines() if l.strip()}
d=[sum(a!=b for a,b in zip(seq[x_],seq[y_])) for x_,y_ in itertools.combinations(seq,2)]
c=collections.Counter(d)
fig,ax=plt.subplots(figsize=(6,3.6))
ax.bar(list(range(0,17)),[c.get(k,0) for k in range(17)],color=BL,width=0.8)
ax.set_xlabel('Pairwise nucleotide differences (1,179 strict D-loop columns)'); ax.set_ylabel('Number of pairs')
ax.set_title('Mismatch distribution, 22 Black Bengal goats (231 pairs, mean 8.0)',fontsize=9,loc='left')
plt.tight_layout(); plt.savefig('figs/fig2_mismatch_strict.png',dpi=220); plt.close()

# Fig 4 phasing
rows={'116f':(.47,.26,.17,.10),'803m':(.46,.17,.09,.26),'Kh_5':(.44,.15,.11,.29),'246f':(.56,.23,.13,.07),'369m':(.64,.14,.15,.07),'584f':(.56,.18,.09,.17),'9ks':(.45,.32,.07,.16)}
labs=['TTGCC','TTATT','TTGTC','CCATT (HiFi / reference)']
cols=[GR,SK,YL,OR]
fig,ax=plt.subplots(figsize=(7,3.6))
names=list(rows)[::-1]
left=[0]*len(names)
for j in range(4):
    vals=[rows[n][j] for n in names]
    ax.barh(names,vals,left=left,color=cols[j],height=0.62,edgecolor='white',linewidth=1.2,label=labs[j])
    for i,(v,l0) in enumerate(zip(vals,left)):
        if v>=0.1: ax.text(l0+v/2,i,f'{v*100:.0f}',ha='center',va='center',fontsize=7,color='white' if j in (0,3) else 'black')
    left=[a+b for a,b in zip(left,vals)]
ax.set_xlim(0,1.0); ax.set_xlabel('Fraction of reads spanning positions 100–201'); ax.set_xticks([0,.25,.5,.75,1]); ax.set_xticklabels(['0','25%','50%','75%','100%'])
ax.set_title('Four read haplotypes at constant proportions; the mitochondrial one (CCATT) is a minority',fontsize=8.5,loc='left')
ax.legend(frameon=False,fontsize=7,ncol=4,loc='upper center',bbox_to_anchor=(0.5,-0.22))
plt.tight_layout(); plt.savefig('figs/fig4_phasing.png',dpi=220); plt.close()
print('ok')
