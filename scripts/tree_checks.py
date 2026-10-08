import dendropy, sys, collections
bbg=[l[1:].split()[0] for l in open('mt_remap/bbg_22S_mt_consensus_v6.fasta') if l.startswith('>')]
for name in ('complete_v2','dloop_v2','world_nomoth'):
    for kind in ('treefile','contree'):
        t=dendropy.Tree.get(path=f'mt_remap/trees_v6/{name}.{kind}',schema='newick',preserve_underscores=True)
        tips=[l.taxon.label for l in t.leaf_node_iter()]
        b=[x for x in tips if x in bbg]
        t.is_rooted=False
        t.encode_bipartitions()
        bb=set(t.taxon_namespace.get_taxon(x) for x in b)
        mask=sum(1<<i for i,tx in enumerate(t.taxon_namespace) if tx in bb)
        # split check: is there an edge whose split == BBG set (either side)?
        full=(1<<len(t.taxon_namespace))-1
        ok=any((e.bipartition.split_bitmask in (mask, full^mask)) for e in t.postorder_edge_iter() if e.bipartition)
        print(f'{name}.{kind}: tips={len(tips)} BBG={len(b)} BBG-forms-a-split={ok}')
    # midpoint-rooted (treefile): BBG clustering
t=dendropy.Tree.get(path='mt_remap/trees_v6/world_nomoth.contree',schema='newick',preserve_underscores=True)
t.reroot_at_midpoint(update_bipartitions=True)
bbgset=set(bbg)
# find maximal clades that are pure BBG
pure=[]
for n in t.preorder_node_iter():
    lv=[l.taxon.label for l in n.leaf_iter()]
    if lv and all(x in bbgset for x in lv):
        if n.parent_node is None or not all(x in bbgset for x in [l.taxon.label for l in n.parent_node.leaf_iter()]):
            pure.append((len(lv),n.label,lv))
print('midpoint-rooted world contree: pure-BBG maximal clades:',[(a,b) for a,b,_ in pure])
# nearest non-BBG reference (patristic) for every BBG on ML treefile, plus clade composition
t=dendropy.Tree.get(path='mt_remap/trees_v6/world_nomoth.treefile',schema='newick',preserve_underscores=True)
pdm=t.phylogenetic_distance_matrix()
tx={x.label:x for x in t.taxon_namespace}
res=[]
for b in bbg:
    ds=sorted((pdm.patristic_distance(tx[b],tx[o]),o) for o in tx if o not in bbgset)
    res.append((b,ds[0][1],round(ds[0][0],5)))
cnt=collections.Counter(r[1] for r in res)
print('nearest ref per BBG:',dict(cnt))
print('range of nearest dist:',min(r[2] for r in res),max(r[2] for r in res))
bd=[pdm.patristic_distance(tx[a],tx[b]) for i,a in enumerate(bbg) for b in bbg[i+1:]]
print('BBG-BBG patristic: mean %.5f max %.5f'%(sum(bd)/len(bd),max(bd)))
