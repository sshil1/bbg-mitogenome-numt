import sys, subprocess, re, collections
bam = sys.argv[1]
SH, L = 8000, 16643
sites = [100, 105, 152, 193, 201]            # v3 coordinates
rot = [(p - 1 - SH) % L + 1 for p in sites]
lo, hi = min(rot), max(rot)
out = subprocess.run(["samtools", "view", bam, f"NC_005044.2:{lo}-{hi}"],
                     capture_output=True, text=True).stdout
haps = collections.defaultdict(list)
for line in out.splitlines():
    f = line.split("\t")
    flag, pos, mapq, cigar, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
    if flag & 0x904:
        continue
    r, q, base = pos, 0, {}
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                base[r + k] = seq[q + k]
            r += n; q += n
        elif op in "IS":
            q += n
        elif op in "DN":
            r += n
    if all(p in base for p in rot):
        haps["".join(base[p] for p in rot)].append(mapq)
tot = sum(len(v) for v in haps.values())
print(f"{bam.split('/')[-1]}  reads spanning all 5 sites: {tot}  (sites {sites})")
for h, v in sorted(haps.items(), key=lambda x: -len(x[1]))[:6]:
    print(f"  {h}  n={len(v):4d}  frac={len(v)/tot:.2f}  meanMAPQ={sum(v)/len(v):.1f}")
