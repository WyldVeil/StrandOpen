#!/usr/bin/env python3
"""Step 6: the maternal and paternal lines.

mtDNA: PhyloTree 17 (Forensic Update, rCRS coordinates). Each haplogroup's
expected profile is built down the tree from rCRS (back mutations, marked
'!', undo the last change at their position); the sample's profile is its
differences from rCRS. Candidates are scored the way Haplogrep does, with a
weighted Kulczynski measure (weights: PhyloTree's phylogenetic weights),
hotspots left out.

Y chromosome: the YFull tree, each branch's SNPs placed on GRCh38 through
YBrowse (name -> position, ancestral and derived allele). From the root,
the walk takes the child branch with the most derived SNPs while derived
calls outnumber ancestral ones.

Tables: hap_mt (path + candidates), hap_y (path), meta mt_/y_ keys."""
import json
import re
import xml.etree.ElementTree as ET
from strandlib import *
from genotype import Genotyper

con = db()
cur = con.cursor()
G = Genotyper(con)
lens = chrom_lengths()


def refbase(chrom, pos):
    """reference base(s) from the FASTA via samtools faidx"""
    return refseq(chrom, pos, pos)


_fai = {}
with open(FA + ".fai") as f:
    for line in f:
        p = line.split("\t")
        _fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
_fa = open(FA, "rb")


def refseq(chrom, a, b):
    ln, off, lb, lw = _fai[chrom]
    out = []
    for p in range(a, b + 1):
        z = p - 1
        _fa.seek(off + z // lb * lw + z % lb)
        out.append(_fa.read(1).decode().upper())
    return "".join(out)


# ======================== mtDNA ========================
PT = os.path.join(RAW, "phylotree")
mt_ref = refseq("MT", 1, 16569)
HOT = {"315.1C", "309.1C", "309.1CC", "523d", "524d", "524.1AC", "524.1ACAC", "3107d", "16182C", "16183C",
       "16193.1C", "16193.1CC", "16519C", "309.2C", "315.2C", "16193.2C", "310C"}
weights = {}
with open(os.path.join(PT, "weights.txt")) as f:
    for line in f:
        p = line.split()
        if len(p) >= 2:
            try:
                weights[p[0]] = float(p[1])
            except ValueError:
                pass


def wt(poly):
    return weights.get(poly, 1.0)


tree = ET.parse(os.path.join(PT, "tree.xml")).getroot()
nodes = []   # (name, parent index, polys)


def walk(el, parent):
    for hg in el.findall("haplogroup"):
        polys = [p.text.strip() for p in hg.findall("details/poly") if p.text]
        idx = len(nodes)
        nodes.append((hg.get("name"), parent, polys))
        walk(hg, idx)


walk(tree, -1)
log("PhyloTree nodes", len(nodes))

POLY = re.compile(r"^(\d+)(\.\d+)?([ACGTdNRYKMSWBDHV]+)(!*)$")


def apply(state, hist, poly):
    """state: {key: value} differences from rCRS; key is 'pos' or 'pos.k'"""
    m = POLY.match(poly.replace("(", "").replace(")", ""))
    if not m:
        return
    pos, ins, base, bang = m.groups()
    key = pos + (ins or "")
    if bang:
        h = hist.get(key)
        if h:
            h.pop()
            if h:
                state[key] = h[-1]
            else:
                state.pop(key, None)
        return
    hist.setdefault(key, []).append(base)
    if not ins and base == mt_ref[int(pos) - 1]:
        state.pop(key, None)
    else:
        state[key] = base


profiles = []
for i, (name, par, polys) in enumerate(nodes):
    if par < 0:
        st, hi = {}, {}
    else:
        st = dict(profiles[par][0])
        hi = {k: list(v) for k, v in profiles[par][1].items()}
    for p in polys:
        apply(st, hi, p)
    profiles.append((st, hi))


def as_polys(state):
    out = set()
    for k, b in state.items():
        out.add(k + b)
    return out


# the sample's profile
sample = set()
for pos, ref, alt, zyg, ps, dp in con.execute("SELECT pos, ref, alt, zyg, pass, dp FROM variant WHERE chrom=25 ORDER BY pos"):
    if len(ref) == 1 and len(alt) == 1:
        sample.add("%d%s" % (pos, alt))
    elif len(alt) > len(ref) and alt.startswith(ref):
        for k, b in enumerate(alt[len(ref):]):
            sample.add("%d.%d%s" % (pos + len(ref) - 1, k + 1, b))
    elif len(ref) > len(alt) and ref.startswith(alt):
        for k in range(len(ref) - len(alt)):
            sample.add("%dd" % (pos + len(alt) + k))
sample_eff = {p for p in sample if p not in HOT}
log("mtDNA sample polymorphisms", sorted(sample, key=lambda s: int(re.match(r"\d+", s).group())))

# A position the data never tested says nothing either way: a chip reads a
# few thousand of the 16,569, so an expected mutation there is neither found
# nor missing (otherwise a branch defined by losing an untested mutation, such
# as V+@72, beats V itself). Sequencing tests nearly every position.
_tested = {}


def tested(poly):
    m = re.match(r"\d+", poly)
    pos = int(m.group()) if m else 0
    if pos not in _tested:
        _tested[pos] = G.at(25, pos)[0] != "nocall"
    return _tested[pos]


# depth from the mitochondrial MRCA (the file is rooted at rCRS's own branch)
_adj = {}
for i, (name, par, polys) in enumerate(nodes):
    if par >= 0:
        _adj.setdefault(i, []).append(par)
        _adj.setdefault(par, []).append(i)
_mrca = next((k for k, n in enumerate(nodes) if n[0] == "mt-MRCA"), 0)
hdepth = {_mrca: 0}
_queue = [_mrca]
for _i in _queue:
    for _j in _adj.get(_i, []):
        if _j not in hdepth:
            hdepth[_j] = hdepth[_i] + 1
            _queue.append(_j)

best = []
sw = sum(wt(p) for p in sample_eff)
for i, (name, par, polys) in enumerate(nodes):
    exp = {p for p in as_polys(profiles[i][0]) if p not in HOT and tested(p)}
    if not exp:
        continue
    found = exp & sample_eff
    fw = sum(wt(p) for p in found)
    ew = sum(wt(p) for p in exp)
    k = 0.5 * (fw / ew + (fw / sw if sw else 0))
    # equal scores (nothing tested tells them apart): the less specific branch
    best.append((round(k, 9), -hdepth.get(i, 0), i))
best.sort(reverse=True)
best = [(k, i) for k, _, i in best]
cur.executescript("""
DROP TABLE IF EXISTS hap_mt; DROP TABLE IF EXISTS hap_mt_cand; DROP TABLE IF EXISTS hap_y;
CREATE TABLE hap_mt(depth INTEGER, name TEXT, polys TEXT, found TEXT, missing TEXT);
CREATE TABLE hap_mt_cand(rank INTEGER, name TEXT, score REAL);
CREATE TABLE hap_y(depth INTEGER, name TEXT, nder INTEGER, nanc INTEGER, der TEXT, anc TEXT, formed TEXT, tmrca TEXT, isogg TEXT);
""")
if not best:
    # no mitochondrial calls (some VCFs leave MT out): the maternal line stays unknown
    log("no mitochondrial calls: maternal line not determined")
    for k in ("mt_haplogroup", "mt_quality", "mt_private", "mt_missing", "mt_sample"):
        set_meta(con, k, "")
else:
    top_k, top_i = best[0]
    log("mtDNA best", nodes[top_i][0], "%.3f" % top_k)

    # The rCRS edition of PhyloTree is rooted at rCRS's own haplogroup
    # (H2a2a1): the ancestors of H2a2a1 hang below it with their edges reversed.
    # The path shown is the human one, from the mitochondrial MRCA down to the
    # sample's haplogroup: up the rCRS rooting from the MRCA to the shared node,
    # then down to the haplogroup. On a reversed edge the defining mutations are
    # the parent's states at the child's positions.
    def up(i):
        out = []
        while i >= 0:
            out.append(i)
            i = nodes[i][1]
        return out


    def state_at(i, key):
        st = profiles[i][0]
        if key in st:
            return st[key]
        return None if "." in key else mt_ref[int(key) - 1]


    root = next((k for k, n in enumerate(nodes) if n[0] == "mt-MRCA"), None)
    a_up, b_up = up(root) if root is not None else [], up(top_i)
    bset = set(b_up)
    lca = next((k for k in a_up if k in bset), b_up[-1])
    canon = []   # (node, defining mutations as written for the canonical direction)
    seq_a = a_up[:a_up.index(lca) + 1] if root is not None else []
    for k in range(len(seq_a) - 1):          # from the MRCA upwards: edges reversed
        child, par = seq_a[k], seq_a[k + 1]
        muts = []
        for p in nodes[child][2]:
            m = POLY.match(p.replace("(", "").replace(")", "").rstrip("!"))
            if not m:
                continue
            key = m.group(1) + (m.group(2) or "")
            b = state_at(par, key)
            muts.append(key + (b if b else "d"))
        canon.append((child, muts if k else []))
        canon[-1] = (child, canon[-1][1])
    # the node order: MRCA ... LCA, then down to the haplogroup
    order = [seq_a[0]] if seq_a else []
    mut_of = {}
    for k in range(1, len(seq_a)):
        child = seq_a[k - 1]
        node = seq_a[k]
        ms = []
        for p in nodes[child][2]:
            m = POLY.match(p.replace("(", "").replace(")", "").rstrip("!"))
            if m:
                key = m.group(1) + (m.group(2) or "")
                b = state_at(node, key)
                ms.append(key + (b if b else "d"))
        mut_of[node] = ms
        order.append(node)
    down = b_up[:b_up.index(lca)][::-1]
    for node in down:
        mut_of[node] = [p.rstrip("!") for p in nodes[node][2]]
        order.append(node)
    if not order:
        order = b_up[::-1]


    def sample_has(mut):
        m = POLY.match(mut)
        if not m:
            return False
        pos, ins, base = m.group(1), m.group(2), m.group(3)
        if mut in sample:
            return True
        if not ins and base == mt_ref[int(pos) - 1]:   # the reference state: present unless the sample differs there
            return not any(re.match(r"^%s[ACGTd]" % pos, s) for s in sample)
        return False


    for d, i in enumerate(order):
        found, missing = [], []
        later = set()
        for j in order[d + 1:]:
            for q in mut_of.get(j, []):
                m = POLY.match(q)
                if m:
                    later.add(m.group(1))
        for q in mut_of.get(i, []):
            if q in HOT or not tested(q):
                continue
            m = POLY.match(q)
            if not sample_has(q) and m and m.group(1) in later:
                continue   # changed again further down the path: not expected now
            (found if sample_has(q) else missing).append(q)
        cur.execute("INSERT INTO hap_mt VALUES(?,?,?,?,?)", (d, nodes[i][0], " ".join(mut_of.get(i, [])), " ".join(found), " ".join(missing)))
    for r, (k, i) in enumerate(best[:12]):
        cur.execute("INSERT INTO hap_mt_cand VALUES(?,?,?)", (r + 1, nodes[i][0], k))
    exp_best = {p for p in as_polys(profiles[top_i][0]) if p not in HOT}
    private = sorted(sample_eff - exp_best, key=lambda s: int(re.match(r"\d+", s).group()))
    missing_all = sorted((p for p in exp_best - sample_eff if tested(p)), key=lambda s: int(re.match(r"\d+", s).group()))
    set_meta(con, "mt_haplogroup", nodes[top_i][0])
    set_meta(con, "mt_quality", "%.3f" % top_k)
    set_meta(con, "mt_private", " ".join(private))
    set_meta(con, "mt_missing", " ".join(missing_all))
    set_meta(con, "mt_sample", " ".join(sorted(sample, key=lambda s: int(re.match(r"\d+", s).group()))))
    set_meta(con, "mt_tree", "PhyloTree 17 Forensic Update (rCRS), via genepi/phylotree-fu-rcrs")

# ======================== Y chromosome ========================
yt = json.load(open(os.path.join(RAW, "yfull_tree.json")))
names_needed = set()


def names_of(snps):
    out = []
    for grp in snps.split(","):
        grp = grp.strip()
        if not grp:
            continue
        alts = [a.strip() for a in grp.split("/") if a.strip()]
        out.append(alts)
        names_needed.update(alts)
    return out


ynodes = []


def ywalk(n, parent, depth):
    idx = len(ynodes)
    ynodes.append(dict(id=n.get("id", ""), parent=parent, depth=depth, snps=names_of(n.get("snps", "")),
                       formed=str(n.get("formed", "")), tmrca=str(n.get("tmrca", "")), kids=[]))
    if parent >= 0:
        ynodes[parent]["kids"].append(idx)
    for c in n.get("children", []):
        ywalk(c, idx, depth + 1)


ywalk(yt, -1, 0)
log("YFull nodes", len(ynodes), "SNP names", len(names_needed))

ypos = {}   # name -> (pos, anc, der, isogg)
for line in opentext(os.path.join(RAW, "snps_hg38.vcf.gz")):
    if line[0] == "#":
        continue
    p = line.split("\t", 8)
    # one SNP is often listed under several names, comma-separated: each is looked up on its own
    names = [x for x in p[2].split(",") if x in names_needed]
    if not names or len(p[3]) != 1 or len(p[4]) != 1:
        continue
    info = p[7]
    m = re.search(r"AA=([ACGT])", info)
    anc = m.group(1) if m else p[3]
    der = p[4] if anc == p[3] else p[3]
    iso = re.search(r'ISOGG="([^"]*)"', info)
    for nm in names:
        ypos[nm] = (int(p[1]), anc, der, iso.group(1) if iso else "")
log("Y SNPs placed on GRCh38", len(ypos))

ycall = {}


def ystate(name):
    """'D' derived, 'A' ancestral, None no call / not placed"""
    if name in ycall:
        return ycall[name]
    v = ypos.get(name)
    s = None
    if v:
        pos, anc, der, _ = v
        al = G.alleles(24, pos, refseq("Y", pos, pos), min_dp=2)
        if al:
            a = al[0]
            s = "D" if a == der else "A" if a == anc else None
    ycall[name] = s
    return s


def node_score(i):
    d = a = 0
    dn, an = [], []
    for alts in ynodes[i]["snps"]:
        st = None
        nm = alts[0]
        for x in alts:
            s = ystate(x)
            if s:
                st, nm = s, x
                break
        if st == "D":
            d += 1
            dn.append(nm)
        elif st == "A":
            a += 1
            an.append(nm)
    return d, a, dn, an


# Many branches on a true path have no callable SNP at all (a few SNPs, all
# in unmappable sequence), so a walk that needs evidence at every step stops
# early. Instead every node gets a path score: its ancestors' derived calls
# minus twice their ancestral calls, plus its own. The terminal is the best
# scoring node that has derived support of its own (ties: the deeper one).
scores = [None] * len(ynodes)
pscore = [0.0] * len(ynodes)
best_i, best_s = 0, -1e9
for i, n in enumerate(ynodes):   # parents come before children (preorder)
    d, a, _, _ = node_score(i)
    scores[i] = (d, a)
    pscore[i] = (pscore[n["parent"]] if n["parent"] >= 0 else 0) + d - 2 * a
    if d >= 1 and d > a and (pscore[i] > best_s or (pscore[i] == best_s and n["depth"] > ynodes[best_i]["depth"])):
        best_i, best_s = i, pscore[i]
if best_s < -1e8:
    # no derived Y call anywhere: a woman, or a file without chromosome Y
    log("no Y chromosome calls: paternal line not determined")
    for k in ("y_haplogroup", "y_terminal_snps", "y_formed", "y_tmrca"):
        set_meta(con, k, "")
else:
    ypath = []
    i = best_i
    while i >= 0:
        ypath.append(i)
        i = ynodes[i]["parent"]
    ypath.reverse()
    for d, i in enumerate(ypath):
        n = ynodes[i]
        nd, na, dn, an = node_score(i)
        iso = ""
        for nm in dn:
            if ypos.get(nm, (0, "", "", ""))[3]:
                iso = ypos[nm][3]
                break
        cur.execute("INSERT INTO hap_y VALUES(?,?,?,?,?,?,?,?,?)", (d, n["id"] or "root", nd, na, " ".join(dn[:60]), " ".join(an[:60]),
                                                                    n["formed"], n["tmrca"], iso))
    term = ynodes[ypath[-1]]
    set_meta(con, "y_haplogroup", term["id"])
    nd, na, dn, an = node_score(ypath[-1])
    set_meta(con, "y_terminal_snps", " ".join(dn[:8]))
    set_meta(con, "y_formed", term["formed"])
    set_meta(con, "y_tmrca", term["tmrca"])
    set_meta(con, "y_tree", "YFull YTree (current_tree.json), SNP positions from YBrowse (GRCh38)")
    log("Y haplogroup", term["id"], "derived", nd, "ancestral", na, "depth", len(ypath))
con.commit()
