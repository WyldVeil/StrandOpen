#!/usr/bin/env python3
"""Step 19: where your two lineages come from.

Your paternal (Y) and maternal (mtDNA) haplogroups from step 6 are placed
next to ancient people from the Allen Ancient DNA Resource (AADR,
ref/raw/aadr/*.anno: dated, located individuals with their own Y and mtDNA
haplogroups). For each ancient person the deepest branch shared with you is
found on the same trees step 6 uses (YFull for Y, PhyloTree for mtDNA), so
"shares your branch R-U106" means their line and yours split at that branch.

Tables:
  lin_ancient  ancient people on your lines: where, when, culture, their
               haplogroup and the branch shared with you
  lin_node     the branches on your path with their ages and the earliest
               ancient carriers found under each
  lin_story    a short narrative per line, built from the above and from
               data/haplogroup_origins.tsv (what is known of the big clades)
"""
import csv
import glob
import json
import re
import xml.etree.ElementTree as ET
from strandlib import *

con = db()
cur = con.cursor()
HERE = os.path.dirname(os.path.abspath(__file__))
anno = sorted(glob.glob(os.path.join(RAW, "aadr", "*1240K*.anno")))
if not anno:
    log("no AADR annotation file in ref/raw/aadr: skipped")
    raise SystemExit(0)
ver = os.path.basename(anno[-1]).split("_")[0]

# ---- the two trees
yt = json.load(open(os.path.join(RAW, "yfull_tree.json")))
yparent, yformed, ytmrca, ysnp = {}, {}, {}, {}


anon = [0]


def ywalk(n, parent):
    nid = n.get("id", "")
    if not nid:
        nid = "root" if parent is None else None
    if nid is None or nid in yparent:          # unnamed or repeated ids get their own key
        anon[0] += 1
        nid = "_node%d" % anon[0]
    yparent[nid] = parent
    yformed[nid], ytmrca[nid] = str(n.get("formed", "")), str(n.get("tmrca", ""))
    for grp in (n.get("snps", "") or "").split(","):
        for a in grp.split("/"):
            a = a.strip()
            if a and a not in ysnp:
                ysnp[a] = nid
    for c in n.get("children", []):
        ywalk(c, nid)


ywalk(yt, None)
# PhyloTree's rCRS edition is rooted at the reference sequence (H2a2a1), so
# along that stretch its parent links point the wrong way: re-root at mt-MRCA
madj = {}


def mwalk(el, parent):
    for hg in el.findall("haplogroup"):
        name = hg.get("name")
        if parent:
            madj.setdefault(parent, set()).add(name)
            madj.setdefault(name, set()).add(parent)
        madj.setdefault(name, set())
        mwalk(hg, name)


mwalk(ET.parse(os.path.join(RAW, "phylotree", "tree.xml")).getroot(), None)
mparent = {"mt-MRCA": None}
todo = ["mt-MRCA"]
while todo:
    n = todo.pop()
    for k in madj.get(n, ()):
        if k not in mparent:
            mparent[k] = n
            todo.append(k)


def chain(node, parent):
    out, seen = [], set()
    while node and node not in seen:
        seen.add(node)
        out.append(node)
        node = parent.get(node)
    return out


def y_node(call):
    """'R-U106' (terminal SNP notation) -> YFull node id"""
    call = call.strip().rstrip("~*")
    if call in yparent:
        return call
    m = re.match(r"^([A-Z][A-Z0-9]*)-(\S+)$", call)
    if m:
        snp = m.group(2).rstrip("~*")
        if snp in ysnp:
            return ysnp[snp]
    return None


def mt_node(call):
    """'K1a+195', 'H1c3', 'U5a1a1' -> a PhyloTree node, trimming the name until one is found"""
    call = call.strip().split(" ")[0]
    if not call or call.startswith("n/a") or call == "..":
        return None
    for c in (call, call.split("+")[0], call.rstrip("*'")):
        if c in mparent:
            return c
    c = call.split("+")[0]
    while len(c) > 1:
        c = c[:-1]
        if c in mparent:
            return c
    return None


# ---- your paths (step 6)
ypath = [r[0] for r in cur.execute("SELECT name FROM hap_y ORDER BY depth")] if \
    cur.execute("SELECT 1 FROM sqlite_master WHERE name='hap_y'").fetchone() else []
mpath = [r[0] for r in cur.execute("SELECT name FROM hap_mt ORDER BY depth")] if \
    cur.execute("SELECT 1 FROM sqlite_master WHERE name='hap_mt'").fetchone() else []
ypath = ["root" if n == "root" else n for n in ypath]
yterm = y_node(ypath[-1]) if len(ypath) > 1 else None
mterm = mpath[-1] if mpath else None
yset = {n: i for i, n in enumerate(chain(yterm, yparent)[::-1])} if yterm else {}
mset = {n: i for i, n in enumerate(chain(mterm, mparent)[::-1])} if mterm and mterm in mparent else {}
log("your lines:", yterm or "no Y", "/", mterm or "no mtDNA")

PERIOD = [("Meso", "Mesolithic"), ("EN", "Early Neolithic"), ("MN", "Middle Neolithic"), ("LN", "Late Neolithic"), ("N", "Neolithic"),
          ("CA", "Copper Age"), ("Eneolithic", "Copper Age"), ("EBA", "Early Bronze Age"), ("MBA", "Middle Bronze Age"),
          ("LBA", "Late Bronze Age"), ("BA", "Bronze Age"), ("EIA", "Early Iron Age"), ("LIA", "Late Iron Age"), ("IA", "Iron Age"),
          ("Viking", "Viking Age"), ("VA", "Viking Age"), ("Roman", "Roman period"), ("Medieval", "medieval"), ("EMedieval", "early medieval"),
          ("Saxon", "Anglo-Saxon"), ("BellBeaker", "Bell Beaker culture"), ("CordedWare", "Corded Ware culture"), ("Yamnaya", "Yamnaya culture"),
          ("UP", "Upper Palaeolithic"), ("HG", "hunter-gatherer"), ("WHG", "western hunter-gatherer")]


def culture(group):
    words = []
    for part in group.split("_")[1:]:
        if part in ("o", "o1", "o2", "contam", "lc", "dup", "rel", "brother", "sister", "son", "daughter") or re.match(r"^o\d*$", part):
            break
        lab = dict(PERIOD).get(part)
        if lab is None and len(part) <= 2 and part.isupper():
            continue          # site letters such as _C, _MA
        words.append(lab or re.sub(r"(?<=[a-z])(?=[A-Z])", " ", part))
    return ", ".join(w for w in words if w)[:60]


def years_ago(bp):
    y = bp + 75                # BP counts from 1950
    if y < 1000:
        return int(round(y, -1))
    return int(round(y, -2))


# ---- the ancient people
rows = list(csv.reader(open(anno[-1], encoding="utf-8", errors="replace"), delimiter="\t"))[1:]
seen = set()
anc = []
for r in rows:
    if len(r) < 39:
        continue
    try:
        bp = int(float(r[10]))
        lat, lon = float(r[17]), float(r[18])
    except ValueError:
        continue
    if bp <= 0:
        continue                # present-day people
    ind = r[2]
    if ind in seen:
        continue
    seen.add(ind)
    anc.append(dict(id=ind, bp=bp, date=r[12], group=r[14], place=r[15], country=r[16], lat=lat, lon=lon, sex=r[30],
                    y=r[34].strip(), mt=r[38].strip(), pub=r[5], doi=r[7].strip()))
log("ancient people with dates and places", len(anc))

cur.executescript("""DROP TABLE IF EXISTS lin_ancient; DROP TABLE IF EXISTS lin_node; DROP TABLE IF EXISTS lin_story;
CREATE TABLE lin_ancient(line TEXT, ind TEXT, culture TEXT, grp TEXT, place TEXT, country TEXT, lat REAL, lon REAL, years INTEGER,
  date_text TEXT, hg TEXT, shared TEXT, shared_depth INTEGER, gap INTEGER, pub TEXT, doi TEXT);
CREATE TABLE lin_node(line TEXT, depth INTEGER, name TEXT, formed TEXT, tmrca TEXT, n_ancient INTEGER, earliest_years INTEGER,
  earliest_place TEXT, earliest_culture TEXT);
CREATE TABLE lin_story(line TEXT, ord INTEGER, text TEXT);""")

def major_depth(pathset):
    """the most specific one-letter clade on your path (I, R, H, V, U...):
    ancient people are listed from there down, not back to the root of humanity"""
    best = 0
    for name, d in pathset.items():
        if re.match(r"^[A-Z]$", name):
            best = max(best, d)
    return best


for line, pathset, parent, term, getnode, field in (("Y", yset, yparent, yterm, y_node, "y"), ("MT", mset, mparent, mterm, mt_node, "mt")):
    if not pathset:
        continue
    depth_term = len(pathset) - 1
    major = major_depth(pathset)
    set_meta(con, "lin_%s_major" % line.lower(), [n for n, d in pathset.items() if d == major][0])
    n = 0
    for a in anc:
        node = getnode(a[field]) if a[field] and not a[field].startswith("n/a") else None
        if not node:
            continue
        shared, sd = None, -1
        for x in chain(node, parent):
            if x in pathset:
                shared, sd = x, pathset[x]
                break
        if shared is None or sd < major:
            continue
        cur.execute("INSERT INTO lin_ancient VALUES(" + ",".join("?" * 16) + ")",
                    (line, a["id"], culture(a["group"]), a["group"], a["place"], a["country"], a["lat"], a["lon"], years_ago(a["bp"]),
                     a["date"], a[field], shared, sd, depth_term - sd, a["pub"], a["doi"]))
        n += 1
    log(line, "ancient people sharing a branch with you", n)
    # each branch on your path: ages, and the earliest ancient carriers at or below it
    for name, d in sorted(pathset.items(), key=lambda t: t[1]):
        if d < major:
            continue
        r = cur.execute("SELECT count(*), max(years) FROM lin_ancient WHERE line=? AND shared_depth>=?", (line, d)).fetchone()
        e = cur.execute("SELECT place, country, culture FROM lin_ancient WHERE line=? AND shared_depth>=? ORDER BY years DESC LIMIT 1",
                        (line, d)).fetchone()
        cur.execute("INSERT INTO lin_node VALUES(?,?,?,?,?,?,?,?,?)",
                    (line, d, name, yformed.get(name, "") if line == "Y" else "", ytmrca.get(name, "") if line == "Y" else "",
                     r[0], r[1], ("%s, %s" % (e[0], e[1])) if e else "", e[2] if e else ""))

# ---- the stories
origins = {}
of = os.path.join(HERE, "..", "data", "haplogroup_origins.tsv")
if os.path.exists(of):
    for ln in open(of, encoding="ascii"):
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 3 and not p[0].startswith("#"):
            origins[(p[0], p[1])] = p[2]


def fmt(n):
    return format(int(n), ",")


for line, pathset in (("Y", yset), ("MT", mset)):
    if not pathset:
        continue
    path = [n for n, _ in sorted(pathset.items(), key=lambda t: t[1])]
    out = []
    # what is known of the big clades on the path, broadest first
    for name in path:
        t = origins.get((line, name))
        if t:
            out.append(t)
    term = path[-1]
    if line == "Y":
        f, t = yformed.get(term, ""), ytmrca.get(term, "")
        if f.isdigit():
            out.append("Your own branch, %s, split from its nearest relatives about %s years ago%s." % (
                term, fmt(f), "" if not t.isdigit() or t == f else ", and the men known to carry it share a common paternal ancestor who lived about %s years ago" % fmt(t)))
    # the closest ancient people
    best = cur.execute("SELECT shared, shared_depth, count(*), max(years), min(years) FROM lin_ancient WHERE line=? "
                       "GROUP BY shared ORDER BY shared_depth DESC LIMIT 1", (line,)).fetchone()
    if best:
        ex = cur.execute("SELECT culture, place, country, years FROM lin_ancient WHERE line=? AND shared=? ORDER BY years DESC LIMIT 1",
                         (line, best[0])).fetchone()
        who = "man" if line == "Y" else "person"
        out.append("The ancient %s closest to your %s line in the database (%s ancient genomes searched) shared your branch %s: %s"
                   % ("people" if best[2] > 1 else who, "paternal" if line == "Y" else "maternal", fmt(len(anc)), best[0],
                      ("%d people, the oldest" % best[2] if best[2] > 1 else "one")
                      + " from %s, %s (%s) about %s years ago." % (ex[1], ex[2], ex[0] or "unknown culture", fmt(ex[3]))))
    for i, t in enumerate(out):
        cur.execute("INSERT INTO lin_story VALUES(?,?,?)", (line, i, t))

set_meta(con, "lineage_source", "Allen Ancient DNA Resource %s (%s ancient individuals with dates and places)" % (ver, format(len(anc), ",")))
con.commit()
log("done")
