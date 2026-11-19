#!/usr/bin/env python3
"""Step 2: genes, their canonical transcript structure, what they do and
which conditions they are linked to.

Tables: chrom, cytoband, gene, exon, gene_disease, gene_fts."""
import re
from strandlib import *

con = db()
cur = con.cursor()
lens = chrom_lengths()

cur.executescript("""
DROP TABLE IF EXISTS chrom; DROP TABLE IF EXISTS cytoband; DROP TABLE IF EXISTS gene;
DROP TABLE IF EXISTS exon; DROP TABLE IF EXISTS gene_disease; DROP TABLE IF EXISTS gene_fts;
CREATE TABLE chrom(idx INTEGER PRIMARY KEY, name TEXT, length INTEGER);
CREATE TABLE cytoband(chrom INTEGER, start INTEGER, end INTEGER, name TEXT, stain TEXT);
CREATE TABLE gene(id INTEGER PRIMARY KEY, ensg TEXT, symbol TEXT, name TEXT, biotype TEXT,
  chrom INTEGER, start INTEGER, end INTEGER, strand INTEGER, tx TEXT, mane INTEGER,
  cds_start INTEGER, cds_end INTEGER, nexon INTEGER, hgnc TEXT,
  protein TEXT, function TEXT, uniprot_disease TEXT, uniprot TEXT,
  pli REAL, loeuf REAL, mis_z REAL);
CREATE TABLE exon(gene INTEGER, start INTEGER, end INTEGER, kind INTEGER);
CREATE TABLE gene_disease(gene INTEGER, disease_id TEXT, disease TEXT, inheritance TEXT);
""")
for c in CHROMS:
    cur.execute("INSERT INTO chrom VALUES(?,?,?)", (CHROM_IDX[c], c, lens[c]))

with opentext(os.path.join(RAW, "cytoBand.txt.gz")) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        c = p[0][3:] if p[0].startswith("chr") else p[0]
        if c == "M":
            c = "MT"
        if c in CHROM_IDX:
            cur.execute("INSERT INTO cytoband VALUES(?,?,?,?,?)", (CHROM_IDX[c], int(p[1]) + 1, int(p[2]), p[3], p[4]))
log("chromosomes and cytobands")

# ---- GFF3: genes, the canonical transcript of each, its exons and CDS
genes = {}      # ensg -> dict
txgene = {}     # canonical enst -> ensg
txparts = {}    # enst -> list of (start, end, kind) kind: 0 exon, 1 CDS
gff = os.path.join(RAW, "Homo_sapiens.GRCh38.115.gff3.gz")


def attrs(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


with opentext(gff) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[0] not in CHROM_IDX:
            continue
        typ = p[2]
        if "gene" in typ and typ not in ("pseudogenic_transcript",) and p[8].startswith("ID=gene:"):
            a = attrs(p[8])
            ensg = a["ID"][5:]
            desc = unquote_gff(a.get("description", ""))
            hgnc = ""
            m = re.search(r"\[Source:([^\]]*)\]", desc)
            if m:
                mm = re.search(r"Acc:(HGNC:\d+)", m.group(1))
                hgnc = mm.group(1) if mm else ""
                desc = desc[:m.start()].strip()
            genes[ensg] = dict(ensg=ensg, symbol=a.get("Name", ensg), name=desc, biotype=a.get("biotype", typ),
                               chrom=CHROM_IDX[p[0]], start=int(p[3]), end=int(p[4]),
                               strand=1 if p[6] == "+" else -1, tx="", mane=0, hgnc=hgnc)
        elif p[8].startswith("ID=transcript:"):
            a = attrs(p[8])
            tags = a.get("tag", "")
            if "Ensembl_canonical" in tags:
                par = a.get("Parent", "")[5:]
                enst = a["ID"][11:]
                txgene[enst] = par
                if par in genes:
                    genes[par]["tx"] = enst
                    genes[par]["mane"] = 1 if "MANE_Select" in tags else 0
        elif typ in ("exon", "CDS"):
            par = p[8].split("Parent=transcript:", 1)
            if len(par) == 2:
                enst = par[1].split(";", 1)[0]
                if enst in txgene:
                    txparts.setdefault(enst, []).append((int(p[3]), int(p[4]), 0 if typ == "exon" else 1))
log("genes", len(genes), "canonical transcripts", len(txgene))

# ---- UniProt (reviewed human): protein name, function, disease involvement
uni = {}
with opentext(os.path.join(RAW, "uniprot_human.tsv.gz")) as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5 or not p[1]:
            continue
        func = re.sub(r"\s*\{ECO:[^}]*\}", "", p[3]).replace("FUNCTION: ", "").strip()
        func = re.sub(r"\(PubMed:[0-9, PubMed:]+\)", "", func)
        func = re.sub(r"\s+", " ", func).strip()
        func = re.sub(r"\s+([.,;:)])", r"\1", func).replace("( ", "(").replace("..", ".")
        dis = re.sub(r"\s*\{ECO:[^}]*\}", "", p[4]).replace("DISEASE: ", "").strip()
        dis = re.sub(r"\[MIM:\d+\]:?", "", dis)
        dis = re.sub(r"\(PubMed:[0-9, PubMed:]+\)", "", dis)
        dis = re.sub(r"\s+", " ", dis).strip()
        dis = re.sub(r"\s+([.,;:)])", r"\1", dis).replace("( ", "(").replace("..", ".")
        prot = re.sub(r"\s*\(EC [^)]*\)", "", p[2]).split(" (")[0]
        if p[1] not in uni:
            uni[p[1]] = (p[0], prot, func, dis)
log("uniprot entries", len(uni))

# ---- gnomAD v4.1 constraint: Ensembl canonical rows first
cons = {}
with open(os.path.join(RAW, "gnomad.v4.1.constraint_metrics.tsv")) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}

    def num(v):
        try:
            return float(v)
        except ValueError:
            return None
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[ix["canonical"]] != "true":
            continue
        gid = p[ix["gene_id"]]
        key = gid if gid.startswith("ENSG") else p[ix["gene"]]
        row = (num(p[ix["lof.pLI"]]), num(p[ix["lof.oe_ci.upper"]]), num(p[ix["mis.z_score"]]))
        if key not in cons or gid.startswith("ENSG"):
            cons[key] = row
log("constraint rows", len(cons))

rows = []
gid_of = {}
for i, g in enumerate(sorted(genes.values(), key=lambda g: (g["chrom"], g["start"]))):
    gid = i + 1
    gid_of[g["ensg"]] = gid
    parts = txparts.get(g["tx"], [])
    cds = [x for x in parts if x[2] == 1]
    u = uni.get(g["symbol"], ("", "", "", ""))
    c = cons.get(g["ensg"]) or cons.get(g["symbol"]) or (None, None, None)
    rows.append((gid, g["ensg"], g["symbol"], g["name"], g["biotype"], g["chrom"], g["start"], g["end"], g["strand"],
                 g["tx"], g["mane"], min((x[0] for x in cds), default=None), max((x[1] for x in cds), default=None),
                 sum(1 for x in parts if x[2] == 0), g["hgnc"], u[1], u[2], u[3], u[0], c[0], c[1], c[2]))
    for s, e, k in sorted(parts):
        cur.execute("INSERT INTO exon VALUES(?,?,?,?)", (gid, s, e, k))
cur.executemany("INSERT INTO gene VALUES(" + ",".join("?" * 22) + ")", rows)
cur.execute("CREATE INDEX gene_pos ON gene(chrom, start)")
cur.execute("CREATE INDEX gene_sym ON gene(symbol COLLATE NOCASE)")
cur.execute("CREATE INDEX exon_gene ON exon(gene)")
log("gene rows", len(rows))

# ---- HPO: gene -> disease, with the disease's mode of inheritance
inh_name = {"HP:0000006": "AD", "HP:0000007": "AR", "HP:0001417": "XL", "HP:0001419": "XLR",
            "HP:0001423": "XLD", "HP:0001427": "MT", "HP:0001426": "MF", "HP:0001442": "SOM",
            "HP:0001450": "YL", "HP:0010984": "DIG", "HP:0032113": "SDOM"}
dname, dinh = {}, {}
with open(os.path.join(RAW, "phenotype.hpoa")) as f:
    for line in f:
        if line.startswith("#") or line.startswith("database_id"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 11:
            continue
        dname.setdefault(p[0], p[1])
        if p[10] == "I" and p[3] in inh_name and not p[2]:
            dinh.setdefault(p[0], set()).add(inh_name[p[3]])
sym_gid = {}
for g in rows:
    sym_gid.setdefault(g[2], g[0])
n = 0
seen = set()
with open(os.path.join(RAW, "genes_to_disease.txt")) as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or p[2] != "MENDELIAN" or p[1] not in sym_gid:
            continue
        key = (p[1], p[3])
        if key in seen:
            continue
        seen.add(key)
        cur.execute("INSERT INTO gene_disease VALUES(?,?,?,?)",
                    (sym_gid[p[1]], p[3], dname.get(p[3], p[3]), ",".join(sorted(dinh.get(p[3], [])))))
        n += 1
cur.execute("CREATE INDEX gd_gene ON gene_disease(gene)")
log("gene-disease links", n)

cur.execute("CREATE VIRTUAL TABLE gene_fts USING fts5(symbol, name, protein, content='gene', content_rowid='id')")
cur.execute("INSERT INTO gene_fts(rowid, symbol, name, protein) SELECT id, symbol, name, protein FROM gene")
set_meta(con, "genes_source", "Ensembl 115 (GRCh38), UniProt reviewed human, gnomAD v4.1 constraint, HPO")
con.commit()
log("done")
