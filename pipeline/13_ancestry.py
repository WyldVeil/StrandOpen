#!/usr/bin/env python3
"""Step 13: ancestry against the 1000 Genomes reference (step 4).

1. A panel of common, independent SNPs (MAF >= 5%, LD-pruned, the long-range
   LD regions left out), exported from the founders' genotypes.
2. PCA of the 2,590 founders; the sample projected onto it. A second PCA
   of the European founders alone, for the finer picture inside Europe.
3. Supervised admixture: the sample as a mixture of 21 reference
   populations (the admixed ones, ASW ACB MXL PUR CLM, left out), by EM on
   the likelihood of its genotypes; 95% intervals from a bootstrap over
   5 Mb blocks.
4. Chromosome painting: windows of 250 SNPs, each given the pair of
   continental groups that best explains it, smoothed by a Viterbi pass
   that charges for a switch.
5. Closest populations: mean log-likelihood per SNP under each population.

Tables: anc_pop, anc_admix, anc_pca, anc_paint, anc_closest."""
import subprocess
from shlex import quote as shq
import numpy as np
from strandlib import *
from genotype import Genotyper

P = KG
A = os.path.join(WORK, "anc")
os.makedirs(A, exist_ok=True)
con = db()
cur = con.cursor()
G = Genotyper(con)

LABEL = {
    "GBR": ("British (England and Scotland)", "EUR"), "CEU": ("Northwest European (Utah, CEPH)", "EUR"), "FIN": ("Finnish", "EUR"),
    "IBS": ("Spanish (Iberian)", "EUR"), "TSI": ("Italian (Tuscany)", "EUR"),
    "YRI": ("Yoruba (Nigeria)", "AFR"), "LWK": ("Luhya (Kenya)", "AFR"), "GWD": ("Gambian Mandinka", "AFR"), "MSL": ("Mende (Sierra Leone)", "AFR"),
    "ESN": ("Esan (Nigeria)", "AFR"), "ASW": ("African American (US)", "AFR"), "ACB": ("African Caribbean (Barbados)", "AFR"),
    "CHB": ("Han Chinese (Beijing)", "EAS"), "JPT": ("Japanese (Tokyo)", "EAS"), "CHS": ("Han Chinese (South)", "EAS"),
    "CDX": ("Dai Chinese (Xishuangbanna)", "EAS"), "KHV": ("Kinh (Vietnam)", "EAS"),
    "GIH": ("Gujarati Indian", "SAS"), "PJL": ("Punjabi (Pakistan)", "SAS"), "BEB": ("Bengali (Bangladesh)", "SAS"),
    "STU": ("Sri Lankan Tamil", "SAS"), "ITU": ("Telugu Indian", "SAS"),
    "MXL": ("Mexican (Los Angeles)", "AMR"), "PUR": ("Puerto Rican", "AMR"), "CLM": ("Colombian (Medellin)", "AMR"), "PEL": ("Peruvian (Lima)", "AMR"),
}
SUPER = {"EUR": "European", "AFR": "African", "EAS": "East Asian", "SAS": "South Asian", "AMR": "Native American"}
ADMIXED = {"ASW", "ACB", "MXL", "PUR", "CLM"}
LONG_LD = [(6, 25_000_000, 35_000_000), (8, 8_000_000, 12_500_000), (17, 45_000_000, 47_000_000), (2, 134_000_000, 138_500_000),
           (5, 44_000_000, 52_000_000), (11, 87_000_000, 91_000_000), (12, 33_000_000, 41_000_000), (3, 47_000_000, 50_000_000),
           (14, 47_000_000, 50_000_000), (20, 33_000_000, 36_000_000)]

# ---- 1. the panel
rng = os.path.join(A, "longld.txt")
with open(rng, "w") as f:
    for c, a, b in LONG_LD:
        f.write("%d\t%d\t%d\tld\n" % (c, a, b))
# A chip reads a fixed few hundred thousand positions, and a panel pruned
# without regard to them keeps only about one in seven where the chip has a
# call. For chip data the panel is pruned from the positions it read, which
# gives several times the markers and much tighter estimates.
_cm = cur.execute("SELECT value FROM meta WHERE key='call_model'").fetchone()
CHIP = bool(_cm) and _cm[0] == "chip"
tag, extract = "", ""
if CHIP:
    tag = "_chip"
    tb = os.path.join(A, "tested.bed")
    with open(tb, "w") as f:
        for c, pos in cur.execute("SELECT chrom, start FROM refblock WHERE chrom<=22 UNION SELECT chrom, pos FROM variant WHERE chrom<=22 "
                                  "AND length(ref)=1 AND length(alt)=1 ORDER BY 1, 2"):
            f.write("%d\t%d\t%d\n" % (c, pos, pos))
    extract = "--extract bed1 " + shq(tb)
beds = []
for c in range(1, 23):
    pf = os.path.join(P, "pgen", "chr%d" % c)
    out = os.path.join(A, "panel%s%d" % (tag, c))
    if not os.path.exists(out + ".bed"):
        subprocess.run(f"plink2 --pfile {shq(pf)} --maf 0.05 --snps-only just-acgt {extract} --exclude range {shq(rng)} --indep-pairwise 500kb 1 0.1 "
                       f"--threads {NPROC} --out {shq(out)} > {shq(out + '.log')} 2>&1", shell=True, check=True)
        subprocess.run(f"plink2 --pfile {shq(pf)} --extract {shq(out + '.prune.in')} --make-bed --threads {NPROC} --out {shq(out)} >> {shq(out + '.log')} 2>&1", shell=True, check=True)
    beds.append(out)
fam = [l.split() for l in open(beds[0] + ".fam")]
samples = [x[1] for x in fam]
pop_of = {}
for line in open(os.path.join(P, "founders.pop.tsv")):
    p = line.split()
    pop_of[p[0]] = p[2]
pops = np.array([pop_of.get(s, "?") for s in samples])
N = len(samples)
log("reference samples", N)


def read_bed(prefix, n):
    bim = [l.split() for l in open(prefix + ".bim")]
    m = len(bim)
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)[3:]
    bpv = (n + 3) // 4
    raw = raw.reshape(m, bpv)
    codes = np.zeros((m, bpv * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (raw >> (2 * k)) & 3
    codes = codes[:, :n]
    # plink1 codes: 0 hom A1, 1 missing, 2 het, 3 hom A2  -> A1 dosage
    lut = np.array([2, 255, 1, 0], dtype=np.uint8)
    return bim, lut[codes]


bims, gts = [], []
for b in beds:
    bim, g = read_bed(b, N)
    bims += bim
    gts.append(g)
Gm = np.concatenate(gts, axis=0)          # SNPs x samples, A1 dosage, 255 missing
M = Gm.shape[0]
log("panel SNPs", M)
# the sample at each panel SNP: dosage of A1
chrom_a = np.array([int(b[0]) for b in bims])
pos_a = np.array([int(b[3]) for b in bims])
a1 = [b[4] for b in bims]
a2 = [b[5] for b in bims]
# The variant ids carry the reference allele (chr:pos:REF:ALT), so the
# sample's dosage is counted for A1 whichever allele plink put there.
ref_of = [b[1].split(":")[2] if b[1].count(":") >= 3 else b[5] for b in bims]
ug = np.full(M, -1.0)
for c in range(1, 23):
    idx = np.where(chrom_a == c)[0]
    for i in idx:
        r = G.at(c, int(pos_a[i]))
        a1_is_ref = a1[i] == ref_of[i]
        alt_n = None
        if r[0] == "ref":
            alt_n = 0
        elif r[0] == "var":
            for ref, alt, zyg, ps, dp in r[1]:
                if len(ref) == 1 and len(alt) == 1 and ref == ref_of[i] and alt in (a1[i], a2[i]):
                    alt_n = 2 if zyg == 2 else 1
                    break
        if alt_n is not None:
            ug[i] = float(2 - alt_n if a1_is_ref else alt_n)
log("A1 is the reference allele at %.1f%% of panel SNPs" % (100 * np.mean([a1[i] == ref_of[i] for i in range(M)])))
# sanity: the reference allele must be A2 for hom-ref calls to mean 0
fa = os.path.join(REF, "GRCh38.fa")
called = ug >= 0
log("sample called at %d of %d panel SNPs" % (called.sum(), M))

# ---- reference allele frequencies per population
miss = Gm == 255
Gf = Gm.astype(np.float32)
Gf[miss] = np.nan
upops = sorted(set(pops) - {"?"})
af = {}
for p in upops:
    cols = pops == p
    af[p] = np.nanmean(Gf[:, cols], axis=1) / 2.0
p_all = np.nanmean(Gf, axis=1) / 2.0

# ---- 2. PCA
def pca(colmask, k=10):
    X = Gf[:, colmask]
    p = np.nanmean(X, axis=1) / 2.0
    keep = (p > 0.01) & (p < 0.99)
    X = X[keep]
    p = p[keep]
    sd = np.sqrt(2 * p * (1 - p))
    Z = (X - 2 * p[:, None]) / sd[:, None]
    Z = np.nan_to_num(Z)
    K = (Z.T @ Z) / Z.shape[0]
    w, U = np.linalg.eigh(K)
    order = np.argsort(w)[::-1][:k]
    w, U = w[order], U[:, order]
    scores = U * np.sqrt(w)[None, :] * np.sqrt(Z.shape[0])
    # loadings to project new samples: v = Z U / sqrt(w M)
    V = (Z @ U) / np.sqrt(w)[None, :] / np.sqrt(Z.shape[0])
    u = np.where(called[keep], (ug[keep] - 2 * p) / sd, 0.0)
    proj = (u @ V) / np.sqrt(Z.shape[0]) * np.sqrt(Z.shape[0]) / (called[keep].mean() if called[keep].mean() > 0 else 1)
    return scores, proj, w / w.sum(), K, U, w, Z.shape[0]


def loo_projections(K, U, w, M, k=4):
    """Each reference sample projected onto a PCA fitted WITHOUT it, the same
    way you are projected. A sample's own in-sample score is inflated by its
    self-similarity (projections shrink by about 1 - 1/eigenvalue), so
    anything that compares you with the references (the polygenic score
    ancestry adjustment) must use these, not the in-sample scores."""
    n = K.shape[0]
    out = np.zeros((n, k))
    for i in range(n):
        keep = np.arange(n) != i
        wi, Ui = np.linalg.eigh(K[np.ix_(keep, keep)])
        order = np.argsort(wi)[::-1][:k]
        wi, Ui = wi[order], Ui[:, order]
        # align each axis with the full fit
        sign = np.sign(np.sum(Ui * U[keep][:, :k], axis=0))
        sign[sign == 0] = 1
        Ui = Ui * sign
        out[i] = np.sqrt(M) * (K[i, keep] @ Ui) / np.sqrt(wi)
    return out


cur.executescript("""
DROP TABLE IF EXISTS anc_pop; DROP TABLE IF EXISTS anc_admix; DROP TABLE IF EXISTS anc_pca; DROP TABLE IF EXISTS anc_paint; DROP TABLE IF EXISTS anc_closest;
CREATE TABLE anc_pop(pop TEXT, label TEXT, super TEXT, super_label TEXT, n INTEGER);
CREATE TABLE anc_admix(level TEXT, pop TEXT, label TEXT, frac REAL, lo REAL, hi REAL);
CREATE TABLE anc_pca(sett TEXT, sample TEXT, pop TEXT, super TEXT, pc1 REAL, pc2 REAL, pc3 REAL, pc4 REAL);
CREATE TABLE anc_paint(chrom INTEGER, start INTEGER, end INTEGER, a TEXT, b TEXT, ll REAL);
CREATE TABLE anc_closest(pop TEXT, label TEXT, score REAL, rank INTEGER);
""")
for p in upops:
    lab, sp = LABEL.get(p, (p, "?"))
    cur.execute("INSERT INTO anc_pop VALUES(?,?,?,?,?)", (p, lab, sp, SUPER.get(sp, sp), int((pops == p).sum())))

for sett, mask in [("global", np.ones(N, dtype=bool)), ("eur", np.isin(pops, ["GBR", "CEU", "FIN", "IBS", "TSI"]))]:
    sc, pr, ve, K, U, w, M = pca(mask)
    idx = np.where(mask)[0]
    rows = [(sett, samples[i], pops[i], LABEL.get(pops[i], ("", "?"))[1], *map(float, sc[j, :4])) for j, i in enumerate(idx)]
    rows.append((sett, "YOU", "YOU", "YOU", *map(float, pr[:4])))
    cur.executemany("INSERT INTO anc_pca VALUES(?,?,?,?,?,?,?,?)", rows)
    if sett == "eur":
        # out-of-sample coordinates for the references (set 'eur_oos'): what
        # the polygenic score step regresses on, so references and you sit
        # on the same scale
        lo = loo_projections(K, U, w, M)
        cur.executemany("INSERT INTO anc_pca VALUES(?,?,?,?,?,?,?,?)",
                        [("eur_oos", samples[i], pops[i], LABEL.get(pops[i], ("", "?"))[1], *map(float, lo[j])) for j, i in enumerate(idx)])
        log("European PCA: leave-one-out projections for %d references" % len(idx))
    set_meta(con, "pca_%s_var" % sett, " ".join("%.4f" % v for v in ve[:4]))
    log("PCA %s: you at %.3f %.3f" % (sett, pr[0], pr[1]))

# ---- 3. supervised admixture
refpops = [p for p in upops if p not in ADMIXED]
F = np.stack([af[p] for p in refpops], axis=1)            # M x K
F = np.clip(np.nan_to_num(F, nan=0.5), 0.001, 0.999)
g = ug[called]
Fc = F[called]


def em_step(q, gv, Fv):
    pa = Fv @ q
    pb = (1 - Fv) @ q
    upd = (gv[:, None] * (Fv * q[None, :]) / pa[:, None] + (2 - gv)[:, None] * ((1 - Fv) * q[None, :]) / pb[:, None]).sum(0)
    return upd / (2 * len(gv))


def em(gv, Fv, iters=3000, tol=1e-9):
    """mixture proportions by EM, accelerated with SQUAREM (Varadhan and
    Roland 2008): close populations make the likelihood nearly flat and plain
    EM crawls, stopping far from the optimum"""
    K = Fv.shape[1]
    q = np.full(K, 1.0 / K)
    for _ in range(iters):
        q1 = em_step(q, gv, Fv)
        q2 = em_step(q1, gv, Fv)
        r = q1 - q
        v = q2 - q1 - r
        nv = np.linalg.norm(v)
        if nv == 0:
            q = q2
            break
        a = min(-np.linalg.norm(r) / nv, -1.0)
        qn = q - 2 * a * r + a * a * v
        qn = np.clip(qn, 1e-10, None)
        qn /= qn.sum()
        qn = em_step(qn, gv, Fv)
        if np.abs(qn - q).max() < tol:
            q = qn
            break
        q = qn
    return q


# Each population's frequencies come from about a hundred people, so they
# carry sampling noise, and a free mixture of 21 populations fits that noise
# as small spurious fractions (a few percent of everything). Two guards:
# frequencies shrunk towards the panel mean with a pseudo-count of one
# person, and components under 2% dropped and the rest refitted until none
# is left below the line.
npop = np.array([(pops == p).sum() for p in refpops], dtype=float)
F = np.clip((F * 2 * npop[None, :] + p_all[:, None] * 2) / (2 * npop[None, :] + 2), 0.001, 0.999)
Fc = F[called]
keep = np.arange(len(refpops))
for _ in range(10):
    qk = em(g, Fc[:, keep])
    small = qk < 0.02
    if not small.any():
        break
    keep = keep[~small]
q = np.zeros(len(refpops))
q[keep] = qk
log("populations kept after refitting:", ", ".join(refpops[k] for k in keep))
# bootstrap over 5 Mb blocks
blk = (chrom_a[called] * 1000 + pos_a[called] // 5_000_000)
ub = np.unique(blk)
rs = np.random.default_rng(7)
boots = []
for b in range(30):
    pick = rs.choice(ub, size=len(ub), replace=True)
    sel = np.concatenate([np.where(blk == x)[0] for x in pick])
    bq = np.zeros(len(refpops))
    bq[keep] = em(g[sel], Fc[sel][:, keep])
    boots.append(bq)
boots = np.array(boots)
lo, hi = np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0)
for k, p in enumerate(refpops):
    cur.execute("INSERT INTO anc_admix VALUES(?,?,?,?,?,?)", ("pop", p, LABEL[p][0], float(q[k]), float(lo[k]), float(hi[k])))
sup = {}
for k, p in enumerate(refpops):
    sp = LABEL[p][1]
    sup.setdefault(sp, [0.0, np.zeros(len(boots))])
    sup[sp][0] += q[k]
    sup[sp][1] += boots[:, k]
for sp, (v, bs) in sup.items():
    cur.execute("INSERT INTO anc_admix VALUES(?,?,?,?,?,?)", ("super", sp, SUPER[sp], float(v), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))))
log("admixture:", ", ".join("%s %.1f%%" % (p, q[k] * 100) for k, p in sorted(enumerate(refpops), key=lambda t: -q[t[0]])[:6]))

# ---- 5. closest populations (mean log-likelihood per SNP, single-population model)
scores = []
for p in upops:
    f = np.clip(np.nan_to_num(af[p][called], nan=0.5), 0.001, 0.999)
    ll = np.where(g == 2, 2 * np.log(f), np.where(g == 1, np.log(2 * f * (1 - f)), 2 * np.log(1 - f)))
    scores.append((float(ll.mean()), p))
scores.sort(reverse=True)
best = scores[0][0]
for r, (s, p) in enumerate(scores):
    cur.execute("INSERT INTO anc_closest VALUES(?,?,?,?)", (p, LABEL[p][0], s - best, r + 1))

# ---- 4. painting: pairs of continental groups per window, Viterbi smoothed
groups = ["EUR", "AFR", "EAS", "SAS", "AMR"]
FG = np.stack([np.clip(np.nan_to_num(np.mean([af[p] for p in refpops if LABEL[p][1] == s], axis=0), nan=0.5), 0.005, 0.995) for s in groups], axis=1)
pairs = [(a, b) for a in range(5) for b in range(a, 5)]
W = 250
for c in range(1, 23):
    idx = np.where((chrom_a == c) & called)[0]
    if len(idx) < W:
        continue
    wins = [idx[i:i + W] for i in range(0, len(idx), W)]
    E = np.zeros((len(wins), len(pairs)))
    for wi, ws in enumerate(wins):
        gg = ug[ws]
        for pi, (a, b) in enumerate(pairs):
            fa_, fb_ = FG[ws, a], FG[ws, b]
            p0, p2 = (1 - fa_) * (1 - fb_), fa_ * fb_
            p1 = 1 - p0 - p2
            E[wi, pi] = np.log(np.where(gg == 0, p0, np.where(gg == 1, p1, p2))).sum()
    # Viterbi with a switch cost
    cost = 25.0
    V = E[0].copy()
    back = np.zeros(E.shape, dtype=int)
    for t in range(1, len(wins)):
        best_prev = V.max()
        arg = V.argmax()
        stay = V
        nv = np.where(stay >= best_prev - cost, stay, best_prev - cost)
        back[t] = np.where(stay >= best_prev - cost, np.arange(len(pairs)), arg)
        V = nv + E[t]
    path = [int(V.argmax())]
    for t in range(len(wins) - 1, 0, -1):
        path.append(int(back[t][path[-1]]))
    path.reverse()
    s0 = 0
    for t in range(1, len(wins) + 1):
        if t == len(wins) or path[t] != path[s0]:
            a, b = pairs[path[s0]]
            st = int(pos_a[wins[s0][0]]) if s0 else 1
            en = int(pos_a[wins[t - 1][-1]]) if t < len(wins) else int(chrom_lengths()[str(c)])
            cur.execute("INSERT INTO anc_paint VALUES(?,?,?,?,?,?)", (c, st, en, groups[a], groups[b], float(E[s0:t, path[s0]].sum())))
            s0 = t

top = sorted(sup.items(), key=lambda kv: -kv[1][0])
summary = ", ".join("%.0f%% %s" % (v[0] * 100, SUPER[k]) for k, v in top if v[0] >= 0.005)
near = " and ".join(LABEL[p][0] for s, p in scores[:2])
set_meta(con, "ancestry_summary", "%s; closest to %s" % (summary, near))
set_meta(con, "ancestry_panel_snps", int(called.sum()))
con.commit()
log("summary:", summary, "| closest:", near)
