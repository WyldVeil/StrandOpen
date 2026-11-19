#!/usr/bin/env python3
"""Setup: the public reference data Strand Open reads, downloaded once and
shared by every person on this computer.

    setup.py --root DATA_ROOT [--population] [--keep-panel]

Essentials (about 3 GB to download, 6 GB on disk): the GRCh38 genome, gene
models, ClinVar, AlphaMissense, gene constraint and function, conditions,
medicine guidelines (PharmCAT, PharmGKB), the GWAS Catalog, the Y and
mitochondrial trees, ancient DNA annotations, the genome-build chain and the
1000 Genomes accessibility mask. Enough for health, traits, medicines and
lineages.

Population reference (--population, about 33 GB to download): the 1000
Genomes 30x panel (3,202 people) and the polygenic score files, turned into
frequencies, genotypes and the African outgroup (about 12 GB kept; the panel
itself is deleted afterwards unless --keep-panel). Needed for ancestry,
polygenic scores, rarity and Neanderthal DNA.

Downloads resume where they stopped and are checked against the size the
server reports; a file that already exists at its full size is not fetched
again. Progress for the app: <root>/shared/setup/progress.kv. Nothing about
anyone's genome is sent: these are plain downloads of public files.
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
UA = "StrandOpen-setup/1.0"
ENS = "https://ftp.ensembl.org/pub/release-115"
KG1 = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
       "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel%s.vcf.gz")
PGS = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs"
PHYLO = "https://raw.githubusercontent.com/genepi/phylotree-fu-rcrs/main/src/"

# (destination under ref/, url, about)
ESSENTIALS = [
    ("raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", ENS + "/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
     "GRCh38 reference genome (Ensembl 115)"),
    ("raw/Homo_sapiens.GRCh38.115.gff3.gz", ENS + "/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz", "Gene models (Ensembl 115)"),
    ("raw/clinvar.vcf.gz", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz", "ClinVar"),
    ("raw/clinvar.vcf.gz.tbi", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi", "ClinVar index"),
    ("raw/AlphaMissense_hg38.tsv.gz", "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz", "AlphaMissense predictions"),
    ("raw/gnomad.v4.1.constraint_metrics.tsv",
     "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv", "gnomAD gene constraint"),
    ("raw/uniprot_human.tsv.gz", "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=tsv"
     "&fields=accession,gene_primary,protein_name,cc_function,cc_disease&compressed=true", "UniProt gene functions"),
    ("raw/genes_to_disease.txt", "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt",
     "HPO genes and conditions"),
    ("raw/phenotype.hpoa", "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa", "HPO annotations"),
    ("raw/cytoBand.txt.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz", "Chromosome bands (UCSC)"),
    ("raw/chains/hg19ToHg38.over.chain.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
     "GRCh37 to GRCh38 chain (UCSC)"),
    ("raw/1kg_mask/strict.bed", "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/"
     "20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed", "1000 Genomes accessibility mask"),
    ("raw/pharmcat/pharmcat-3.4.0-all.jar", "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat-3.4.0-all.jar", "PharmCAT"),
    ("raw/pharmcat/pharmcat_positions_3.4.0.vcf", "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat_positions_3.4.0.vcf",
     "PharmCAT positions"),
    ("raw/clinicalAnnotations.zip", "https://api.clinpgx.org/v1/download/file/data/clinicalAnnotations.zip", "PharmGKB clinical annotations"),
    ("raw/variants.zip", "https://api.clinpgx.org/v1/download/file/data/variants.zip", "PharmGKB variants"),
    ("raw/gwas-catalog-associations_ontology-annotated-full.zip",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip", "GWAS Catalog"),
    ("raw/yfull_tree.json", "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json", "YFull Y tree"),
    ("raw/snps_hg38.vcf.gz", "http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz", "Y SNP positions (YBrowse)"),
    ("raw/phylotree/tree.xml", PHYLO + "tree.xml", "PhyloTree (mtDNA)"),
    ("raw/phylotree/tree.yaml", PHYLO + "tree.yaml", "PhyloTree"),
    ("raw/phylotree/rules.csv", PHYLO + "rules.csv", "PhyloTree"),
    ("raw/phylotree/weights.txt", PHYLO + "weights.txt", "PhyloTree weights"),
    ("raw/aadr/v66.p1_1240K.aadr.PUB.anno", "https://dataverse.harvard.edu/api/access/datafile/13994515", "Ancient DNA annotations (AADR v66)"),
]
POPULATION = [("raw/20130606_g1k_3202_samples_ped_population.txt",
               "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
               "1000 Genomes sample sheet")]
for _c in [str(i) for i in range(1, 23)] + ["X"]:
    _v = ".v2" if _c == "X" else ""
    POPULATION.append(("1kg/" + os.path.basename(KG1 % (_c, _v)), KG1 % (_c, _v), "1000 Genomes chromosome %s" % _c))
    POPULATION.append(("1kg/" + os.path.basename(KG1 % (_c, _v)) + ".tbi", KG1 % (_c, _v) + ".tbi", "index"))
for _f in ("pgs_all_metadata_scores.csv", "pgs_all_metadata_performance_metrics.csv", "pgs_all_metadata_evaluation_sample_sets.csv"):
    POPULATION.append(("raw/" + _f, PGS + "/metadata/" + _f, "PGS Catalog metadata"))


def pgs_ids():
    ids = []
    for fn in ("pgs_selected.tsv", "pgs_selected_mind.tsv"):
        p = os.path.join(APP, "data", fn)
        if os.path.exists(p):
            for line in open(p):
                x = line.split("\t")[0].strip()
                if x.startswith("PGS") and x not in ids:
                    ids.append(x)
    return ids


class Progress:
    def __init__(self, path):
        self.path = path
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.st = {"state": "running", "phase": "", "item": "", "done_bytes": "0", "total_bytes": "0", "step": "0", "steps": "0",
                   "pid": str(os.getpid()), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
        self.last = 0.0
        self.write()

    def write(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            for k, v in self.st.items():
                f.write("%s=%s\n" % (k, str(v).replace("\n", " ")))
        os.replace(tmp, self.path)

    def set(self, force=True, **kw):
        self.st.update({k: str(v) for k, v in kw.items()})
        now = time.time()
        if force or now - self.last > 0.5:
            self.last = now
            self.write()


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def remote_size(url):
    """the file's size from the server: HEAD, or a one-byte ranged GET for
    servers that refuse HEAD (Dataverse); None when it will not say"""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            if n and n.isdigit() and int(n) > 0:       # a generated file (UniProt) says 0: unknown
                return int(n)
    except (urllib.error.URLError, OSError, ValueError):
        pass
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Range": "bytes=0-0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            cr = r.headers.get("Content-Range", "")
            if r.status == 206 and "/" in cr and cr.rsplit("/", 1)[1].isdigit():
                return int(cr.rsplit("/", 1)[1])
    except (urllib.error.URLError, OSError, ValueError):
        pass
    return None


def fetch(url, dest, prog, done_before, total):
    """download url to dest, resuming a .part file; returns bytes now on disk"""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    size = remote_size(url)
    if os.path.exists(dest) and (size is None or os.path.getsize(dest) == size):
        return os.path.getsize(dest)
    part = dest + ".part"
    for attempt in range(8):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if size is not None and have > size:
            os.remove(part)
            have = 0
        if size is not None and have == size:
            break
        hdr = {"User-Agent": UA}
        if have:
            hdr["Range"] = "bytes=%d-" % have
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=120) as r:
                if have and r.status != 206:      # the server ignored the range: start again
                    have = 0
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        buf = r.read(1 << 20)
                        if not buf:
                            break
                        f.write(buf)
                        have += len(buf)
                        prog.set(force=False, done_bytes=done_before + have)
            if size is None or have == size:
                break
        except (urllib.error.URLError, OSError) as e:
            log("download interrupted (%s), retrying in %d s" % (e, 5 * (attempt + 1)))
            time.sleep(5 * (attempt + 1))
    else:
        raise SystemExit("DOWNLOAD_FAILED: %s" % url)
    if not os.path.exists(part):
        raise SystemExit("DOWNLOAD_FAILED: nothing came back from %s" % url)
    if size is not None and os.path.getsize(part) != size:
        raise SystemExit("DOWNLOAD_INCOMPLETE: %s (%d of %d bytes)" % (url, os.path.getsize(part), size))
    os.replace(part, dest)
    return os.path.getsize(dest)


def run(cmd, env, logf):
    log("run:", cmd if isinstance(cmd, str) else " ".join(cmd))
    r = subprocess.run(cmd, shell=isinstance(cmd, str), env=env, stdout=logf, stderr=subprocess.STDOUT, executable="/bin/bash" if isinstance(cmd, str) else None)
    if r.returncode:
        raise SystemExit("SETUP_STEP_FAILED: %s" % (cmd if isinstance(cmd, str) else cmd[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--population", action="store_true")
    ap.add_argument("--keep-panel", action="store_true")
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    shared = os.path.join(root, "shared")
    ref = os.path.join(shared, "ref")
    kg = os.path.join(shared, "1kg")
    sdir = os.path.join(shared, "setup")
    os.makedirs(sdir, exist_ok=True)
    prog = Progress(os.path.join(sdir, "progress.kv"))
    logf = open(os.path.join(sdir, "log.txt"), "a")
    env = dict(os.environ)
    rb = os.path.join(APP, "runtime", "bin")
    if os.path.isdir(rb):
        env["PATH"] = rb + os.pathsep + env.get("PATH", "")
    # a profile-less environment for the shell steps: they only need the shared paths
    env.update(STRANDOPEN_PROFILE=os.path.join(sdir, "noprofile"), STRANDOPEN_REF=ref, STRANDOPEN_KG=kg)
    items = list(ESSENTIALS)
    if a.population:
        items += POPULATION
        for i in pgs_ids():
            items.append(("raw/pgs/%s_hmPOS_GRCh38.txt.gz" % i, "%s/scores/%s/ScoringFiles/Harmonized/%s_hmPOS_GRCh38.txt.gz" % (PGS, i, i),
                          "Polygenic score %s" % i))
    try:
        prog.set(phase="Checking sizes", steps=len(items))
        sizes = []
        for dest, url, about in items:
            p = os.path.join(ref, dest)
            sizes.append(os.path.getsize(p) if os.path.exists(p) else (remote_size(url) or 0))
        total = sum(sizes)
        log("to fetch: %d files, %.1f GB" % (len(items), total / 1e9))
        done = 0
        prog.set(phase="Downloading", total_bytes=total)
        for k, (dest, url, about) in enumerate(items):
            prog.set(item=about, step=k + 1)
            got = fetch(url, os.path.join(ref, dest), prog, done, total)
            done += got
            prog.set(done_bytes=done)
        # ---- preparing what was downloaded
        prog.set(phase="Preparing the reference genome", item="")
        fa = os.path.join(ref, "GRCh38.fa")
        if not os.path.exists(fa + ".fai"):
            with gzip.open(os.path.join(ref, "raw", "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"), "rb") as src, open(fa + ".tmp", "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 22)
            os.replace(fa + ".tmp", fa)
            run(["samtools", "faidx", fa], env, logf)
        with open(os.path.join(ref, "ready_essentials"), "w") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
        if a.population:
            prog.set(phase="Building the 1000 Genomes reference (about an hour)", item="")
            # Windows has no bash: the same work in Python (identical output, checked on Linux)
            run([sys.executable, "-B", os.path.join(HERE, "kg_build.py")] if os.name == "nt" else ["bash", os.path.join(HERE, "04_1kg.sh")], env, logf)
            prog.set(phase="Building the African outgroup for Neanderthal DNA", item="")
            run([sys.executable, "-B", os.path.join(HERE, "archaic_build.py")] if os.name == "nt" else ["bash", os.path.join(HERE, "setup_archaic.sh")], env, logf)
            if not a.keep_panel:
                for dest, url, about in POPULATION:
                    if dest.startswith("1kg/"):
                        p = os.path.join(ref, dest)
                        if os.path.exists(p):
                            os.remove(p)
            with open(os.path.join(kg, "ready"), "w") as f:
                f.write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
        with open(os.path.join(root, "shared.kv"), "w") as f:
            f.write("ref=%s\nkg=%s\n" % (ref, kg))
        prog.set(state="done", phase="", item="", finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    except SystemExit as e:
        prog.set(state="failed", error=str(e))
        raise
    except Exception as e:                                   # anything unexpected still reaches the app
        prog.set(state="failed", error="%s: %s" % (type(e).__name__, e))
        raise


if __name__ == "__main__":
    main()
