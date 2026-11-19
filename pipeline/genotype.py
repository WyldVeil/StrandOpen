"""Genotype lookup for the import steps: what the sample carries at a
position, from the variant table and the gVCF reference blocks.

    g = Genotyper(con)
    g.at(chrom_idx, pos) -> ("var", [(ref, alt, zyg, pass, dp), ...])
                          | ("ref", dp, passed) | ("nocall",)
    g.alleles(chrom_idx, pos, ref_base) -> list of the sample's alleles at
        a single base (two for diploid, one for X/Y/MT in a male), or None

Reference blocks are loaded per chromosome on first use (sorted starts and
ends in numpy arrays), so a few hundred thousand lookups take seconds."""
import bisect
import numpy as np


class Genotyper:
    def __init__(self, con, haploid=("24", "25")):
        self.con = con
        self.blocks = {}
        self.vars = {}
        self.haploid_chroms = set(int(c) for c in haploid)
        sex = con.execute("SELECT value FROM meta WHERE key='sex'").fetchone()
        self.male = (sex[0] == "male") if sex else True
        # depth means nothing for a chip (a genotype, no reads) or for the reference a
        # differences-only VCF leaves implied (blocks with MinDP 0): every such position counts as called
        cm = con.execute("SELECT value FROM meta WHERE key='call_model'").fetchone()
        self.chip = bool(cm) and cm[0] in ("chip", "vcf")
        if self.male:
            self.haploid_chroms.add(23)

    def _load(self, chrom):
        if chrom in self.blocks:
            return
        rows = self.con.execute("SELECT start, end, dp, pass FROM refblock WHERE chrom=? ORDER BY start", (chrom,)).fetchall()
        a = np.array(rows, dtype=np.int64) if rows else np.zeros((0, 4), dtype=np.int64)
        # running maximum of ends, so overlapping blocks are found
        ends = np.maximum.accumulate(a[:, 1]) if len(a) else a[:, 1]
        self.blocks[chrom] = (a[:, 0], a[:, 1], a[:, 2], a[:, 3], ends)
        vs = {}
        # chip calls set aside as rare (step 5c, FILTER ChipRare) are read as untested
        for pos, ref, alt, zyg, ps, dp in self.con.execute(
                "SELECT pos, ref, alt, zyg, pass, dp FROM variant WHERE chrom=? AND coalesce(filter,'') != 'ChipRare'", (chrom,)):
            vs.setdefault(pos, []).append((ref, alt, zyg, ps, dp))
        self.vars[chrom] = vs

    def at(self, chrom, pos):
        self._load(chrom)
        v = self.vars[chrom].get(pos)
        if v:
            return ("var", v)
        s, e, dp, ps, ends = self.blocks[chrom]
        i = int(np.searchsorted(s, pos, side="right")) - 1
        # walk back over blocks that may still reach pos
        while i >= 0 and ends[i] >= pos:
            if s[i] <= pos <= e[i]:
                return ("ref", int(dp[i]), bool(ps[i]))
            i -= 1
        return ("nocall",)

    def alleles(self, chrom, pos, ref_base, min_dp=1):
        """the sample's alleles at one base; None when not called"""
        hap = chrom in self.haploid_chroms
        r = self.at(chrom, pos)
        if r[0] == "ref":
            if r[1] < min_dp and not self.chip:
                return None
            return [ref_base] if hap else [ref_base, ref_base]
        if r[0] == "var":
            out = []
            for ref, alt, zyg, ps, dp in r[1]:
                if len(ref) != 1 or len(alt) != 1:
                    continue
                if zyg == 3:
                    return [alt]
                if zyg == 2:
                    return [alt] if hap else [alt, alt]
                out.append(alt)
            if not out:
                return None
            if hap:
                return [out[0]]
            if len(out) == 1:
                return [ref_base, out[0]]
            return out[:2]
        return None
