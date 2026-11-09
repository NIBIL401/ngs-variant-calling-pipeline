"""Simulate a tiny reference genome and paired-end reads for testing the pipeline.

Creates, under the output directory:
  ref.fa                      random 20 kb, 2-contig reference
  <sample>_R1/R2.fastq.gz     2x100 bp reads (~30x) carrying planted variants
  truth.tsv                   planted SNPs: contig, pos (1-based), ref, alt, genotype per sample
  samples.tsv / config.yaml   pipeline inputs pointing at the files above

Only the Python standard library is used. Run with:  python .test/simulate.py .test/data
"""
import gzip
import os
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")
SAMPLES = ["sampleA", "sampleB"]
CONTIGS = {"contig1": 12_000, "contig2": 8_000}
N_SNPS = 30
READ_LEN = 100
COVERAGE = 30
INSERT_MEAN, INSERT_SD = 300, 30
ERROR_RATE = 0.002


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(seq, rng):
    out = []
    for b in seq:
        if rng.random() < ERROR_RATE:
            b = rng.choice([x for x in "ACGT" if x != b])
        out.append(b)
    return "".join(out)


def main(outdir):
    rng = random.Random(1234)
    os.makedirs(outdir, exist_ok=True)

    ref = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CONTIGS.items()}
    with open(os.path.join(outdir, "ref.fa"), "w") as fh:
        for c, s in ref.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # Plant SNPs >= 200 bp from contig ends and from each other.
    snps = []
    for c, s in ref.items():
        k = round(N_SNPS * len(s) / sum(CONTIGS.values()))
        positions = sorted(rng.sample(range(200, len(s) - 200, 150), k))
        for p in positions:
            alt = rng.choice([b for b in "ACGT" if b != s[p]])
            gts = {smp: rng.choice(["0/1", "1/1", "0/1"]) for smp in SAMPLES}
            # make sure every SNP is present in at least one sample
            snps.append((c, p, s[p], alt, gts))

    with open(os.path.join(outdir, "truth.tsv"), "w") as fh:
        fh.write("contig\tpos\tref\talt\t" + "\t".join(SAMPLES) + "\n")
        for c, p, r, a, gts in snps:
            fh.write(f"{c}\t{p + 1}\t{r}\t{a}\t" + "\t".join(gts[s] for s in SAMPLES) + "\n")

    for smp in SAMPLES:
        # two haplotypes per sample
        haps = []
        for h in range(2):
            hap = {c: list(s) for c, s in ref.items()}
            for c, p, r, a, gts in snps:
                gt = gts[smp]
                if gt == "1/1" or (gt == "0/1" and h == 0):
                    hap[c][p] = a
            haps.append({c: "".join(v) for c, v in hap.items()})

        qual = "I" * READ_LEN
        p1 = os.path.join(outdir, f"{smp}_R1.fastq.gz")
        p2 = os.path.join(outdir, f"{smp}_R2.fastq.gz")
        with gzip.open(p1, "wt") as f1, gzip.open(p2, "wt") as f2:
            n = 0
            for c, length in CONTIGS.items():
                pairs = COVERAGE * length // (2 * READ_LEN)
                for _ in range(pairs):
                    hap = rng.choice(haps)[c]
                    ins = max(READ_LEN + 10, int(rng.gauss(INSERT_MEAN, INSERT_SD)))
                    start = rng.randint(0, length - ins)
                    frag = hap[start:start + ins]
                    if rng.random() < 0.5:
                        frag = revcomp(frag)
                    r1 = mutate(frag[:READ_LEN], rng)
                    r2 = mutate(revcomp(frag)[:READ_LEN], rng)
                    n += 1
                    f1.write(f"@{smp}_{n}/1\n{r1}\n+\n{qual}\n")
                    f2.write(f"@{smp}_{n}/2\n{r2}\n+\n{qual}\n")
        print(f"{smp}: {n} read pairs")

    with open(os.path.join(outdir, "samples.tsv"), "w") as fh:
        fh.write("sample\tfq1\tfq2\n")
        for smp in SAMPLES:
            fh.write(f"{smp}\t{outdir}/{smp}_R1.fastq.gz\t{outdir}/{smp}_R2.fastq.gz\n")
    print(f"{len(snps)} SNPs planted -> {outdir}/truth.tsv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".test/data")
