"""Compare pipeline SNP calls against the simulated truth set (recall, precision, genotype concordance).

Usage: python .test/check_truth.py <truth.tsv> <filtered.vcf.gz> [min_recall]
Exits non-zero if recall or precision is below the threshold.
"""
import csv
import gzip
import sys


def load_truth(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = [k for k in rows[0] if k not in ("contig", "pos", "ref", "alt")]
    truth = {(r["contig"], int(r["pos"])): r for r in rows}
    return truth, samples


def load_calls(path):
    calls, samples = {}, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = cols[9:]
                continue
            if cols[6] not in ("PASS", "."):
                continue
            ref, alts = cols[3], cols[4].split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue
            fmt = cols[8].split(":")
            gts = {s: dict(zip(fmt, v.split(":"))).get("GT", "./.").replace("|", "/")
                   for s, v in zip(samples, cols[9:])}
            calls[(cols[0], int(cols[1]))] = {"alt": alts[0], "gts": gts}
    return calls


def main(truth_path, vcf_path, min_score=0.9):
    truth, samples = load_truth(truth_path)
    calls = load_calls(vcf_path)
    tp = [k for k in truth if k in calls and calls[k]["alt"] == truth[k]["alt"]]
    recall = len(tp) / len(truth)
    precision = len(tp) / len(calls) if calls else 0.0

    gt_total = gt_match = 0
    for k in tp:
        for s in samples:
            gt_total += 1
            gt_match += calls[k]["gts"].get(s) == truth[k][s]
    concordance = gt_match / gt_total if gt_total else 0.0

    print(f"truth SNPs:            {len(truth)}")
    print(f"PASS SNP calls:        {len(calls)}")
    print(f"recall:                {recall:.3f}")
    print(f"precision:             {precision:.3f}")
    print(f"genotype concordance:  {concordance:.3f}")
    missed = sorted(set(truth) - set(tp))
    if missed:
        print("missed:", ", ".join(f"{c}:{p}" for c, p in missed))
    if recall < min_score or precision < min_score:
        sys.exit(f"FAIL: recall/precision below {min_score}")
    print("OK")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 0.9)
