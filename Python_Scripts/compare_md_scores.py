#!/usr/bin/env python3

import pysam
import csv
from collections import defaultdict

BAM1 = "/home/tmhaxs421/brannanlab/tmhaxs421/riboSTAMP_mouse/results/bam_calmd/1-Brain.calmd.bam"
BAM2 = "/home/tmhaxs421/brannanlab/tmhaxs421/riboSTAMP_mouse/results/star_output/1-Brain.Aligned.sortedByCoord.out.bam"
OUT  = "/home/tmhaxs421/brannanlab/tmhaxs421/riboSTAMP_mouse/results/1-brain_md_compare_python.tsv"

def md_map(bam_path):
    d = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            try:
                md = aln.get_tag("MD")
            except KeyError:
                continue
            mate = "1" if aln.is_read1 else ("2" if aln.is_read2 else ".")
            key = f"{aln.query_name}/{mate}"
            # keep the first occurrence if multiple appear
            if key not in d:
                d[key] = md
    return d

md1 = md_map(BAM1)
md2 = md_map(BAM2)

all_keys = sorted(set(md1) | set(md2))

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["read_id", f"md_{BAM1.split('/')[-1]}", f"md_{BAM2.split('/')[-1]}"])
    for k in all_keys:
        w.writerow([k, md1.get(k, "NA"), md2.get(k, "NA")])

print(f"Wrote: {OUT}")