#!/usr/bin/env python3
import argparse
import glob
import os
import shlex
import subprocess
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class SampleRow:
    acc: str
    sample: str
    layout: Optional[str] = None  # optional: SE / PE


def run_cmd(cmd: str, cwd: Optional[str] = None, module: Optional[str] = None) -> None:
    """
    Run a shell command. If module is provided, runs via `bash -lc` with `module load`.
    """
    if module:
        wrapped = f"module load {shlex.quote(module)} && {cmd}"
        subprocess.run(["bash", "-lc", wrapped], cwd=cwd, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, shell=True, check=True)


def parse_acc_list(tsv_path: str) -> List[SampleRow]:
    rows: List[SampleRow] = []
    with open(tsv_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # allow extra trailing tabs
            parts = [p for p in parts if p != ""]
            if len(parts) < 2:
                continue
            acc, sample = parts[0].strip(), parts[1].strip()
            layout = parts[2].strip().upper() if len(parts) >= 3 else None
            if layout == "":
                layout = None
            rows.append(SampleRow(acc=acc, sample=sample, layout=layout))
    return rows


def rename_and_compress_fastqs(fastq_dir: str, acc: str, sample: str, threads: int) -> None:
    """
    After fasterq-dump, rename output FASTQs for this accession to sample-based names and pigz compress them.
    Handles:
      - acc.fastq (SE or unpaired/orphan)
      - acc_1.fastq / acc_2.fastq (PE)
    """
    acc_fastq = os.path.join(fastq_dir, f"{acc}.fastq")
    acc_fq1 = os.path.join(fastq_dir, f"{acc}_1.fastq")
    acc_fq2 = os.path.join(fastq_dir, f"{acc}_2.fastq")

    out_files = []

    has_pe = os.path.exists(acc_fq1) or os.path.exists(acc_fq2)

    if has_pe:
        # paired-end outputs
        if os.path.exists(acc_fq1):
            dst1 = os.path.join(fastq_dir, f"{sample}_1.fastq")
            os.replace(acc_fq1, dst1)
            out_files.append(dst1)
        if os.path.exists(acc_fq2):
            dst2 = os.path.join(fastq_dir, f"{sample}_2.fastq")
            os.replace(acc_fq2, dst2)
            out_files.append(dst2)

        # orphan/unpaired reads sometimes go to acc.fastq with --split-3
        if os.path.exists(acc_fastq) and os.path.getsize(acc_fastq) > 0:
            dstu = os.path.join(fastq_dir, f"{sample}.unpaired.fastq")
            os.replace(acc_fastq, dstu)
            out_files.append(dstu)
        else:
            # remove empty orphan file if present
            if os.path.exists(acc_fastq):
                os.remove(acc_fastq)

    else:
        # single-end output: usually acc.fastq, but sometimes acc_1.fastq depending on toolkit/run
        if os.path.exists(acc_fastq):
            dst = os.path.join(fastq_dir, f"{sample}.fastq")
            os.replace(acc_fastq, dst)
            out_files.append(dst)
        elif os.path.exists(acc_fq1) and not os.path.exists(acc_fq2):
            # treat lone _1 as single-end
            dst = os.path.join(fastq_dir, f"{sample}.fastq")
            os.replace(acc_fq1, dst)
            out_files.append(dst)

    if out_files:
        # compress renamed fastqs
        # pigz accepts multiple files; keep simple
        files_quoted = " ".join(shlex.quote(p) for p in out_files)
        run_cmd(f"pigz -p {threads} {files_quoted}", cwd=None)


def main():
    ap = argparse.ArgumentParser(
        description="Download SRA accessions via prefetch + fasterq-dump; auto-handle SE/PE; rename + pigz."
    )
    ap.add_argument("--out-dir", required=True, help="Base output directory (will create sra/, fastq/, tmp/ inside)")
    ap.add_argument("--acc-list", required=True, help="TSV: accession<TAB>sample[<TAB>SE|PE]")
    ap.add_argument("--threads", type=int, default=32, help="Threads for fasterq-dump (-e) and pigz (-p)")
    # ap.add_argument("--module", default=None, help="Optional module name to load (e.g., sratoolkit)")
    args = ap.parse_args()

    out_dir = args.out_dir
    acc_list = args.acc_list
    threads = args.threads
    # module = args.module

    sra_dir = os.path.join(out_dir, "sra")
    fastq_dir = os.path.join(out_dir, "fastq")
    tmp_dir = os.path.join(out_dir, "tmp")

    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    rows = parse_acc_list(acc_list)
    accs = [r.acc for r in rows]

    # 1) prefetch all accessions (like your ACC_ONLY temp file)
    # acc_only_path = os.path.join(tmp_dir, "acc_only.txt")
    # with open(acc_only_path, "w", encoding="utf-8") as f:
    #     for a in accs:
    #         f.write(a + "\n")

    # Run prefetch from within sra_dir (same behavior as your bash script)
    # run_cmd(f"prefetch --option-file {shlex.quote(acc_only_path)}", cwd=sra_dir, module="sratoolkit/3.2.0")

    # os.remove(acc_only_path)

    # 2) fasterq-dump each accession; auto-detect SE/PE by output files created
    for r in rows:
        cmd = (
            f"fastq-dump {shlex.quote(r.acc)} "
            f"--outdir {shlex.quote(fastq_dir)} "
            f"--split-3 "
            f"--origfmt "
        )
        # Run from within sra_dir so fasterq-dump finds the prefetched .sra files
        run_cmd(cmd, cwd=sra_dir, module="sratoolkit/3.2.0")

        rename_and_compress_fastqs(fastq_dir, r.acc, r.sample, threads)

    # Optional: cleanup any leftover accession-named fastqs (shouldn't be any, but safe)
    leftovers = glob.glob(os.path.join(fastq_dir, "SRR*.fastq"))
    for lf in leftovers:
        # if any exist, compress them rather than deleting
        run_cmd(f"pigz -p {threads} {shlex.quote(lf)}", cwd=None)


if __name__ == "__main__":
    main()
