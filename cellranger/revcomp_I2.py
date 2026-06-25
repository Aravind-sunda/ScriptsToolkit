#!/usr/bin/env python3
import sys
import pandas as pd
from io import StringIO

# usage : revcomp_I2.py <input.csv> <output.csv>
# usage to copy and paste : 
# mamba activate datascience & python3 /home/tmhaxs421/brannanlab/tmhaxs421/scripts/cellranger/revcomp_I2.py samplesheet.csv samplesheet_revcomp.csv

def revcomp(seq: str) -> str:
    """Reverse-complement for simple A/C/G/T (case-insensitive). Keeps unknowns as-is."""
    comp = {'A':'T','T':'A','C':'G','G':'C',
            'a':'t','t':'a','c':'g','g':'c'}
    return ''.join(comp.get(b, b) for b in seq[::-1])

def main(in_path: str, out_path: str):
    # Read the file, find header starting at line with 'Sample_ID'
    with open(in_path, 'r') as f:
        lines = f.read().splitlines()

    # Keep [Data] if present
    has_data_tag = len(lines) > 0 and lines[0].strip() == "[Data]"

    # Find header line index
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Sample_ID"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find header line starting with 'Sample_ID'.")

    # Build a CSV string from header to end and load with pandas
    csv_str = "\n".join(lines[header_idx:])
    df = pd.read_csv(StringIO(csv_str))

    # Apply reverse complement to 3rd column (Index2). Use column name if present.
    # Prefer 'Index2' if it exists; else use position 2.
    target_col = 'Index2' if 'Index2' in df.columns else df.columns[2]
    df[target_col] = df[target_col].astype(str).map(revcomp)

    # Write back out with [Data] header retained at the top
    with open(out_path, 'w') as out:
        if has_data_tag:
            out.write("[Data]\n")
        df.to_csv(out, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <input.csv> <output.csv>\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
