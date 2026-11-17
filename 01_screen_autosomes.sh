#!/usr/bin/env bash
# Purpose: screen autosomes for all N476 samples before cohort QC.
# Phenotype: none. Input: the vcf column of the workspace sample.tsv.
# Assumption: previously inspected single-sample GRCh37 VCFs with GT defined.
# Output: one GT-only VCF and tabix index per table row, outside the code repo.
# Rules: chr1-22; exclude DUPE/TRIALLELIC/ZEROED_OUT_ASSAY; distinct A/C/G/T
# REF and ALT; exclude ALL records at duplicated candidate coordinates.
# Preserve original sample IDs, allele order and GT. Process samples sequentially.

set -euo pipefail

# Part 1. Configuration
TABLE_URI="gs://fc-secure-6c44e649-ddc0-4b8b-a14d-3ae8f9538b09/sample.tsv"
OUT_DIR="/home/jupyter/CHIP_CAD_N476/result/01_screen_autosomes"

if [[ -e "$OUT_DIR" ]]; then
    echo "ERROR: Output directory already exists; nothing overwritten."
    exit 1
fi

echo "Purpose: autosomal screening for all table rows; GRCh37"
echo "Chromosomes: 1-22; distinct single-base A/C/G/T REF and ALT"
echo "Excluded FILTER tags: DUPE,TRIALLELIC,ZEROED_OUT_ASSAY"
echo "Duplicate candidate coordinates: exclude all records within each sample"
echo "Sample fields: original GT only; INFO removed"
echo "Output directory: $OUT_DIR"
bcftools --version

# Part 2. Download the path table
# The result/ parent directory must already exist.
mkdir "$OUT_DIR"
TMP_DIR=$(mktemp -d "$OUT_DIR/tmp.XXXXXX")
trap 'rm -rf -- "$TMP_DIR"' EXIT
gcloud storage cp "$TABLE_URI" "$TMP_DIR/sample.tsv"

# Part 3. Screen each sample, compress and index
python3 - "$TMP_DIR/sample.tsv" "$OUT_DIR" "$TMP_DIR" <<'PY'
import csv
import gzip
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

table_path, output_dir, temporary_dir = sys.argv[1:]
autosomes = {str(i) for i in range(1, 23)}
excluded_filters = {"DUPE", "TRIALLELIC", "ZEROED_OUT_ASSAY"}
bases = {"A", "C", "G", "T"}


def screen_autosomes(input_path, output_path):
    counts = Counter(input=0, excluded_non_autosomal=0,
                     excluded_FILTER=0, excluded_non_SNP=0,
                     excluded_duplicate_coordinate_records=0, retained=0)
    coordinate_counts = Counter()
    candidates = []
    headers = []

    # Step 1. Keep the needed headers and collect candidate SNPs.
    with gzip.open(input_path, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO="):
                continue
            if line.startswith("##FORMAT=") and not line.startswith("##FORMAT=<ID=GT,"):
                continue
            if line.startswith("#"):
                headers.append(line)
                continue

            fields = line.rstrip("\r\n").split("\t")
            counts["input"] += 1
            chrom = fields[0]
            if chrom not in autosomes:
                counts["excluded_non_autosomal"] += 1
                continue
            if excluded_filters.intersection(fields[6].split(";")):
                counts["excluded_FILTER"] += 1
                continue

            ref, alt = fields[3].upper(), fields[4].upper()
            if ref not in bases or alt not in bases or ref == alt:
                counts["excluded_non_SNP"] += 1
                continue

            # Preserve original GT; replace INFO with '.' and FORMAT with 'GT'.
            coordinate = (int(chrom), int(fields[1]))
            gt_index = fields[8].split(":").index("GT")
            gt = fields[9].split(":")[gt_index]
            record = "\t".join(fields[:7] + [".", "GT", gt]) + "\n"
            candidates.append((coordinate, record))
            coordinate_counts[coordinate] += 1

    # Step 2. Write only records whose candidate coordinate occurs once.
    with open(output_path, "x") as handle:
        handle.writelines(headers)
        for coordinate, record in candidates:
            if coordinate_counts[coordinate] == 1:
                handle.write(record)
                counts["retained"] += 1
            else:
                counts["excluded_duplicate_coordinate_records"] += 1
    return counts


with open(table_path, encoding="utf-8-sig", newline="") as handle:
    vcf_paths = [row["vcf"] for row in csv.DictReader(handle, delimiter="\t")]

print(f"Samples_to_process: {len(vcf_paths)}", flush=True)
for row_number, vcf_uri in enumerate(vcf_paths, start=1):
    print(f"Processing_row: {row_number}/{len(vcf_paths)}", flush=True)
    output = Path(output_dir) / f"sample_row_{row_number:04d}.autosomal.screened.vcf.gz"

    # Only one raw VCF and one uncompressed working VCF occupy disk at a time.
    with tempfile.TemporaryDirectory(dir=temporary_dir) as sample_tmp:
        raw_vcf = Path(sample_tmp) / "input.vcf.gz"
        screened_vcf = Path(sample_tmp) / "screened.vcf"
        subprocess.run(["gcloud", "storage", "cp", vcf_uri, raw_vcf], check=True)
        counts = screen_autosomes(raw_vcf, screened_vcf)
        subprocess.run([
            "bcftools", "view", "--no-version", "-Oz",
            "-o", output, screened_vcf
        ], check=True)
        subprocess.run(["bcftools", "index", "--tbi", output], check=True)

    print(f"ROW {row_number} COMPLETE {json.dumps(counts)}", flush=True)

print(f"Samples_processed: {len(vcf_paths)}", flush=True)
PY

# Part 4. Completion
echo "Autosomal_screening: complete"
echo "Output directory: $OUT_DIR"
