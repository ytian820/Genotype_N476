#!/usr/bin/env bash
# Purpose: screen the first N476 sample before cohort QC; no phenotype used.
# Input: original single-sample GRCh37 VCF in Google Cloud Storage.
# Output: GT-only autosomal VCF and tabix index, outside the code repository.
# Rules: chr1-22; exclude DUPE/TRIALLELIC/ZEROED_OUT_ASSAY; distinct A/C/G/T
# REF and ALT; exclude ALL records at duplicated candidate coordinates.
# Original sample ID, allele order and GT are preserved. No sample QC yet.

set -euo pipefail

# Part 1. Configuration and preflight checks
VCF_URI="gs://datarepo-79224062-bucket/45b05d1b-4a96-4cb5-8712-0924201d71bc/b11890ac-13de-405c-b5e6-83dca67df2ba/209190550070_R01C01.vcf.gz"
OUT_DIR="/home/jupyter/CHIP_CAD_N476/result/01_screen_first_sample"
OUTPUT="$OUT_DIR/sample_row_1.autosomal.screened.vcf.gz"

if [[ -e "$OUT_DIR" ]]; then
    echo "ERROR: Output directory already exists; nothing overwritten."
    exit 1
fi

echo "Purpose: first-sample autosomal screening; GRCh37"
echo "Chromosomes: 1-22; distinct single-base A/C/G/T REF and ALT"
echo "Excluded FILTER tags: DUPE,TRIALLELIC,ZEROED_OUT_ASSAY"
echo "Duplicate candidate coordinates: exclude all records"
echo "Sample fields: original GT only; INFO removed"
echo "Output directory: $OUT_DIR"
bcftools --version | head -n 2

# Part 2. Download a temporary copy
mkdir "$OUT_DIR"
TMP_DIR=$(mktemp -d "$OUT_DIR/tmp.XXXXXX")
trap 'rm -rf -- "$TMP_DIR"' EXIT

gcloud storage cp "$VCF_URI" "$TMP_DIR/input.vcf.gz"

# Part 3. Screen records and remove duplicate coordinates
python3 - "$TMP_DIR/input.vcf.gz" "$TMP_DIR/screened.vcf" <<'PY'
import gzip
import re
import sys
from collections import Counter

input_path, output_path = sys.argv[1:]
autosomes = {str(i) for i in range(1, 23)}
excluded_filters = {"DUPE", "TRIALLELIC", "ZEROED_OUT_ASSAY"}
bases = {"A", "C", "G", "T"}

counts = Counter()
coordinates = Counter()
candidates = []
headers = []
sample_header = None
gt_definition_found = False

with gzip.open(input_path, "rt") as handle:
    for line in handle:
        if line.startswith("##"):
            if sample_header is not None:
                sys.exit("ERROR: Metadata found after sample header.")
            # Keep GT's definition; omit unused INFO and FORMAT definitions.
            if line.startswith("##FORMAT="):
                if line.startswith("##FORMAT=<ID=GT,"):
                    headers.append(line)
                    gt_definition_found = True
            elif not line.startswith("##INFO="):
                headers.append(line)
            continue

        fields = line.rstrip("\r\n").split("\t")
        if line.startswith("#CHROM\t"):
            fixed_columns = [
                "#CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT"
            ]
            if (sample_header is not None or len(fields) != 10
                    or fields[:9] != fixed_columns or not fields[9]):
                sys.exit("ERROR: Expected one valid single-sample header.")
            sample_header = line
            continue

        if sample_header is None or len(fields) != 10:
            sys.exit("ERROR: Missing header or invalid record column count.")

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

        if not fields[1].isdigit() or int(fields[1]) < 1:
            sys.exit("ERROR: Invalid candidate position.")
        coordinate = (int(chrom), int(fields[1]))

        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        if format_keys.count("GT") != 1:
            sys.exit("ERROR: Candidate lacks a unique GT field.")
        gt_index = format_keys.index("GT")
        if gt_index >= len(sample_values):
            sys.exit("ERROR: Candidate has no GT value.")
        gt = sample_values[gt_index]
        if not re.fullmatch(r"[01.](?:[/|][01.])?", gt):
            sys.exit("ERROR: Unexpected candidate GT encoding.")

        # Keep CHROM through FILTER, replace INFO with '.', and keep only GT.
        record = "\t".join(fields[:7] + [".", "GT", gt]) + "\n"
        candidates.append((coordinate, record))
        coordinates[coordinate] += 1

if sample_header is None or not gt_definition_found:
    sys.exit("ERROR: Required sample header or GT definition missing.")

# Count duplicates AFTER the preceding filters; never arbitrarily keep one.
duplicate_sizes = [n for n in coordinates.values() if n > 1]
counts["candidate"] = len(candidates)
counts["duplicate_coordinates"] = len(duplicate_sizes)
counts["excluded_duplicate_coordinate_records"] = sum(duplicate_sizes)
counts["retained"] = counts["candidate"] - counts["excluded_duplicate_coordinate_records"]

with open(output_path, "x") as handle:
    handle.writelines(headers)
    handle.write(sample_header)
    for coordinate, record in candidates:
        if coordinates[coordinate] > 1:
            continue
        handle.write(record)

for key in [
    "input", "excluded_non_autosomal", "excluded_FILTER", "excluded_non_SNP",
    "candidate", "duplicate_coordinates", "excluded_duplicate_coordinate_records",
    "retained"
]:
    print(f"{key}: {counts[key]}")
PY

# Part 4. Compress, index and report output paths
bcftools view --no-version -Oz -o "$OUTPUT" "$TMP_DIR/screened.vcf"
bcftools index --tbi "$OUTPUT"

ls -lh "$OUTPUT" "$OUTPUT.tbi"
echo "First_sample_screening: complete"
